#!/usr/bin/env bash
# analyze_regions_normalized.sh
#
# Objetivo: 
#   1. Gera BAM Mestre.
#   2. Itera sobre BED.
#   3. Gera gráficos normalizados (0-100%) para facilitar comparação visual.
#
# Uso:
#   ./analyze_regions_normalized.sh -b input.bam -R regioes.bed -f assembly.fasta -o outdir

set -euo pipefail

# --- Configurações Padrão ---
CUTOFF_CLIP=200
CUTOFF_MAPQ=30
THREADS=4
BED_FILE=""
OUTDIR="analysis_out"
INPUT_BAM=""
GENOME_FASTA=""
MISMATCHES=0.02

# --- Parse Arguments ---
while getopts "b:R:f:o:c:q:t:m:h" opt; do
  case ${opt} in
    b) INPUT_BAM="$OPTARG" ;;
    R) BED_FILE="$OPTARG" ;;
    f) GENOME_FASTA="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    c) CUTOFF_CLIP="$OPTARG" ;;
    q) CUTOFF_MAPQ="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MISMATCHES="$OPTARG" ;;
    h) echo "Usage: $0 -b input.bam -R regions.bed ..."; exit 0 ;;
    *) echo "Invalid option"; exit 1 ;;
  esac
done

if [[ -z "$INPUT_BAM" || -z "$BED_FILE" || -z "$GENOME_FASTA" ]]; then
  echo "ERROr: Missing arguments. Use -b BAM, -R BED e -f FASTA."
  exit 1
fi

mkdir -p "$OUTDIR"
MAIN_LOG="$OUTDIR/pipeline_global.log"
echo "Pipeline started: $(date)" > "$MAIN_LOG"

if ! python3 -c "import pysam, matplotlib, numpy, Bio" 2>/dev/null; then
    echo "CRITIC ERROR: Python libraries missing." | tee -a "$MAIN_LOG"
    exit 1
fi

# ==============================================================================
#  1. PREPARAÇÃO DO SCRIPT PYTHON (NORMALIZADO %)
# ==============================================================================
PYTHON_SCRIPT="$OUTDIR/classify_quality.py"
cat << 'EOF' > "$PYTHON_SCRIPT"
import argparse, os, pysam, numpy as np, matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def analyze_region(bam_path, region_str, fasta_path):
    try:
        if ':' in region_str and '-' in region_str:
            chrom, coords = region_str.split(':')
            start, end = map(int, coords.split('-'))
        else: return None, "Invalid Format"
    except: return None, "Parse Error"

    samfile = pysam.AlignmentFile(bam_path, "rb")
    fasta_ref = pysam.FastaFile(fasta_path)
    region_len = end - start
    
    # Coverage Array
    coverage_array = np.zeros(region_len)
    
    for col in samfile.pileup(chrom, start, end, truncate=True):
        rel_pos = col.pos - start
        if 0 <= rel_pos < region_len:
            coverage_array[rel_pos] = col.n

    # Read Stats
    total_reads = 0; high_clip = 0; total_mm = 0; total_aln = 0
    for read in samfile.fetch(chrom, start, end):
        total_reads += 1
        if read.query_length > 0:
            soft = read.get_cigar_stats()[0][4]
            if (soft / read.query_length) > 0.20: high_clip += 1
        if read.has_tag("NM"):
            total_mm += read.get_tag("NM")
            total_aln += read.query_alignment_length
    samfile.close()

    # Classification Stats
    zeros = np.count_nonzero(coverage_array == 0)
    mean_cov = np.mean(coverage_array)
    max_cov = np.max(coverage_array) if len(coverage_array) > 0 else 0
    
    pct_clip = (high_clip/total_reads*100) if total_reads > 0 else 0
    err_rate = (total_mm/total_aln*100) if total_aln > 0 else 0

    # Logic
    status = "good"; reasons = []
    if zeros > 0: status="bad"; reasons.append(f"Zero cov ({zeros}bp)")
    if pct_clip > 20.0: status="bad"; reasons.append(f"High clip ({pct_clip:.1f}%)")
    if err_rate > 5.0: status="bad"; reasons.append(f"High mismatch ({err_rate:.1f}%)")
    
    if status == "good":
        if mean_cov < 5.0: status="doubtful"; reasons.append(f"Low cov ({mean_cov:.1f}x)")
        elif pct_clip > 10.0: status="doubtful"; reasons.append(f"Mod. clip ({pct_clip:.1f}%)")

    if not reasons: reasons.append("Clean")
    
    try: seq_str = fasta_ref.fetch(chrom, start, end)
    except: seq_str = "N" * region_len

    return {
        "status": status, "reasons": "; ".join(reasons),
        "coverage": coverage_array, "max_cov": max_cov, 
        "sequence": seq_str, "region_name": region_str.replace(":", "_")
    }, None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", required=True); parser.add_argument("-r", required=True)
    parser.add_argument("-f", required=True); parser.add_argument("-o", required=True)
    args = parser.parse_args()

    res, err = analyze_region(args.b, args.r, args.f)
    if err: print(f"PY_ERR:{err}"); return

    # Save Fasta
    fa_path = os.path.join(args.o, f"{res['region_name']}_{res['status']}.fasta")
    with open(fa_path, "w") as f:
        SeqIO.write(SeqRecord(Seq(res['sequence']), id=res['region_name'], description=res['reasons']), f, "fasta")

    # --- PLOT NORMALIZADO (PORCENTAGEM) ---
    plt.figure(figsize=(10,4))
    
    max_val = res['max_cov']
    y_data = res['coverage']
    
    # Se houver dados, normaliza para porcentagem (0-100%)
    if max_val > 0:
        y_norm = (y_data / max_val) * 100
        ylabel_txt = "Relative Coverage (%)"
        
        # Calcula onde fica a linha de 5x na escala de porcentagem
        thresh_5x_pct = (5 / max_val) * 100
    else:
        y_norm = y_data
        ylabel_txt = "Coverage (Counts)"
        thresh_5x_pct = 0

    # Plotar área preenchida
    plt.fill_between(range(len(y_norm)), y_norm, color='#1f77b4', alpha=0.6, label='Depth')
    plt.plot(y_norm, color='#1f77b4', alpha=0.9) # Linha de borda
    
    # Títulos e Eixos
    #plt.title(f"{res['region_name']} [{res['status'].upper()}]\nMax Depth: {int(max_val)}x | {res['reasons']}")
    plt.title(f"{res['region_name']}")
    plt.ylabel(ylabel_txt)
    plt.xlabel("Position (bp)")
    
    # Configura Eixo Y fixo 0-110%
    plt.ylim(0, 110)
    
    # Linha de 5x (vermelha)
    if max_val > 0 and thresh_5x_pct <= 100:
        plt.axhline(y=thresh_5x_pct, color='red', linestyle='--', alpha=0.8, label=f'5x Threshold ({thresh_5x_pct:.1f}%)')
        plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(os.path.join(args.o, "coverage_plot_norm.png"))
    plt.close()

    print(f"CLASS:{res['status'].upper()}")
    print(f"REASON:{res['reasons']}")

if __name__=="__main__": main()
EOF


# ==============================================================================
#  2. PREPARAÇÃO DO BAM GLOBAL
# ==============================================================================
BASENAME=$(basename "${INPUT_BAM%.*}")
MASTER_BAM="$OUTDIR/${BASENAME}.GLOBAL_MISMATCH_FILTERED.bam"

echo "--- Etapa Global: Preparando BAM Mestre ---" | tee -a "$MAIN_LOG"

if [[ -f "$MASTER_BAM" && -f "${MASTER_BAM}.bai" ]]; then
    echo "BAM Mestre já existe." | tee -a "$MAIN_LOG"
else
    echo "Gerando BAM Mestre (Filtro NM)..." | tee -a "$MAIN_LOG"
    if ! samtools quickcheck "$INPUT_BAM"; then echo "Input BAM inválido!"; exit 1; fi
    
    samtools view -h "$INPUT_BAM" | awk -v mismatch_limit="$MISMATCHES" '{
        if ($1 ~ /^@/) { print $0; next }
        nm_tag = 0;
        for (i=12; i<=NF; i++) {
            if ($i ~ /^NM:i:/) { split($i, a, ":"); nm_tag = a[3]; break; }
        }
        read_len = length($10);
        if (read_len > 0 && (nm_tag / read_len) <= mismatch_limit) { print $0; }
    }' | samtools sort -@ "$THREADS" -o "$MASTER_BAM" -
    samtools index "$MASTER_BAM"
fi


# ==============================================================================
#  3. LOOP SOBRE O ARQUIVO BED
# ==============================================================================
echo "--- Iniciando processamento das regiões ---" | tee -a "$MAIN_LOG"

TOTAL_REGIONS=$(grep -c "^[^#]" "$BED_FILE" || true)
COUNT=0

grep "^[^#]" "$BED_FILE" | tr -d '\r' | while read -r chrom start end rest; do
    
    if [[ -z "$chrom" || -z "$start" || -z "$end" ]]; then continue; fi
    COUNT=$((COUNT+1))
    
    chrom=$(echo "$chrom" | xargs)
    start=$(echo "$start" | xargs)
    end=$(echo "$end" | xargs)

    REGION_STR="${chrom}:${start}-${end}"
    REGION_ID="${chrom}_${start}_${end}"
    
    if [[ -n "${rest:-}" ]]; then
        REGION_NAME=$(echo "$rest" | awk '{print $1}') 
        SUBDIR="$OUTDIR/${REGION_NAME}_${REGION_ID}"
    else
        SUBDIR="$OUTDIR/${REGION_ID}"
    fi
    
    mkdir -p "$SUBDIR"
    LOG="$SUBDIR/region.log"
    
    echo "[$COUNT/$TOTAL_REGIONS] $REGION_STR -> $SUBDIR" | tee -a "$MAIN_LOG"
    
    # 3.1 Depth Raw e Stats
    DEPTH_RAW="$SUBDIR/depth.raw.txt"
    if ! samtools depth -r "$REGION_STR" "$MASTER_BAM" > "$DEPTH_RAW" 2>>"$LOG"; then
        echo "ERRO: Falha depth $REGION_STR" >> "$LOG"; continue
    fi
    awk '{sum+=$3; n++} END{if(n>0) print "Mean_Cov: "sum/n; else print "Mean_Cov: 0"}' "$DEPTH_RAW" >> "$LOG"

    # 3.2 Soft-Clips
    CLIPPED_LIST="$SUBDIR/clipped_reads.txt"
    samtools view "$MASTER_BAM" "$REGION_STR" | python3 -c '
import sys, re
reg = re.compile(r"(\d+)S")
limit = '$CUTOFF_CLIP'
for line in sys.stdin:
    p = line.split("\t")
    if len(p)<6: continue
    s = sum(int(m.group(1)) for m in reg.finditer(p[5]))
    if s > limit: print(p[0])
' | sort -u > "$CLIPPED_LIST"
    
    # 3.3 Final Filtered BAM
    FINAL_BAM="$SUBDIR/${REGION_ID}.final_filtered.bam"
    samtools view -h "$MASTER_BAM" "$REGION_STR" | python3 -c '
import sys, re
try: max_c=int(sys.argv[1]); min_q=int(sys.argv[2])
except: sys.exit(1)
pat = re.compile(r"(\d+)[SH]")
for line in sys.stdin:
    if line.startswith("@"): sys.stdout.write(line); continue
    p = line.split("\t")
    if len(p)<6: continue
    try:
        q = int(p[4])
        c = sum(int(m.group(1)) for m in pat.finditer(p[5]))
        if c <= max_c and q >= min_q: sys.stdout.write(line)
    except: continue
' "$CUTOFF_CLIP" "$CUTOFF_MAPQ" | samtools sort -o "$FINAL_BAM" -
    samtools index "$FINAL_BAM"
    
    # 3.4 BigWig (Opcional)
    if command -v bamCoverage >/dev/null 2>&1; then
        if [[ $(samtools view -c "$FINAL_BAM") -gt 0 ]]; then
            bamCoverage -b "$FINAL_BAM" -o "$SUBDIR/coverage.filtered.bw" --binSize 10 2>/dev/null
        fi
    fi

    # 3.5 Classificação e Plot Normalizado
    # Nota: Removemos o --ymax pois agora é porcentagem fixa
    python3 "$PYTHON_SCRIPT" \
        -b "$MASTER_BAM" \
        -r "$REGION_STR" \
        -f "$GENOME_FASTA" \
        -o "$SUBDIR" >> "$LOG" 2>&1
    
    CLASS=$(grep "CLASS:" "$LOG" | cut -d: -f2 || echo "N/A")
    
    # 3.6 IGV
    cat > "$SUBDIR/igv_session.txt" <<EOF
genome $GENOME_FASTA
load $MASTER_BAM
load $FINAL_BAM
goto $REGION_STR
collapse
snapshotDirectory $SUBDIR
snapshot region_snapshot.png
exit
EOF

    echo "Status: $CLASS" | tee -a "$MAIN_LOG"
done

echo "Completed." | tee -a "$MAIN_LOG"
