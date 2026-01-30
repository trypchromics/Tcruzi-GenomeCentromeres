#!/usr/bin/env bash
# centro_minimap2_sampling.sh
# Usage: ./script.sh ref.fa centromere.bed target_genome FLANK N_REPEATS THREADS OUTDIR

REF_FA="$1"
CENT_BED="$2"
TARGET_GENOME="$3"
FLANK="$4"
N_REPEATS="${5:-100}" # Padrão: 100 repetições se não informado
THREADS="${6:-8}"
OUTDIR="${7:-OUTPUT_DIR}"

mkdir -p "$OUTDIR"

# Dependências
command -v bedtools >/dev/null || { echo "bedtools required"; exit 1; }
command -v minimap2 >/dev/null || { echo "minimap2 required"; exit 1; }
command -v samtools >/dev/null || { echo "samtools required"; exit 1; }

echo "--- Iniciando Pipeline de Amostragem ($N_REPEATS iterações) ---"

# 1) Preparar arquivos estáticos (feitos apenas uma vez)
# Região de exclusão (Original + Flank)
awk -v FL="$FLANK" 'BEGIN{OFS="\t"}{chr=$1; start=$2-FL; if(start<0) start=0; end=$3+FL; name=($4==""?NR:$4); print chr,start,end,name}' "$CENT_BED" > "$OUTDIR/centromere.flank.${FLANK}.bed"

# Indexar genoma
samtools faidx "$REF_FA" -o "$OUTDIR/genome.fai"
cut -f1,2 "$OUTDIR/genome.fai" > "$OUTDIR/genome.lengths"

G=$TARGET_GENOME
GBASE=$(basename "$G")
[ -e "$G" ] || { echo "Target genome not found"; exit 1; }

# Cabeçalho do CSV (Adicionei a coluna 'iteration')
echo "iteration,peak_id,ref_chr,ref_start,ref_end,other_genome,other_chr,other_start,other_end,alen,matches,mapq,qcov" > "$OUTDIR/summary_sampling.csv"

# ==========================================
# LOOP DE AMOSTRAGEM
# ==========================================
for i in $(seq 1 "$N_REPEATS"); do
    
    # Gera uma semente aleatória baseada no tempo + PID + contador
    SEED=$(date +%s%N | cut -b1-9)
    # Ou simplesmente $RANDOM se preferir, mas date é mais seguro para loops rápidos
    
    echo "Iteração $i de $N_REPEATS (Seed: $i$SEED)..."

    # A. Bedtools Shuffle
    # -seed $SEED: Garante aleatoriedade controlada
    # -chrom: Mantém mesmo cromossomo
    # -excl: Fica fora da região original
    bedtools shuffle -i "$CENT_BED" \
                     -g "$OUTDIR/genome.lengths" \
                     -excl "$OUTDIR/centromere.flank.${FLANK}.bed" \
                     -maxTries 100000 \
                     -seed "$i$SEED" \
                     -chrom > "$OUTDIR/random_iter_${i}.bed"

    # B. GetFasta (Header será >chr:start-end)
    bedtools getfasta -fi "$REF_FA" -bed "$OUTDIR/random_iter_${i}.bed" -fo "$OUTDIR/random_iter_${i}.fa"

    # C. Minimap2
    # Usamos arquivos temporários com o ID da iteração
    minimap2 -t "$THREADS" -c -x asm20 --eqx "$G" "$OUTDIR/random_iter_${i}.fa" > "$OUTDIR/iter_${i}.paf" 2>/dev/null

    # D. Filtragem e Parsing para o CSV Acumulado
    awk -v G="$GBASE" -v ITER="$i" 'BEGIN{FS="\t"; OFS=","}{
      # Filtros básicos: mapq >= 10, len >= 100 (ajuste se necessário)
      qlen=$2; 
      if (qlen > 0) qcov=($4-$3)/qlen*100; else qcov=0;
      
      if($11>=100 && $12>=10 && qcov>=10) {
          # Parsing das variáveis
          q=$1; tname=$6; tstart=$8; tend=$9; matches=$10; alen=$11; mapq=$12;
          
          # Parsing do nome da query "chr:start-end" para extrair coords
          split(q, parts, ":");
          ref_chr=parts[1];
          split(parts[2], coords, "-");
          ref_s=coords[1];
          ref_e=coords[2];

          # Imprime direto no formato CSV
          print ITER, q, ref_chr, ref_s, ref_e, G, tname, tstart, tend, alen, matches, mapq, qcov
      }
    }' "$OUTDIR/iter_${i}.paf" >> "$OUTDIR/summary_sampling.csv"

    # E. Limpeza (Opcional, mas recomendado para economizar espaço)
    #rm "$OUTDIR/random_iter_${i}.bed" "$OUTDIR/random_iter_${i}.fa" "$OUTDIR/iter_${i}.paf"

done

echo "Amostragem finalizada. Resultados em $OUTDIR/summary_sampling.csv"
