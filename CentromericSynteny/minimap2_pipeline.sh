#!/usr/bin/env bash
# centro_minimap2_pipeline.sh
# Usage: ./minimap2_pipeline.sh ref.fa centromere.bed target_genome FLANK THREADS OUTDIR

REF_FA="$1"
CENT_BED="$2"
TARGET_GENOME="$3"
FLANK="$4"
THREADS="${5:-8}"
OUTDIR="${6:-OUTPUT_DIR}"
mkdir -p "$OUTDIR"

# Requirements: bedtools, samtools, minimap2, awk, python (pandas optional)
command -v bedtools >/dev/null || { echo "bedtools required"; exit 1; }
command -v minimap2 >/dev/null || { echo "minimap2 required"; exit 1; }
command -v samtools >/dev/null || { echo "samtools required"; exit 1; }

# 1) make flanked BED and fasta
awk -v FL="$FLANK" 'BEGIN{OFS="\t"}{chr=$1; start=$2-FL; if(start<0) start=0; end=$3+FL; name=($4==""?NR:$4); print chr,start,end,name}' "$CENT_BED" > "$OUTDIR/centromere.flank.${FLANK}.bed"

bedtools getfasta -fi "$REF_FA" -bed "$OUTDIR/centromere.flank.${FLANK}.bed" -name -fo "$OUTDIR/centromere_with_flank.fa" 

# 2) loop genomes
echo "peak_id,ref_chr,ref_start,ref_end,other_genome,other_chr,other_start,other_end,alen,match,mapq,qcov" > "$OUTDIR/summary.csv"
#for G in "$GENOMES_DIR"/*.fa "$GENOMES_DIR"/*.fasta; do

G=$TARGET_GENOME

[ -e "$G" ] || continue

GBASE=$(basename "$G")

echo "Processing $GBASE ..."

# map query (centromere) -> target (other genome)
#/usr/local/bin/minimap2 -t "$THREADS" -k19 -w10 -U50,500 --rmq -r1k,100k -g10k -A1 -B4 -O6,26 -E2,1 -s200 -z200 -N50 --eqx "$G" "$OUTDIR/centromere_with_flank.fa" > "$OUTDIR/${GBASE}.sam"
/usr/local/bin/minimap2 -t "$THREADS" -ax asm20 --eqx "$G" "$OUTDIR/centromere_with_flank.fa" > "$OUTDIR/${GBASE}.sam"
/usr/local/bin/minimap2 -t "$THREADS" -x asm20 -c "$G" "$OUTDIR/centromere_with_flank.fa" > "$OUTDIR/${GBASE}.paf"
#grep -v "^@" "$OUTDIR/${GBASE}.sam" > "$OUTDIR/${GBASE}_noHeader.sam"
# filter: require aln_len (column 11) >= 200 and mapq (col12) >= 10 (adjustable)
awk '{
  qlen=$2; qcov=($4-$3)/qlen*100; if($11>=100 && $12>=10 && qcov>=10) print $0
}' "$OUTDIR/${GBASE}.paf" > "$OUTDIR/${GBASE}.filtered.paf"

# format filtered PAF to CSV lines (one line per alignment)
# PAF cols: qname qlen qstart qend strand tname tlen tstart tend matches alen mapq [opt...]
awk -v G="$GBASE" 'BEGIN{FS="\t"; OFS=","}{
  q=$1; qlen=$2; qstart=$3; qend=$4; tname=$6; tstart=$8; tend=$9; matches=$10; alen=$11; mapq=$12;
  qcov= (qend-qstart)/qlen*100;
  print q, tname, tstart, tend, alen, matches, mapq, qcov
}' "$OUTDIR/${GBASE}.filtered.paf" | \

while IFS=, read -r q t tstart tend alen matches mapq qcov; do
  # find ref coords for q in centromere.flank bed (match name)
#  echo "${q},${GBASE},${t},${tstart},${tend},${alen},${matches},${mapq},${qcov}"
  #echo "${qcov} - ok"
  chr_q=$(echo $q|cut -d":" -f3)
  ref_line=$(grep -w -F -m 1 "${chr_q}" "$OUTDIR/centromere.flank.${FLANK}.bed" || true)
  #echo "${ref_line}"
  if [ -z "$ref_line" ]; then
    ref_chr="NA"; ref_s="NA"; ref_e="NA";
  else
    ref_chr=$(echo "$ref_line" | awk '{print $1}'); ref_s=$(echo "$ref_line" | awk '{print $2}'); ref_e=$(echo "$ref_line" | awk '{print $3}');
    echo "$ref_chr,${ref_s},${ref_e}"
  fi
  echo "${q},${ref_chr},${ref_s},${ref_e},${GBASE},${t},${tstart},${tend},${alen},${matches},${mapq},${qcov}" >> "$OUTDIR/summary.csv"
  #echo "${q},${ref_chr},${ref_s},${ref_e},${GBASE},${t},${tstart},${tend},${alen},${matches},${mapq},${qcov}"
done
#done

echo "Done. Results in $OUTDIR/summary.csv and individual .sam and .paf files."
