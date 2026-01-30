#!/usr/bin/env bash
# synteny_scaffold_pipeline.sh

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 genome1.fa genome2.fa [minLength] [mapq] [threads] [outdir]"
    exit 1
fi

REF="$1"
QUERY="$2"
MINLEN="${3:-5000}"
MAPQ="${4:-10}"
THREADS="${5:-8}"
OUTDIR="${6:-OUTPUT_DIR}"

mkdir -p "$OUTDIR"

# Checar se ferramentas existem
command -v minimap2 >/dev/null || { echo "minimap2 required"; exit 1; }

RBASE=$(basename "$REF" .fa)
RBASE=$(basename "$RBASE" .fasta)

GBASE=$(basename "$QUERY" .fa)
GBASE=$(basename "$GBASE" .fasta)


echo "Alinhando $QUERY contra $REF (Threads: $THREADS)..."

# 1) Alinhamento Whole-Genome 
# Gerando apenas o PAF para ser mais rápido e economizar espaço em disco
# -x asm20 é bom para divergência de ~5-10%. Se for a mesma espécie, use asm5.
#minimap2 -t "$THREADS" -x asm20 -c --eqx "$REF" "$QUERY" > "$OUTDIR/${GBASE}_vs_ref.paf"
/usr/local/bin/minimap2 -t "$THREADS" -ax asm20 --eqx "$REF" "$QUERY" > "$OUTDIR/${GBASE}_vs_${RBASE}.sam"
/usr/local/bin/minimap2 -t "$THREADS" -x asm20 -c "$REF" "$QUERY" > "$OUTDIR/${GBASE}_vs_${RBASE}.paf"

# 2) Filtragem para Sintenia Real
# Corrigido: Passando variáveis para o awk com -v
echo "Filtrando alinhamentos (MinLen: $MINLEN, MapQ: $MAPQ)..."
awk -v ML="$MINLEN" -v MQ="$MAPQ" '$11 >= ML && $12 >= MQ' "$OUTDIR/${GBASE}_vs_${RBASE}.paf" > "$OUTDIR/${GBASE}_filtered.paf"

# 3) Criar um resumo CSV simples
echo "query_scaffold,query_len,query_start,query_end,strand,target_scaffold,target_len,target_start,target_end,matches,aln_len,mapq" > "$OUTDIR/summary_synteny.csv"

awk 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' "$OUTDIR/${GBASE}_filtered.paf" >> "$OUTDIR/summary_synteny.csv"

echo "Pronto! Resultados em $OUTDIR"
echo "Arquivo de resumo: $OUTDIR/summary_synteny.csv"
