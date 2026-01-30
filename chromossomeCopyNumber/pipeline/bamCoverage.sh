#!/usr/bin/env bash

# Checks if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <bam_directory> <path_to_reference_genome>"
    exit 1
fi

# Input and output directories
INPUT_DIR="$1"
GENOME="$2"
OUTPUT_DIR="output"
LOG_DIR="log"
JOB_DIR="job"

#echo "$INPUT_DIR and $GENOME and $LOG_DIR"

# Calculates the effective genome size
GEN_LEN=$(seqkit stats "$GENOME" | tail -1 | awk '{print $5}' | tr -d ',')

echo "BAM files found in $INPUT_DIR:"
ls ${INPUT_DIR}/*.bam

# Creates the script for bamCoverage
cat > "$JOB_DIR/bamCoverage.slurm" << EOI
#!/usr/bin/env bash

# Define variables
BIN_SIZE=500  # Bin size in bases
MAP_QUALITY=30  # Minimum mapping quality

# Input BAM files
for BAM_FILE in ${INPUT_DIR}/*.bam; do
  OUT_FILE=${OUTPUT_DIR}/\$(basename \${BAM_FILE%.*}).bedgraph

  # bamCoverage command
  bamCoverage \\
      -p 20 \\
      --bam \$BAM_FILE \\
      --binSize \$BIN_SIZE \\
      --effectiveGenomeSize ${GEN_LEN} \\
      --normalizeUsing RPKM \\
      --outFileFormat bedgraph \\
      --exactScaling \\
      --ignoreDuplicates \\
      --extendReads \\
      --centerReads \\
      --minMappingQuality \$MAP_QUALITY \\
      --skipNAs \\
      --outFileName \$OUT_FILE \\
      > $LOG_DIR/\$(basename \${BAM_FILE%.*})_bamCoverage.out 2> $LOG_DIR/\$(basename \${BAM_FILE%.*})_bamCoverage.err

  echo "Chromosome coverage completed for \$(basename \${BAM_FILE}). Output: \$OUT_FILE"
done

exit 0
EOI
chmod +x "$JOB_DIR/bamCoverage.slurm"

exit 0