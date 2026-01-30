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

# Creation of working directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$JOB_DIR"

# Calculates the effective genome size
GEN_LEN=$(seqkit stats "$GENOME" | tail -1 | awk '{print $5}' | tr -d ',')

# Creates the script for bamCompare
cat > "$JOB_DIR/bamCompare.slurm" << EOI
#!/usr/bin/env bash

# Define variables
BIN_SIZE=500
SMOOTH_LENGTH=1500
PSEUDOCOUNT=1
MAP_QUALITY=30

# BAM file names
DNAg_1="${INPUT_DIR}/DNAg-DM28c-WT-1_S13.bam"
DNAg_2="${INPUT_DIR}/DNAg-DM28c-WT-2_S14.bam"
OUTPUT_FILE="$OUTPUT_DIR/DNAg_1vs2_bamCompare_log2.bw"

# bamCompare command
bamCompare \\
    --bamfile1 \$DNAg_2 \\
    --bamfile2 \$DNAg_1 \\
    --binSize \$BIN_SIZE \\
    --numberOfProcessors 20 \\
    --scaleFactorsMethod readCount \\
    --operation log2 \\
    --pseudocount \$PSEUDOCOUNT \\
    --exactScaling \\
    --smoothLength \$SMOOTH_LENGTH \\
    --extendReads \\
    --ignoreDuplicates \\
    --centerReads \\
    --minMappingQuality \$MAP_QUALITY \\
    --skipNAs \\
    --skipZeroOverZero \\
    --outFileFormat bigwig \\
    --outFileFormat bigwig \\
    -o \$OUTPUT_FILE > $LOG_DIR/DNAg_1vs2_bamCompare_log2.out 2> $LOG_DIR/DNAg_1vs2_bamCompare_log2.err

wait

exit 0
EOI

chmod +x "$JOB_DIR/bamCompare.slurm"
exit 0