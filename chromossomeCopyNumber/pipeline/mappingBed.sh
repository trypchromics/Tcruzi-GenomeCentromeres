#!/usr/bin/env bash

# Checks if the number of arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <path_to_single_copy_genes_BED> <path_to_reference_genome>"
    exit 1
fi

# Directories and files
SINGLE_COPY_GENES_BED="$1"
GENOME_DIR="genome"
GENOME_FILE="$2"
OUTPUT_DIR="output"
LOG_DIR="log"
JOB_DIR="job"

# Generates the .fai and .genome files
samtools faidx "$GENOME_FILE"
cut -f1,2 "${GENOME_FILE}.fai" > "${GENOME_DIR}/tryCru-Dm28c-lcc2024.genome"

# Creates the SLURM script for mapBed with the generated windows
cat > "$JOB_DIR/mapBed.slurm" << EOI
#!/usr/bin/env bash

# Single copy genes BED file
BED_FILE="${SINGLE_COPY_GENES_BED}"

# Previously generated genome file
GENOME_FILE="${GENOME_DIR}/tryCru-Dm28c-lcc2024.genome"

# Process each BEDGRAPH file in the output folder
for BEDGRAPH in ${OUTPUT_DIR}/*.bedgraph; do
  BASENAME=\$(basename \${BEDGRAPH%.*})
  bedtools map -a \${BED_FILE} -b \${BEDGRAPH} -c 4 -o mean -g \${GENOME_FILE} > ${OUTPUT_DIR}/\${BASENAME}_mapbed_output.bedgraph
done

exit 0
EOI

# Makes the mapBed.slurm script executable
chmod +x "$JOB_DIR/mapBed.slurm"

exit 0