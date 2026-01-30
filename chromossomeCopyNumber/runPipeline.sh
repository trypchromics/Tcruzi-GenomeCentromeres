#!/usr/bin/env bash

# Create necessary directories
mkdir -p genome input log job final output

# Function to get the BAM directory
get_bam_directory() {
    if [ -n "$1" ]; then
        bam_dir="$1"
    else
        read -p "Please provide the path to the directory containing BAM files: " bam_dir
    fi

    if [ ! -d "$bam_dir" ]; then
        echo "[ERROR] The directory '$bam_dir' does not exist or is not accessible."
        exit 1
    fi

    echo "$bam_dir"
}

# Function to get the reference genome path
get_reference_genome() {
    if [ -n "$2" ]; then
        ref_genome=$(echo "$2" | tr -d '\n' | tr -d '\r')  # Remove newlines and carriage returns
    else
        read -p "Please provide the path to the reference genome: " ref_genome
    fi

    # File verification
    if [ ! -f "$ref_genome" ]; then
        printf "[ERROR] The genome file '%s' was not found or is not accessible.\n" "$ref_genome"
        exit 1
    fi

    echo "$ref_genome"
}

# Function to create symbolic links for BAMs in the input folder
link_bam_files() {
    local bam_dir="$1"
    for bam_file in "$bam_dir"/*.bam*; do
        if [ -e "$bam_file" ]; then
            # Check if the symbolic link already exists
            if [ ! -L "input/$(basename "$bam_file")" ]; then
                ln -s "$bam_file" input/
                echo "Symbolic link created for $bam_file in the input folder"
            else
                echo "Symbolic link for $bam_file already exists, skipping."
            fi
        else
            echo "[ERROR] No BAM files found in $bam_dir"
            exit 1
        fi
    done
}

# Function to create a symbolic link for the reference genome in the genome folder
link_reference_genome() {
    local ref_genome="$1"
    
    # Check if the genome file exists
    if [ -f "$ref_genome" ]; then
        # Create the full path for the symbolic link
        local link_path="genome/$(basename "$ref_genome")"
        
        # Check if the symbolic link already exists
        if [ ! -L "$link_path" ]; then
            ln -s "$ref_genome" genome/
            echo "Symbolic link created for the reference genome in the genome folder."

            # Confirm if the symbolic link was actually created
            if [ -L "$link_path" ]; then
                echo "Confirmation: Symbolic link successfully created at '$link_path'."
            else
                echo "[ERROR] Failed to create symbolic link at '$link_path'."
                exit 1
            fi
        else
            echo "Symbolic link for the reference genome already exists, skipping."
        fi
    else
        echo "[ERROR] The genome file '$ref_genome' was not found or is not accessible."
        exit 1
    fi
}

# Function to submit the job and wait for completion
submit_and_wait() {
    local job_script="$1"
    shift

    if [ ! -f "$job_script" ]; then
        echo "[ERROR] The SLURM job script '$job_script' was not found."
        exit 1
    fi

    # Prompt user for submission options
    read -p "Enter the number of tasks (-n): " num_tasks
    read -p "Enter memory (--mem): " mem

    # Assemble the command
    local sbatch_command="sbatch -n $num_tasks --mem=$mem $job_script"
    echo "The command to be executed is: $sbatch_command"

    # Ask user if they want to proceed
    read -p "Do you want to proceed with job submission? (y/n): " confirm
    if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
        echo "Submission cancelled."
        exit 1  # Exit script only if user cancels
    fi

    echo "Submitting job: $job_script"
    job_id=$(eval $sbatch_command | awk '{print $4}')
    echo "Job ID: $job_id"

    # Verify if Job ID was generated correctly
    if [ -z "$job_id" ]; then
        echo "[ERROR] Failed to submit job. Check the SLURM script: $job_script"
        exit 1
    fi

    # Wait until the job is completed
    echo "Waiting for job $job_id to complete..."
    while squeue -j "$job_id" > /dev/null 2>&1; do
        echo -n "."
        sleep 30
    done
    echo "Job $job_id completed."
}

# Function to create symbolic link for the single copy genes BED file
link_single_copy_genes() {
    read -p "Please provide the path to the single copy genes BED file: " bed_file

    if [ ! -f "$bed_file" ]; then
        echo "[ERROR] The BED file '$bed_file' was not found or is not accessible."
        exit 1
    fi

    # Create the symbolic link in the input folder
    ln -s "$bed_file" input/
    echo "Symbolic link created for the single copy genes BED file in the input folder."
}

# Main function
main() {
    # Get BAM directory and reference genome path
    bam_directory=$(get_bam_directory "$1")
    ref_genome=$(get_reference_genome "$2")
    
    # Display obtained information for debugging
    echo "BAM Directory: $bam_directory"
    echo "Reference Genome Path: $ref_genome"

    # Create symbolic links for BAMs in the input folder
    link_bam_files "$bam_directory"

    wait

    # Create symbolic link for reference genome in the genome folder
    link_reference_genome "$ref_genome"
    
    wait

    # Submit SLURM jobs for bamCoverage and bamCompare
    ./pipeline/bamCoverage.sh input genome/*.fasta
    submit_and_wait "job/bamCoverage.slurm"
    
    ./pipeline/bamCompare.sh input genome/*.fasta
    submit_and_wait "job/bamCompare.slurm"
    
    # Link single copy genes
    link_single_copy_genes # Removed $bed_file arg as the function prompts for it internally

    wait

    # Call the MapBed.sh script
    ./pipeline/MapBed.sh input/*.bed genome/*.fasta
    submit_and_wait "job/mapBed.slurm"

    # Run Python script
    python ./pipeline/ShowingDataAsHeatmaps.py

}

# Execute the main script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$1" "$2"
fi