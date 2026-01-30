# Chromosome Copy Number (CCN) Analysis Pipeline

This repository contains a suite of automated scripts designed to analyze and visualize **Chromosome Copy Number (CCN)** variations from Next-Generation Sequencing (NGS) data. The pipeline handles everything from raw BAM processing to the generation of publication-ready heatmaps.

## 🧬 Overview

The pipeline automates a bioinformatics workflow that includes:
* **Genome Coverage Calculation**: Generates RPKM-normalized coverage files (BedGraph) from BAM files using `bamCoverage` from deeptools.
* **Coverage Comparative Analysis**: Computes $log_2$ ratios between specific samples (e.g., Treatment vs. Control) using `bamCompare` from deeptools.
* **Genome Coverage Mapping**: Maps coverage scores to specific genomic coordinates (such as single-copy genes) using `bedtools`.
* **Normalization & Visualization**: Normalizes scores by the genome median to calculate CCN (Chromossome Copy Number) and produces a visual heatmap.

---

## 🛠 Ensure the following tools are installed in your environment:

### Bioinformatics Tools
* **DeepTools**: For `bamCoverage` and `bamCompare`.
* **BedTools**: For `mappingBed` operations.
* **Samtools**: For indexing and genome file generation.
* **SeqKit**: For calculating effective genome size.

### Programming & Environment
* **Python 3**: With `pandas`, `seaborn`, and `matplotlib` libraries.
* **Execution Options**: 
    * **SLURM**: The pipeline is integrated with `sbatch` for high-performance computing clusters.
    * **Direct Execution**: Users can manually execute the generated `.slurm` scripts (which are standard Bash scripts) directly in the terminal if a cluster is not available.

## 📂 Repository Structure

* `runPipeline.sh`: The main script that orchestrates the entire workflow and manages symbolic links.
* `pipeline/`:
    * `bamCoverage.sh`: Script to generate BedGraph coverage files.
    * `bamCompare.sh`: Script to compare two BAM files.
    * `mappingBed.sh`: Script to map BedGraph values to a BED file of interest.
    * `showingDataAsHeatmaps.py`: Python script for CCN calculation and heatmap generation.

---

## 🚀 Getting Started

### 1. Setup
Copy the repository and ensure all `.sh` scripts have execution permissions:
```bash
# Usage
./runPipeline.sh <path_to_bam_directory> <path_to_reference_genome>

# Example
./runPipeline.sh ./data/bams ./genome/reference.fasta
```

---

## Citations
- Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike Dündar, and Thomas Manke. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.
- Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
- Wei Shen*, Botond Sipos, and Liuyang Zhao. 2024. SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta e191. doi:10.1002/imt2.191.
- Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008
