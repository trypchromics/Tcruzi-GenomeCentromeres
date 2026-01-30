#!/usr/bin/bash

#It is necessary remove header from BED_FILE with coordinates
#Example: cat input/40peaks-kkt_sorted.bed3 |tail -n +2 > input/40peaks-kkt_sorted_noHeader.bed3

GENOME_FILE=$1
BED_FILE=$2
OUTPUT_FILE=$3

echo -e "Reading the ${1} genome file and ${2} bed file to create ${3} fasta file"

bedtools getfasta -fi $GENOME_FILE -bed $BED_FILE -fo $OUTPUT_FILE
