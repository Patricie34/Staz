#!/bin/bash

# Input/output directories
RAW_DATA_DIR="/path/to/raw_data/"
RESULTS_DIR="/path/to/results/"
REFERENCE="/path/to/reference/Homo_sapiens_assembly38.fasta"

# File names
INPUT_FASTQ="input.fastq.gz"
OUTPUT_SAM="output.sam"
OUTPUT_BAM="output.bam"
OUTPUT_SORTED_BAM="output_sorted.bam"
OUTPUT_SORTED_MAPPED_BAM="output_sorted_mapped.bam"

# Run FASTQC
fastqc -o "$RESULTS_DIR/fastqc" "$RAW_DATA_DIR$INPUT_FASTQ" &
fastqc --outdir "$RESULTS_DIR/fastqc" "$RAW_DATA_DIR$INPUT_FASTQ" &
set -e

# Run BWA
bwa mem -R "@RG\\tID:HG002_HiSeq\\tSM:HG002_HiSeq\\tPL:illumina" "$REFERENCE" \
"$RAW_DATA_DIR$INPUT_FASTQ" \
"$RAW_DATA_DIR$INPUT_FASTQ" \
> "$RESULTS_DIR/bam/$OUTPUT_SAM"

# Convert .sam file to .bam
samtools view -bS "$RESULTS_DIR/bam/$OUTPUT_SAM" \
> "$RESULTS_DIR/bam/$OUTPUT_BAM"

# Further processing on .bam file
rm "$RESULTS_DIR/bam/$OUTPUT_SAM"

samtools sort "$RESULTS_DIR/bam/$OUTPUT_BAM" \
-o "$RESULTS_DIR/bam/$OUTPUT_SORTED_BAM"

# Extract only mapped reads
samtools view -bS -F 4 "$RESULTS_DIR/bam/$OUTPUT_SORTED_BAM" \
> "$RESULTS_DIR/bam/$OUTPUT_SORTED_MAPPED_BAM"

# Indexing the sorted and mapped .bam file
samtools index "$RESULTS_DIR/bam/$OUTPUT_SORTED_MAPPED_BAM"
