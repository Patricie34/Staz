#!/bin/bash

set -e

# Input/output directories
RESULTS_DIR="/path/to/results/"
REFERENCE_DIR="/path/to/reference/"
REF_SEQ="Homo_sapiens_assembly38.fasta"
REF_ANNOT="1000G_phase1.snps.high_confidence.hg38.vcf.gz"

# Mark Duplicates
input_file_MD="$RESULTS_DIR/bam/input_MD.bam"

# Check if input file exists
if test -e "$input_file_MD"; then
    # Run GATK MarkDuplicates
    java -jar gatk-package-4.2.3.0-local.jar MarkDuplicates \
        -I "$input_file_MD" \
        -O "$RESULTS_DIR/bam/output_MD.bam" \
        -M "$RESULTS_DIR/bam/marked_dup_metrics.txt"
else
    # Print error message
    echo "Error: Input file for MarkDuplicates does not exist."
    # Exit with error code
    exit 1
fi

# Base Recalibrator
input_file_BR="$RESULTS_DIR/bam/output_MD.bam"

# Check if input file exists
if test -e "$input_file_BR"; then
    # Run Base Recalibrator
    java -jar gatk-package-4.2.3.0-local.jar BaseRecalibrator \
       -I "$input_file_BR" \
       -R "$REFERENCE_DIR/seq/$REF_SEQ" \
       --known-sites "$REFERENCE_DIR/annot/$REF_ANNOT" \
       -O "$RESULTS_DIR/bam/output_BR.table"
else
    echo "Error: Input file for Base Recalibrator does not exist."
    exit 1
fi

# ApplyBQSR
input_file_AB="$RESULTS_DIR/bam/output_MD.bam"
input_recal_file="$RESULTS_DIR/bam/output_BR.table"

# Check if input files exist
if test -e "$input_file_AB" && test -e "$input_recal_file"; then
    # Run ApplyBQSR
    java -jar gatk-package-4.2.3.0-local.jar ApplyBQSR \
       -R "$REFERENCE_DIR/seq/$REF_SEQ" \
       -I "$input_file_AB" \
       --bqsr-recal-file "$input_recal_file" \
       -O "$RESULTS_DIR/bam/output_AB.bam"
else
    echo "Error: One or both input files do not exist."
    exit 1
fi

# HaplotypeCaller
input_file_HC="$RESULTS_DIR/bam/output_AB.bam"

# Check if input files exist
if test -e "$REFERENCE_DIR/seq/$REF_SEQ" && test -e "$RESULTS_DIR/bam/output_AB.bam"; then
    # Run HaplotypeCaller
    java -Xmx4g -jar gatk-package-4.2.3.0-local.jar HaplotypeCaller \
        -R "$REFERENCE_DIR/seq/$REF_SEQ" \
        -I "$RESULTS_DIR/bam/output_AB.bam" \
        -O "$RESULTS_DIR/bam/output_HC.g.vcf.gz" \
        -ERC GVCF
else
    echo "Error: One or both input files do not exist."
    exit 1
fi
