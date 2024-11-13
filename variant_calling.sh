#!/bin/bash

# variant_calling.sh
# This script performs variant calling on the aligned BAM file using Medaka.
# I've added comments to explain each step, including input/output files and parameters.

# Define input and output directories
OUTPUT_DIR="path/to/output"  # Replace with the path where output variant files will be saved
ALIGNED_BAM="${OUTPUT_DIR}/aligned_reads.sorted.bam"  # Replace with your aligned and sorted BAM file
REFERENCE_GENOME="path/to/reference/genome.fasta"  # Replace with the path to your reference genome

# Define output VCF file
OUTPUT_VCF="${OUTPUT_DIR}/variant_calls.vcf"  # Output VCF file with called variants

# Check if input BAM file and reference genome exist
if [ ! -f "$ALIGNED_BAM" ]; then
    echo "Error: Aligned BAM file not found at $ALIGNED_BAM"
    exit 1
fi

if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "Error: Reference genome not found at $REFERENCE_GENOME"
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory doesn't exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi

# Run Medaka for variant calling
echo "Calling variants using Medaka..."

# I've added the command to call variants with Medaka.
# -i: Input BAM file
# -r: Reference genome
# -o: Output directory for Medaka results
# --model: Medaka model to use; 'r941_min_high_g360' is a good choice for high accuracy ONT data
medaka call -i "$ALIGNED_BAM" -r "$REFERENCE_GENOME" -o "$OUTPUT_DIR/medaka_output" --model r941_min_high_g360

# Check if Medaka variant calling was successful
if [ $? -eq 0 ]; then
    echo "Variant calling completed successfully. Output VCF will be saved to $OUTPUT_VCF"
else
    echo "Error: Variant calling with Medaka failed!"
    exit 1
fi

# Convert Medaka results to VCF (if necessary)
echo "Converting Medaka output to VCF format..."

# Medaka outputs a consensus FASTA file by default, but I can use bcftools to convert it to VCF
bcftools mpileup -f "$REFERENCE_GENOME" "$ALIGNED_BAM" | bcftools call -mv -Ob -o "$OUTPUT_VCF"

# Check if conversion to VCF was successful
if [ $? -eq 0 ]; then
    echo "Variant calling and VCF conversion completed successfully. VCF saved to $OUTPUT_VCF"
else
    echo "Error: VCF conversion failed!"
    exit 1
fi

echo "Variant calling complete!"
