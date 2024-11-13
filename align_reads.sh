#!/bin/bash

# align_reads.sh
# This script aligns ONT reads to a reference genome using minimap2.

# Define input and output directories
RAW_DATA_DIR="path/to/raw/data"   # Replace with the path to your preprocessed (filtered/trimmed) FASTQ files
OUTPUT_DIR="path/to/output"       # Replace with the path where BAM files and logs will be saved

# Define input FASTQ file and reference genome
FASTQ_FILE="${RAW_DATA_DIR}/filtered_reads.fastq"  # Replace with your filtered or adapter-trimmed FASTQ file
REFERENCE_GENOME="path/to/reference/genome.fasta"  # Replace with the path to your reference genome

# Define output BAM file
OUTPUT_BAM="${OUTPUT_DIR}/aligned_reads.bam"  # Output BAM file for aligned reads

# Check if input FASTQ file and reference genome exist
if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file not found at $FASTQ_FILE"
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

# Run minimap2 to align reads to the reference genome
echo "Aligning reads from $FASTQ_FILE to $REFERENCE_GENOME..."

# I've added the parameters for minimap2:
# -ax map-ont: ONT read mapping mode
# --secondary=no: Disables secondary alignments (e.g., secondary hits)
minimap2 -ax map-ont "$REFERENCE_GENOME" "$FASTQ_FILE" | samtools view -Sb - > "$OUTPUT_BAM"

# Check if minimap2 and samtools ran successfully
if [ $? -eq 0 ]; then
    echo "Alignment completed successfully. Output BAM saved to $OUTPUT_BAM"
else
    echo "Error: Alignment failed!"
    exit 1
fi

# Sort and index the BAM file using samtools
echo "Sorting and indexing the BAM file..."

samtools sort "$OUTPUT_BAM" -o "${OUTPUT_BAM%.bam}.sorted.bam"
samtools index "${OUTPUT_BAM%.bam}.sorted.bam"

# Check if sorting and indexing were successful
if [ $? -eq 0 ]; then
    echo "BAM file sorted and indexed successfully."
else
    echo "Error: Sorting and indexing failed!"
    exit 1
fi

echo "Alignment complete!"
