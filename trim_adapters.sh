#!/bin/bash

# trim_adapters.sh
# This script trims adapter sequences from ONT (Oxford Nanopore Technology) reads using Porechop.
# I've added comments throughout to explain each step.

# Define input and output directories
RAW_DATA_DIR="path/to/raw/data"  # Replace with the path to your raw FASTQ files
OUTPUT_DIR="path/to/output"      # Replace with the path where the trimmed files will be saved

# Define input and output FASTQ files
RAW_FASTQ="${RAW_DATA_DIR}/filtered_reads.fastq"  # Replace with your input file name (filtered data from the previous step)
TRIMMED_FASTQ="${OUTPUT_DIR}/trimmed_reads.fastq"  # Output trimmed FASTQ file

# Check if input file exists
if [ ! -f "$RAW_FASTQ" ]; then
    echo "Error: Filtered FASTQ file not found at $RAW_FASTQ"
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory doesn't exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi
# only once 

# Run Porechop to trim adapters from the reads
echo "Starting adapter trimming of $RAW_FASTQ..."

# I've set the parameters for Porechop. 
# I'm using default parameters for trimming, but you could adjust them depending on the type of adapter sequences used in your experiment.
porechop -i "$RAW_FASTQ" -o "$TRIMMED_FASTQ"

# Check if the trimming was successful
if [ $? -eq 0 ]; then
    echo "Adapter trimming completed successfully. Output saved to $TRIMMED_FASTQ"
else
    echo "Error: Adapter trimming failed!"
    exit 1
fi

# Optionally, I could add further post-processing steps here, like removing low-quality reads, 
# but for now, the adapter-trimmed reads are ready for the next step.

echo "Adapter trimming complete."
