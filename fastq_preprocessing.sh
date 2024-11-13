#!/bin/bash

# fastq_preprocessing.sh
# This script is used for preprocessing raw ONT (Oxford Nanopore Technology) FASTQ files.
# It performs quality filtering using NanoFilt and prepares the reads for further analysis.
# I'm using NanoFilt to trim low-quality reads and filter by length.

# Define input and output directories
RAW_DATA_DIR="path/to/raw/data"  # Replace with the path to your raw FASTQ files
OUTPUT_DIR="path/to/output"      # Replace with the path where preprocessed files will be saved

# Define input and output FASTQ files
RAW_FASTQ="${RAW_DATA_DIR}/raw_reads.fastq"  # Replace with your raw input file name
FILTERED_FASTQ="${OUTPUT_DIR}/filtered_reads.fastq"  # Output filtered FASTQ file

# Set quality threshold and minimum read length
QUALITY_THRESHOLD=10  # I'm using a quality threshold of 10; adjust if needed
MIN_LENGTH=100        # I'm setting a minimum length of 100 bases; adjust based on your needs

# Check if input file exists
if [ ! -f "$RAW_FASTQ" ]; then
    echo "Error: Raw FASTQ file not found at $RAW_FASTQ"
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory doesn't exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi

# Preprocess the FASTQ file: quality filtering using NanoFilt
echo "Starting quality filtering of $RAW_FASTQ..."

# I've added NanoFilt here to filter out low-quality reads. 
# I'm setting the quality threshold to $QUALITY_THRESHOLD and the minimum read length to $MIN_LENGTH.
nanoFilt -q $QUALITY_THRESHOLD -l $MIN_LENGTH "$RAW_FASTQ" > "$FILTERED_FASTQ"

# Check if the filtering was successful
if [ $? -eq 0 ]; then
    echo "Filtering completed successfully. Output saved to $FILTERED_FASTQ"
else
    echo "Error: Filtering failed!"
    exit 1
fi

# Optionally, I could add any other preprocessing steps, like adapter trimming or read normalization.
# But for now, the filtered reads are ready for the next step in the pipeline.

echo "Preprocessing complete!"
