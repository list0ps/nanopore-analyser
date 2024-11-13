#!/bin/bash

# qc_analysis.sh
# This script performs quality control (QC) on the ONT reads using NanoPlot and FastQC.
# I use NanoPlot to visualize the read length distributions and FastQC to generate detailed quality reports.
# This helps ensure that the reads are of high quality before downstream analysis.

# Define input and output directories
RAW_DATA_DIR="path/to/raw/data"   # Replace with the path to your raw data (filtered or adapter-trimmed FASTQ files)
OUTPUT_DIR="path/to/output"       # Replace with the path where QC reports and plots will be saved

# Define input FASTQ files (filtered or adapter-trimmed reads)
FASTQ_FILE="${RAW_DATA_DIR}/filtered_reads.fastq"  # Replace with your input file name
FASTQC_DIR="${OUTPUT_DIR}/fastqc_reports"  # Directory to save FastQC reports
PLOTTING_DIR="${OUTPUT_DIR}/nanoplot"  # Directory to save NanoPlot visualizations

# Check if the input file exists
if [ ! -f "$FASTQ_FILE" ]; then
    echo "Error: FASTQ file not found at $FASTQ_FILE"
    exit 1
fi

# Create output directories if they don't exist
if [ ! -d "$FASTQC_DIR" ]; then
    echo "FastQC report directory doesn't exist. Creating it now..."
    mkdir -p "$FASTQC_DIR"
fi

if [ ! -d "$PLOTTING_DIR" ]; then
    echo "NanoPlot directory doesn't exist. Creating it now..."
    mkdir -p "$PLOTTING_DIR"
fi

# Run FastQC for quality control
echo "Running FastQC on $FASTQ_FILE..."

# I've specified the output directory for FastQC reports and added the `-o` option to direct the output there.
fastqc "$FASTQ_FILE" -o "$FASTQC_DIR"

# Check if FastQC ran successfully
if [ $? -eq 0 ]; then
    echo "FastQC completed successfully. Reports saved to $FASTQC_DIR"
else
    echo "Error: FastQC failed!"
    exit 1
fi

# Run NanoPlot for visualizing read lengths
echo "Running NanoPlot on $FASTQ_FILE..."

# I've added the `--outdir` flag to specify the output directory for the plot.
NanoPlot --fastq "$FASTQ_FILE" --outdir "$PLOTTING_DIR" --title "Read Length Distribution" --dpi 300

# Check if NanoPlot ran successfully
if [ $? -eq 0 ]; then
    echo "NanoPlot completed successfully. Plot saved to $PLOTTING_DIR"
else
    echo "Error: NanoPlot failed!"
    exit 1
fi

echo "Quality control analysis complete!"
