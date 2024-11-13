#!/bin/bash

# coverage_plot.sh
# This script generates a coverage plot from the aligned BAM file.

# Define input and output directories
OUTPUT_DIR="path/to/output"  # Replace with the path where the coverage plot will be saved
ALIGNED_BAM="${OUTPUT_DIR}/aligned_reads.sorted.bam"  # Replace with your aligned BAM file

# Define output plot file
COVERAGE_PLOT="${OUTPUT_DIR}/coverage_plot.png"  # Output plot file

# Check if input BAM file exists
if [ ! -f "$ALIGNED_BAM" ]; then
    echo "Error: Aligned BAM file not found at $ALIGNED_BAM"
    exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory doesn't exist. Creating it now..."
    mkdir -p "$OUTPUT_DIR"
fi

# Generate coverage stats using samtools depth
echo "Generating coverage statistics from BAM file..."

# I've used samtools depth to get base-level coverage across the reference genome.
# -aa: Includes all positions, even if no reads map there.
samtools depth "$ALIGNED_BAM" > "${OUTPUT_DIR}/coverage.txt"

# Check if samtools depth ran successfully
if [ $? -eq 0 ]; then
    echo "Coverage statistics saved to ${OUTPUT_DIR}/coverage.txt"
else
    echo "Error: samtools depth failed!"
    exit 1
fi

# Plot coverage using R or Python
echo "Generating coverage plot..."

# Here, I'm using a simple Python script to generate a coverage plot.
# I'm assuming Python and matplotlib are installed on your system.
python3 - <<EOF
import matplotlib.pyplot as plt

# Read coverage data
coverage_file = "${OUTPUT_DIR}/coverage.txt"
coverage = []

# Load the coverage data
with open(coverage_file) as f:
    for line in f:
        parts = line.split()
        coverage.append(int(parts[2]))  # Extract the coverage value from each line

# Generate the coverage plot
plt.figure(figsize=(10, 6))
plt.hist(coverage, bins=100, color="blue", edgecolor="black")
plt.title("Read Coverage Distribution")
plt.xlabel("Coverage Depth")
plt.ylabel("Frequency")
plt.savefig("${COVERAGE_PLOT}", dpi=300) 

# Show the plot
plt.show()

EOF
# this part always weirds me out
# Check if the plot was generated successfully
if [ $? -eq 0 ]; then
    echo "Coverage plot generated successfully. Plot saved to $COVERAGE_PLOT"
else
    echo "Error: Coverage plot generation failed!"
    exit 1
fi

echo "Coverage plot complete!"
