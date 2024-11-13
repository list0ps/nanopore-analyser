#!/usr/bin/env python3

# variant_viz.py
# This Python script visualizes the variants from a VCF file using matplotlib.
# The script generates a bar plot showing the distribution of variants across chromosomes.

import matplotlib.pyplot as plt
import vcf  # PyVCF library to parse VCF files

# Define input and output paths
VCF_FILE = "path/to/output/variant_calls.vcf"  # Replace with your VCF file path
OUTPUT_PLOT = "path/to/output/variant_distribution.png"  # Output plot file

# Check if the VCF file exists
try:
    vcf_reader = vcf.Reader(open(VCF_FILE, 'r'))
except FileNotFoundError:
    print(f"Error: VCF file not found at {VCF_FILE}")
    exit(1)

# Collect the number of variants per chromosome
chrom_variant_count = {}

# Iterate over the VCF records and count variants by chromosome
for record in vcf_reader:
    chrom = record.CHROM  # Get chromosome name (e.g., 'chr1', 'chr2', etc.)
    chrom_variant_count[chrom] = chrom_variant_count.get(chrom, 0) + 1

# Sort chromosomes by name (you may want to sort by length or other criteria)
sorted_chromosomes = sorted(chrom_variant_count.items())

# Prepare data for plotting
chromosomes, variant_counts = zip(*sorted_chromosomes)

# Plot the variant distribution across chromosomes
plt.figure(figsize=(12, 6))
plt.bar(chromosomes, variant_counts, color='skyblue', edgecolor='black')
plt.title("Variant Distribution Across Chromosomes")
plt.xlabel("Chromosome")
plt.ylabel("Number of Variants")
plt.xticks(rotation=90)  # Rotate chromosome labels for better readability

# Save the plot as a PNG file
plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300)

# Show the plot
plt.show()

# Check if the plot was saved successfully
print(f"Variant distribution plot saved to {OUTPUT_PLOT}")
