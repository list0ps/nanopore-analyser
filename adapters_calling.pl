#!/usr/bin/env perl

# combines variant calling and adapter trimming (i hope it does at least)
# This Perl script combines adapter trimming and variant calling into one.
# First, it trims adapters from the raw reads using Porechop, then it calls variants using Medaka.

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);

# Define directories and files
my $RAW_DATA_DIR = "path/to/raw/data";   # Replace with the path to your raw FASTQ files
my $OUTPUT_DIR   = "path/to/output";     # Replace with the path where output files will be saved

my $RAW_FASTQ    = "$RAW_DATA_DIR/raw_reads.fastq";  # Input raw FASTQ file
my $FILTERED_FASTQ = "$OUTPUT_DIR/filtered_reads.fastq";  # Output filtered FASTQ file

my $REFERENCE_GENOME = "path/to/reference/genome.fasta";  # Path to the reference genome
my $ALIGNED_BAM      = "$OUTPUT_DIR/aligned_reads.sorted.bam";  # Output BAM file from alignment
my $VCF_OUTPUT       = "$OUTPUT_DIR/variant_calls.vcf";  # Output VCF file for variant calls

# Create output directory if it doesn't exist
make_path($OUTPUT_DIR) unless -d $OUTPUT_DIR;

# Step 1: Trim adapters using Porechop
print "Step 1: Trimming adapters using Porechop...\n";

# Check if raw FASTQ file exists
unless (-e $RAW_FASTQ) {
    die "Error: Raw FASTQ file $RAW_FASTQ not found.\n";
}

# Run Porechop to trim adapters
my $porechop_cmd = "porechop -i $RAW_FASTQ -o $FILTERED_FASTQ";
print "Running: $porechop_cmd\n";
my $porechop_output = `$porechop_cmd`;

# Check if Porechop ran successfully
if ($? == 0) {
    print "Adapter trimming completed successfully. Filtered reads saved to $FILTERED_FASTQ\n";
} else {
    die "Error: Porechop failed!\n";
}

# Step 2: Align reads using minimap2
print "Step 2: Aligning reads using minimap2...\n";

# Check if filtered FASTQ file exists
unless (-e $FILTERED_FASTQ) {
    die "Error: Filtered FASTQ file $FILTERED_FASTQ not found.\n";
}

# Run minimap2 for read alignment
my $minimap2_cmd = "minimap2 -ax map-ont $REFERENCE_GENOME $FILTERED_FASTQ | samtools view -Sb - > $ALIGNED_BAM";
print "Running: $minimap2_cmd\n";
my $minimap2_output = `$minimap2_cmd`;

# Check if minimap2 ran successfully
if ($? == 0) {
    print "Alignment completed successfully. Output BAM saved to $ALIGNED_BAM\n";
} else {
    die "Error: minimap2 failed!\n";
}

# Step 3: Sort and index the BAM file using samtools
print "Step 3: Sorting and indexing the BAM file...\n";

# Sort the BAM file
my $samtools_sort_cmd = "samtools sort $ALIGNED_BAM -o ${ALIGNED_BAM%.bam}.sorted.bam";
print "Running: $samtools_sort_cmd\n";
my $samtools_sort_output = `$samtools_sort_cmd`;

# Index the BAM file
my $samtools_index_cmd = "samtools index ${ALIGNED_BAM%.bam}.sorted.bam";
print "Running: $samtools_index_cmd\n";
my $samtools_index_output = `$samtools_index_cmd`;

# Check if sorting and indexing were successful
if ($? == 0) {
    print "BAM file sorted and indexed successfully.\n";
} else {
    die "Error: Sorting and indexing failed!\n";
}

# Step 4: Call variants using Medaka
print "Step 4: Calling variants using Medaka...\n";

# Check if the aligned BAM file exists
unless (-e $ALIGNED_BAM) {
    die "Error: Aligned BAM file $ALIGNED_BAM not found.\n";
}

# Run Medaka for variant calling
my $medaka_cmd = "medaka call -i ${ALIGNED_BAM%.bam}.sorted.bam -r $REFERENCE_GENOME -o $OUTPUT_DIR/medaka_output --model r941_min_high_g360";
print "Running: $medaka_cmd\n";
my $medaka_output = `$medaka_cmd`;

# Check if Medaka ran successfully
if ($? == 0) {
    print "Variant calling completed successfully. Output VCF saved to $VCF_OUTPUT\n";
} else {
    die "Error: Medaka failed!\n";
}
#hopefully this works.
print "Variant calling and adapter trimming complete.\n";
