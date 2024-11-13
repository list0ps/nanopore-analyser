# ONT Data Processing and Analysis for Blowflies
**under construction**
This repository contains a set of scripts and tools for processing and analyzing Oxford Nanopore Technology (ONT) sequencing data from blowflies. The pipeline includes data preprocessing, quality control (QC), alignment to a reference genome, variant calling, and visualization of results. This project is designed to help researchers efficiently analyze and interpret ONT sequencing data from blowflies.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Data Preprocessing](#data-preprocessing)
- [Quality Control](#quality-control)
- [Alignment](#alignment)
- [Variant Calling](#variant-calling)
- [Visualization](#visualization)
- [Results](#results)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Dependencies
The analysis pipeline depends on several bioinformatics tools and libraries. To ensure compatibility, it is recommended to use a conda environment. The following tools are required:

- `minimap2` (for read alignment)
- `medaka` (for variant calling)
- `NanoPlot` (for quality control)
- `porechop` (for adapter trimming)
- `samtools` (for BAM file handling)

You can set up the environment using the provided `environment.yml` file:

```bash
conda env create -f environment.yml
conda activate blowfly-ont
