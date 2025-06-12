# SLiM Simulation Pipeline for Population Genetics Analysis

This repository contains a Snakemake pipeline for running population genetics simulations using SLiM. The pipeline is designed to simulate evolutionary scenarios and analyze the resulting data.

## Repository Structure

```
.
├── Snakefile              # Main Snakemake workflow file
├── config.yaml           # Configuration parameters for simulations
├── scripts/              # Python and R scripts for simulation and analysis
├── treeseq/             # Scripts for VCF to tree sequence conversion
├── analysis/            # Analysis scripts and notebooks
└── results/             # Simulation results 
```

## Workflow Overview

1. **Initial Setup**
   - Convert VCF files to tree sequence format using scripts in the `treeseq/` folder
   - This tree sequence file serves as the starting point for simulations

2. **Simulation Pipeline**
   - The main workflow is defined in `Snakefile`
   - All simulation parameters are configured in `config.yaml`
   - The pipeline uses scripts from the `scripts/` folder to:
     - Build initial populations
     - Run SLiM simulations
     - Process tree sequences
     - Calculate allele frequencies
     - Perform statistical analyses

3. **Analysis**
   - Results are collected and processed in the `analysis/` folder
   - Jupyter notebooks are provided for reproducing the analyses presented in the paper
   - The `scraping_results` folder contains scripts to aggregate simulation results

## Configuration

All simulation parameters are defined in `config.yaml`, including:
- Heritability values
- Selection strengths
- Number of replicates
- Environmental optima
- Population parameters

## Usage

1. **Setup Environment**
   ```bash
   conda env create -f envs/base_env.yaml
   conda activate base_env
   ```

2. **Run Pipeline**
   ```bash
   snakemake --use-conda
   ```

3. **Analysis**
   - Navigate to the `analysis/` folder
   - Run the provided Jupyter notebooks to reproduce the paper's analyses

## Dependencies

- Python 3.x
- R
- SLiM
- Snakemake
- Additional dependencies listed in `envs/base_env.yaml`

## Paper Reference

[Link to the paper will be added here]

## Citation

If you use this pipeline in your research, please cite:
[Citation information will be added here] 