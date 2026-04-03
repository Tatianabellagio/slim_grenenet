# slim_grenenet

Snakemake pipeline for simulating and analyzing the [GrENE-Net](http://grene-net.org/) evolve-and-resequence experiment in *Arabidopsis thaliana* using [SLiM](https://messerlab.org/slim/). The pipeline models 231 ecotypes evolving across 32 environments over 6 generations, with 12 replicates per environment, tracking allele frequency dynamics under 30 different genetic architectures (varying polygenicity, effect sizes, and initial allele frequencies). Results are used to benchmark methods for detecting causal loci underlying local adaptation.

## Publication

Bellagio T et al. (2025). *Simulating the GrENE-Net experiment.* bioRxiv. [doi: 10.1101/2025.06.13.659553](https://www.biorxiv.org/content/10.1101/2025.06.13.659553v1)

---

## Repository structure

```
Snakefile              # Main pipeline definition
config.yaml            # Pipeline parameters (environments, replicates, architectures)
envs/                  # Conda environment specifications
scripts/               # Analysis and helper scripts
results/               # Pipeline outputs (generated)
log/                   # Snakemake logs
profiles/              # Cluster execution profiles (SLURM, local)
extra.slim             # Additional SLiM model components
```

---

## Usage

1. Clone the repository

2. Install Snakemake (recommended via Mamba):
   ```bash
   mamba create -n snakemake -c bioconda -c conda-forge snakemake
   conda activate snakemake
   ```

3. Edit `config.yaml` to set parameters for your run

4. Run the pipeline:
   ```bash
   snakemake --use-conda --cores <N>
   ```
   On a cluster (SLURM), use the provided profile:
   ```bash
   snakemake --use-conda --profile profiles/slurm
   ```

---

## Tools & dependencies

- [SLiM](https://messerlab.org/slim/) — forward-time population genetics simulator
- [Snakemake](https://snakemake.readthedocs.io/) — workflow management
- Python 3, R
