# slim_grenenet

Snakemake + SLiM pipeline for forward-in-time population genetic simulations studying how genetic architecture influences evolutionary rescue under rapid environmental change.

The central question: does polygenic or monogenic trait architecture better enable populations to adapt and avoid extinction when the environment shifts? The simulations are seeded with real genomic variation from 231 *Arabidopsis thaliana* ecotypes (GrENE-Net founder population, ~3.2M SNPs) to ensure a biologically realistic genomic structure, but the question being tested is about the general relationship between polygenicity, heritability, and population survival, not a direct replay of the GrENE-Net experiment.

The simulations are also fast, despite this realistic genomic background, thanks to tree-sequence structures: the real founder VCF is converted once into a tree sequence (via `tsinfer`) carrying only the adaptive mutations, so SLiM tracks just those loci forward in time rather than the full ~3.2M-SNP genome.

## Publication

Bellagio T & Exposito-Alonso M (2026). *Polygenic and monogenic adaptation drive evolutionary rescue at different magnitudes of environmental change.* Journal of Heredity, esag026. [doi: 10.1093/jhered/esag026](https://doi.org/10.1093/jhered/esag026)

---

## What the simulations test

A quantitative trait under stabilizing selection is simulated across a full factorial design:

| Parameter | Levels |
|-----------|--------|
| Polygenicity (number of loci) | 1, 2, 5, 10, 20, 50, 100, 500, 1000 |
| Heritability (h²) | 0.1, 0.3, 0.5, 0.7, 0.9 |
| Magnitude of environmental shift | 0–7 phenotypic standard deviations |
| Independent genetic architectures per polygenicity level | 30 |
| Replicates per combination | 5 |

Total: **47,250 simulations**, each tracking a population for 10 generations. At each generation, population size, genotypes, phenotypes, fitness values, and full tree sequences are recorded.

Key findings:
- Under small to moderate environmental shifts, high polygenicity increases evolutionary rescue probability.
- Under extreme shifts, polygenic traits lead predictably to extinction, while monogenic traits can occasionally produce a single winning adaptive genotype.

---

## Biological grounding

- Population size: 2,310 individuals (carrying capacity 900), based on GrENE-Net field observations
- 97% selfing rate (*A. thaliana*)
- Recombination rate ~4 cM/Mb
- Initial allele frequencies sampled from the *A. thaliana* site frequency spectrum
- Effect sizes inversely proportional to allele frequency (or drawn from N(0,2) — both give the same qualitative results)

---

## Repository structure

```
Snakefile              # Pipeline definition — runs all simulation combinations
config.yaml            # Parameters (polygenicity levels, heritability, shifts, replicates)
envs/                  # Conda environment specifications
scripts/               # Analysis scripts and helpers
extra.slim             # Additional SLiM model components
results/               # Pipeline outputs (generated)
log/                   # Snakemake logs
profiles/              # Cluster execution profiles (SLURM, local)
```

---

## Usage

1. Clone the repository

2. Install Snakemake (recommended via Mamba):
   ```bash
   mamba create -n snakemake -c bioconda -c conda-forge snakemake
   conda activate snakemake
   ```

3. Provide the starting VCF (GrENE-Net founder population) and edit `config.yaml`

4. Run:
   ```bash
   snakemake --use-conda --cores <N>
   # On a cluster:
   snakemake --use-conda --profile profiles/slurm
   ```

---

## Tools & dependencies

- [SLiM](https://messerlab.org/slim/) — forward-time population genetics simulator with tree-sequence recording
- [Snakemake](https://snakemake.readthedocs.io/) — workflow management
- [tskit](https://tskit.dev/) / [msprime](https://msprime.readthedocs.io/) — tree sequence analysis
- Python 3, R
