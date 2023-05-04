#!/bin/bash
#SBATCH --partition=DPB
#SBATCH --job-name=haplotype_selection
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32g
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xwu@carnegiescience.edu
#SBATCH --error=/home/xwu/scratch/grenet/haplotype_selection/haplotype_selection_Gtest.err
# run the command

cd $SLURM_SUBMIT_DIR
module load SLiM/4.0.1

slim arabidopsis_evolve.slim