import allel
import numpy as np
import pandas as pd

vcf_file_path = snakemake.input["vcf_file"]
bed_file_path = snakemake.input["bed_sc"]
output_pheno_og_path = snakemake.output["output_pheno"]

def phenos(bed_file, vcf_og):
    geno_og = vcf_og["calldata/GT"]
    alt_al_per_pos = geno_og.sum(axis=2) #.sum(axis=1)
    phenotypes = []

    for i in range(alt_al_per_pos.shape[1]):
        gen_effectsize = np.multiply(alt_al_per_pos[:, i] , np.array(bed_file[3]))
        phenotypes.append(gen_effectsize.sum())

    return phenotypes

vcf_og = allel.read_vcf(vcf_file_path)
bed_file = pd.read_csv(bed_file_path, sep = '\t', header=None,usecols=[2,3])
phenotypes_og = phenos(bed_file, vcf_og)
phenotypes_og = pd.DataFrame(phenotypes_og)

# write csv fro p=original phenotypes
phenotypes_og.to_csv(output_pheno_og_path, index=False)