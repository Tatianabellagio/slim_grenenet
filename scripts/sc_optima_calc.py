import numpy as np
import pandas as pd
import allel
import random


#### params 
# get the option of pi beta and alle freq chosen 
chr_number = snakemake.params['chr_number'] 
n_ecotypes = snakemake.params['n_ecotypes'] 
vcf_file = snakemake.input['base_vcf'] 

pi_option =  snakemake.params['pi']
beta = snakemake.params['beta']
print(beta)
alelle_freq_option = snakemake.params['allele_freq']

#get the actual values
n_optima = str(snakemake.params['n_optima']) 

pi = float(snakemake.params[pi_option])

allele_freq = snakemake.params[alelle_freq_option]
lower_bound = float(allele_freq[0])
upper_bound = float(allele_freq[1])

bed_sc = snakemake.output["bed_sc"]
optima_slim = snakemake.output["slim"]
######

## get vcf file
vcf_og = allel.read_vcf(vcf_file)

def scale(beta_values):
    min_value = np.min(beta_values)
    max_value = np.max(beta_values)
    beta_values = ((beta_values - min_value) / (max_value - min_value)) * (5 - 1) + 1
    return beta_values

def getting_selection_coef(vcf_og, pi, beta, lower_bound, upper_bound):
    ### calculation of the number of contributing positions and effect sizes of each ###################
    pos = vcf_og['variants/POS']
    ## get all the pos from the vcf file
    n_pos = len(pos)
    geno_og = vcf_og["calldata/GT"]
    ## freq of alternative alleles
    alt_al_count = geno_og.sum(axis=2).sum(axis=1)
    n_pos_contr = round(pi*n_pos)
    
    alelle_dist = pd.DataFrame({'alt_al_count':alt_al_count, 'pos':pos})

    alelle_dist['alt_al_freq'] = alelle_dist['alt_al_count'] / (n_ecotypes*2)

    if beta == 'betaprop': 
        alelle_dist = alelle_dist.copy()
        beta_values = alelle_dist.loc[:,('alt_al_freq')] * 5
        beta_values = scale(beta_values)
        alelle_dist.loc[:,'beta'] = beta_values

    elif beta == 'betanonprop':
        alelle_dist = alelle_dist.copy()
        beta_values = 5 - (alelle_dist.loc[:,('alt_al_freq')] * 5 )
        beta_values = scale(beta_values)
        alelle_dist.loc[:,'beta'] = beta_values

    alelle_dist = alelle_dist[(alelle_dist['alt_al_freq'] > lower_bound) & (alelle_dist['alt_al_freq'] <= upper_bound)]

    pos_contr = alelle_dist['pos'].sample(n_pos_contr)

    alelle_dist = alelle_dist[alelle_dist['pos'].isin(pos_contr)]

    effect_sizes = [np.random.normal(0, beta, 1) for beta in alelle_dist['beta']]
    effect_sizes = effect_sizes.copy()

    alelle_dist.loc[:,'effect_sizes'] = np.ravel(effect_sizes).copy()
    alelle_dist.loc[:,'effect_sizes'] = alelle_dist.loc[:,'effect_sizes'].round(4)
    effect_sizes = np.round(alelle_dist['effect_sizes'], 4)

    pos = pd.Series(pos)
    pos.name = 'pos'
    bed_file = pd.merge(pos, alelle_dist[['pos', 'effect_sizes']], on='pos', how='outer').fillna(0)
    bed_file['positions_from'] = bed_file['pos'] - 1
    bed_file['chr'] = chr_number
    bed_file = bed_file[['chr', 'positions_from', 'pos', 'effect_sizes']]

    return bed_file

def calc_optima(bed_file, vcf_og):
    geno_og = vcf_og["calldata/GT"]
    alt_al_per_pos = geno_og.sum(axis=2) #.sum(axis=1)
    phenotypes = []

    for i in range(alt_al_per_pos.shape[1]):
        gen_effectsize = np.multiply(alt_al_per_pos[:, i] , np.array(bed_file['effect_sizes']))
        phenotypes.append(gen_effectsize.sum())

    max_pheno = max(phenotypes)
    min_pheno = min(phenotypes)

    length = max_pheno - min_pheno
    step = length/(int(n_optima) - 1)
    optima = [min_pheno + i * step for i in range(0, int(n_optima))]
    return optima

bed_file = getting_selection_coef(vcf_og, pi, beta, lower_bound, upper_bound)
bed_file.to_csv(bed_sc, sep='\t', header=False, index=False) # selection_coef_chr5.bed

pheno_optima = calc_optima(bed_file, vcf_og)
pd.Series(pheno_optima).to_csv(optima_slim, header=None, sep='\n', index=None)
