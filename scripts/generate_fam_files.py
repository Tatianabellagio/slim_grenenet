import numpy as np
import pandas as pd
import os

fam_file_input = snakemake.input['og_famfile'] 
ecotype_counts_file = snakemake.input['og_famfile'] 
fam_file_ouput = snakemake.output['fam_file_ouput'] 
optima = snakemake.params['optima']

optima_full = 'optima' + optima

og_fam = pd.read_csv(fam_file_input, sep = ' ', header=None)
## eliminate the 'fakephenotpye'
og_fam = og_fam.drop(5,axis=1)

ecotypes = pd.read_csv(ecotype_counts_file).drop('Unnamed: 0',axis=1)
ecotypes = ecotypes.set_index('ecotype')
## keep only the ones from this optima 
ecotypes = ecotypes[[i for i in ecotypes.columns if optima_full in i]]
## imgonna check the variance among subpop to see if i should run separate or one gwa
ecotypes = ecotypes.fillna(0)
## calculate the freq of extyopes in each subp 
ecotype_freq = pd.DataFrame()
for i in ecotypes.columns:
    ecotype_freq[i] = ecotypes[i] / ecotypes[i].sum()

ecotype_freq_mean = ecotype_freq.groupby(ecotypes.columns.str.split('/').str[4], axis=1).mean()
ecotype_freq_mean = ecotype_freq_mean.round(6)
## elimiante other since it is not on hte og vcf
ecotype_freq_mean = ecotype_freq_mean[ecotype_freq_mean.index != 'other']
ecotype_freq_mean.index = ecotype_freq_mean.index.astype(int)

final_fam = pd.concat([og_fam.set_index(0), ecotype_freq_mean],axis=1).reset_index()

final_fam.to_csv(fam_file_ouput, header=None, index=None, sep = ' ')