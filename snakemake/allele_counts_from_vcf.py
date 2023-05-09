
import numpy as np
import pandas as pd
import allel

### create a dictionary where 1 dataframe containing all the replicates for each optima is saved 

allele_counts = pd.DataFrame()
indiv_number = pd.DataFrame()
#  take all the vcf files in the folder 
for optima in ['optima0', 'optima1', 'optima2']: #os.listdir(dir_cluster):
    optima_subp = os.listdir(dir_cluster + optima + '/')
    subps = [i for i in optima_subp if '.vcf' in i ]
    
    for subp in subps:
        vcf = allel.read_vcf(dir_cluster + optima + '/' + subp)
        print(len(vcf['variants/POS']))
        genotypes = allel.GenotypeArray(vcf['calldata/GT'])
        ref_allele_counts = genotypes.count_alleles()[:, 0]
        alt_allele_counts = genotypes.count_alleles()[:, 1]

        subp0 = pd.DataFrame(data ={'ref_' + optima+subp[:-4]:ref_allele_counts, 'alt_'+ optima + subp[:-4]: alt_allele_counts}, index = vcf['variants/POS'], )
        number_of_ind = (subp0.iloc[0, 0] + subp0.iloc[0, 1]) / 2 
        print(number_of_ind)
        indiv_number[optima + subp[:-4]] = [number_of_ind]
        allele_counts = pd.concat([allele_counts, subp0], axis=1)


### if the alt is nan is because it is 0, and if the ref is nan is because it is the max number (number of individuals)

for col in allele_counts.columns:
    if 'ref_' in col: 
        allele_counts.loc[:,col] = allele_counts.loc[:, col].fillna(indiv_number.loc[0, col[4:]])
    elif 'alt_' in col:
        allele_counts.loc[:,col] = allele_counts.loc[:, col].fillna(0)


allele_counts = allele_counts.reset_index().rename(columns={'index': 'pos'})


allele_counts.to_csv('results_slim/pi0.001beta5/allele_counts.csv')

indiv_number.to_csv('results_slim/pi0.001beta5/pool_sizes.csv', index=None)