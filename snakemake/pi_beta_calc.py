import numpy as np
import pandas as pd
import allel
import random
import os
import sys
import snakemake

#### params 

chr_number = snakemake.input['chr_number'] 
vcf_file = snakemake.input['base_vcf'] 
optimum_ecotypes = snakemake.input['optimum_ecotypes']


### pi is 0.0001 since length_chr5 = 1702174 * 0.0001 = 170
pi =  snakemake.params['pi'] 
beta = snakemake.params['beta'] 
bed_sc = snakemake.output["bed_sc"]
######

### calculation of the number of contributing positions and effect sizes of each ###################
file_p = allel.read_vcf(vcf_file, fields='variants/POS')
## get all the pos from the vcf file
file_p = allel.read_vcf(vcf_file, fields='variants/POS')
pos = file_p['variants/POS']
n_pos = len(pos)
n_pos_contr =round(pi*n_pos)
## randomly select the positions that are going to contribute 
pos_contr = random.sample(pos.tolist(), k=n_pos_contr)
## now im gonna sample the number of nps contributing from a normal distribution with mean 0 and sd defined by beta 
print('For this run the value of pi is ' + str(pi) + ' and the value of beta is ' + str(beta) + ' based on the number of positions being ' + str(n_pos) + ' the number of contributing SNPs is ' + str(n_pos_contr))
mu, sigma = 0,beta # mean and standard deviation
effect_sizes = np.random.normal(mu, sigma, n_pos_contr)
effect_sizes = np.round(effect_sizes, 4)
## now i have to create a bed file to add this infromation to the vcf file 
## create dataframe of the contributing pos 
## to caluclate the positions FROM in a bed file, the fist value is not included
positions_from = [i-1 for i in pos]
# Define the positions and values
chromosome = [chr_number] * len(pos)
# Create a DataFrame from the positions and values
all_pos = pd.DataFrame({'chromosome':chromosome,'positions_from': positions_from, 'positions_to': pos})
contrib_pos = pd.DataFrame({'positions_to': pos_contr, 'sel_coef':effect_sizes})
# merge and fillna with 0 since that is the selection coefficient for all the non contributing snps 
all_pos = all_pos.merge(contrib_pos, left_on= 'positions_to', right_on= 'positions_to' ,how='left').fillna(0)
# Write the DataFrame to a BED file

## output 
all_pos.to_csv(bed_sc, sep='\t', header=False, index=False) # selection_coef_chr5.bed
print('Finished creating bed file to annotate vcf with effect sizes at each contributing position')


############################### filling the optimum_ecotypes file with phenotpyes based on the snps contributing and their effect sizes previously calculated
print('Starting optimum phenotype calculation based on ecotypes from an env gradient and the effect sized and SNPs contributing')
samples = pd.read_csv(optimum_ecotypes,usecols = ['ecotypeid', 'bio1'])
selection_coef = pd.read_csv(bed_sc, sep = '\t',header = None)[3]
vcf = allel.read_vcf(vcf_file)
phenotypes = []
for ecotype in samples['ecotypeid']:
    vcf_ecotype = allel.read_vcf(vcf_file, samples=[str(ecotype)])
    ## for each position this is the number of alterantive variants 
    alt_alleles_per_pos = vcf_ecotype['calldata/GT'].sum(axis=2)
    gen_effectsize = np.multiply(alt_alleles_per_pos.flatten(), np.array(selection_coef))  ## select coef are actually effect sizes
    phenotypes.append(gen_effectsize.sum())
samples['phenotype'] = phenotypes
samples['phenotype'] = samples['phenotype'].round(4)

# write optimas and also create the folder where the vcf fiels will be dumped 
# Check if the folder already exists
if os.path.exists("vcf_slim"):
    # If it does, delete the folder and all its contents
    os.system("rm -r vcf_slim")
    os.makedirs("vcf_slim")
    os.system("rm -r optima_files")
    os.makedirs('optima_files')
else: 
    os.makedirs("vcf_slim")
    os.makedirs('optima_files')

samples.to_csv('optima_files/optimas.csv')
samples['phenotype'].to_csv('optima_files/optima_slim.txt', header=None, sep='\n', index=None)
print(f"Folders: optima_files with optimas.csv and optimas for slim created successfully.")
for i in range(0,34):
    path = os.path.join("vcf_slim", str(f'optima{i}'))
    os.makedirs(path)
print(f"Folders for saving vcf files generated by slim created successfully.")
### end of calculating phenotypes and creating folders of optimas and folders where future vcf file gen by slim will be saved 


