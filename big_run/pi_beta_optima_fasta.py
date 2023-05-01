import numpy as np
import pandas as pd
import allel
import random
import os
import sys

#### default params 
## params for annotating vcf file with sc
chr_number = 5
vcf_file = 'chr5_grenenet.vcf'

### pi is 0.0001 since length_chr5 = 1702174 * 0.0001 = 170
pi =  0.01
beta = 5

#params for calculating optima 
grenenet_1001g_ecotypes = 'ecotypes_grenenet_1001g.txt'
# cluster
# safedata/ath_evo/grenephase1/data/worldclim_ecotypesdata.csv
# local 
# worldclim_ecotypesdata.csv
path_worldclim_ecotypesdata = 'worldclim_ecotypesdata.csv'
bed_sc = 'selection_coef_chr5.bed'

## params for gen of fasta from vcf 
or_fasta_file = 'chr5.fasta'
slim_fasta_file = 'slim_' + or_fasta_file

########### params defined in command line ############################
## the 2 params that i am going to be iterating through are pi and beta 

# Check if a command-line argument was provided

if len(sys.argv) > 1:
    # If an argument was provided, use it to set the variable
    pi = sys.argv[1]

if len(sys.argv) > 2:
    # If an argument was provided, use it to set the variable
    beta = sys.argv[2]

print('For this run the value of pi is ' + str(pi) + ' and the value of beta is ' + str(beta))
########### params defined in command line ############################

file_p = allel.read_vcf(vcf_file, fields='variants/POS')
## get all the pos from the vcf file
file_p = allel.read_vcf(vcf_file, fields='variants/POS')
pos = file_p['variants/POS']
n_pos = len(pos)
n_pos_contr =round(pi*n_pos)
## randomly select the positions that are going to contribute 
pos_contr = random.sample(pos.tolist(), k=n_pos_contr)
## now im gonna sample the number of nps contributing from a normal distribution with mean 0 and sd defined by beta 
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
all_pos.to_csv(bed_sc, sep='\t', header=False, index=False)


ecotypes_grenenet_1001g = pd.read_csv(grenenet_1001g_ecotypes).columns.astype(int)
ecotypes_temp = pd.read_csv(path_worldclim_ecotypesdata, usecols=['ecotypeid', 'bio1'])
#BIO1 = Annual Mean Temperature
## filter by the ones coming from the 1001 genomes
ecotypes_temp = ecotypes_temp[ecotypes_temp['ecotypeid'].isin(ecotypes_grenenet_1001g)]
### get the min and max values just for completeness of the range and then sample the rest 

# Find the minimum and maximum elements
ecot_max = ecotypes_temp[ecotypes_temp['bio1'] == ecotypes_temp['bio1'].max()]
ecot_min = ecotypes_temp[ecotypes_temp['bio1'] == ecotypes_temp['bio1'].min()]
# Remove the minimum and maximum elements from the dataframe 
ecotypes_temp = ecotypes_temp[~ecotypes_temp['bio1'].isin([ecotypes_temp['bio1'].max(),ecotypes_temp['bio1'].min()])].copy()

## and then subsample the rest 
samples = ecotypes_temp.sample(32)
## concat all 
samples = pd.concat([samples,ecot_max,ecot_min])
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
samples = samples.reset_index(drop=True)



# write optimas and also create the folder where the vcf fiels will be dumped 


# Check if the folder already exists
if os.path.exists("vcf_slim"):
    # If it does, delete the folder and all its contents
    os.system("rm -r vcf_slim")
    os.makedirs("vcf_slim")
else: 
    os.makedirs("vcf_slim")


#os.makedirs("vcf_slim")

os.makedirs('optima_files')
samples.to_csv('optima_files/optimas.csv')
samples['phenotype'].to_csv('optima_files/optima_slim.txt', header=None, sep='\n', index=None)

for i in range(0,34):
    path = os.path.join("vcf_slim", str(f'optima{i}'))
    os.makedirs(path)
    print(f"Folder '{folder_name}' and subfolder 'optima{i}' created successfully.")



### generation of reference fasta from tair fasta with all chr

# Execute this block of code only if the file does not exist aka if we dont have the fasta file in slim format for the chr
if not os.path.exists(slim_fasta_file):
    ## save the first line (chr name) and all the other lines (seq) separatedly 
    with open(or_fasta_file, 'r') as file:
        chro = ""
        seq = ""
        for i, line in enumerate(file):
            if i == 0:
                chro += line.strip()
            elif i != 0:
                seq += line.strip()
    ## replace unknown variant with acgt so slim can read it 
    replacement_dict = {"M": "A", "R": "A", "W": "A", "S": "C", "Y": "C", "K": "G", "V": "A", "H": "A", "D": "A", "B": "C", "N": "A"}

    for key in replacement_dict:
        seq = seq.replace(key, replacement_dict[key])

    with open(slim_fasta_file, "w") as f:
        f.write(f'{chro}\n')
        f.write(seq + '\n')
        # Write the string to the file