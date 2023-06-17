import numpy as np
import pandas as pd
import os

pheno_og = snakemake.input['pheno_og'] 
fam_file_path = snakemake.input['og_famfile'] 

allele_freq = snakemake.params['allele_freq'] 
pi = snakemake.params['pi'] 
beta = snakemake.params['beta'] 

optima = snakemake.params['optima'] 
generation = snakemake.params['generation'] 

output_fam = snakemake.output[0]


fam_file = pd.read_csv(fam_file_path, sep = ' ',header=None)
pheno_og = pd.read_csv(pheno_og)
path = 'results/arq_' +  allele_freq + '_' + pi + '_' + beta + '/optima' + str(optima)

## get valid subpaths 
gen_name = 'generation' + str(generation) + '_phenotypes' 
subps = [i for i in os.listdir(path) if gen_name in i]


## calculate the mean freq for all populations 

all_subp = pd.DataFrame()
for i in subps:
    full_path = path + '/' + i
    supb_number = i.split('_')[0]
    ## check that the file is not empty 
    if os.path.getsize(full_path) > 1:
        
        eco_evolved = pd.read_csv(full_path,header=None).T.dropna()
        #print(eco_evolved)
        num_ind = len(eco_evolved)
        print(num_ind)
        eco_evolved = eco_evolved.value_counts() / num_ind
        eco_evolved.name = 'freq_' + supb_number
        eco_evolved = eco_evolved.reset_index()
    ## if it is empty create a dataframe where there is 0 freq of ecotypes
    else:
        data = {0: 0.0, 'freq_' + supb_number: 0.0}
        index = [0]  # Specify the desired index
        eco_evolved = pd.DataFrame(data, index=index)
    
    if all_subp.empty:
        all_subp =  eco_evolved
    else:
        all_subp = all_subp.merge(eco_evolved, on =0, how='outer')

all_subp = all_subp.set_index(0).mean(axis=1)
all_subp.name = 'freq'
all_subp = all_subp.reset_index()


eco_freq_og = pheno_og.merge(all_subp,how='left', left_on='0', right_on = 0 )

eco_freq_og['freq'] = eco_freq_og['freq'].fillna(0)

pheno = eco_freq_og['freq']

fam_file[5] = pheno

fam_file.to_csv(output_fam, sep = ' ',header=None, index= False)