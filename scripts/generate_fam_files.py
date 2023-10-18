import numpy as np
import pandas as pd
import os

fam_file_input = snakemake.input['og_famfile'] 
ecotype_counts_file = snakemake.input['og_famfile'] 
fam_file_ouput = snakemake.output['fam_file_ouput'] 


og_fam = pd.read_csv(fam_file_input, sep = ' ', header=None)

## eliminate the 'fakephenotpye'
og_fam = og_fam.drop(5,axis=1)

ecotypes = pd.read_csv(ecotype_counts_file).drop('Unnamed: 0',axis=1)

ecotypes_counts = ecotypes.set_index('ecotype').sum(axis=1)

ecotypes_counts = ecotypes_counts[ecotypes_counts.index != 'other']

ecotypes_counts.index = ecotypes_counts.index.astype(int)

final_fam = pd.concat([og_fam.set_index(0), ecotypes_counts],axis=1).reset_index()

final_fam.to_csv(fam_file_ouput, header=None, index=None, sep = ' ')