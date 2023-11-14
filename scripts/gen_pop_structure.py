import numpy as np
import pandas as pd

ecotype_counts = snakemake.input['ecotype_counts'] 
pc_founders = snakemake.input['pc_founders'] 

### ecotype freq normalized 
ecotype_counts = pd.read_csv(ecotype_counts)
## from ectoype counts to 
ecotype_counts = ecotype_counts.drop(columns = 'Unnamed: 0')
## calculate inital ecotype freq 
ecotypep0 = ecotype_counts['ecotype'].copy()
## all thee cotypes had the same initial freq 
initial_freq = 11/2541
ecotypep0 = pd.DataFrame({'ecotypep0': ecotypep0, 'freq': initial_freq})
## other has initial freq. 0 
ecotypep0.loc[ecotypep0['ecotypep0'] == 'other', 'freq'] = 0
ecotype_counts = ecotype_counts.set_index('ecotype')
total = ecotype_counts.sum(axis=0)
ecotypep0 = ecotypep0.set_index('ecotypep0')
ecotypep0['deno_norm'] = ecotypep0['freq'] * (1- ecotypep0['freq'])
ecotype_freq = ecotype_counts.div(total, axis=1)
ecotype_freq = ecotype_freq.fillna(0)
ecotype_freq = pd.concat([ecotypep0, ecotype_freq],axis=1)
ecotype_deltapn = pd.DataFrame()
## calculate delta p 
for col in ecotype_freq.columns:
    ecotype_deltapn[col] = (ecotype_freq[col] - ecotype_freq['freq']) / ecotype_freq['deno_norm']
ecotype_deltapn = ecotype_deltapn.drop(['freq', 'deno_norm'],axis=1)
ecotype_deltapn = ecotype_deltapn.drop('other')
ecotype_deltapn_t = ecotype_deltapn.T

### matrix multiplciation between the  ecotype freq and the pcs, constrcut population structure correction 
pc_founders = pd.read_csv(pc_founders, sep = ' ', header=None)

pc_founders = pc_founders[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]].rename(columns = {1: 'ecotype'})
pc_founders = pc_founders.set_index('ecotype')

pop_structure = np.dot(ecotype_deltapn_t, pc_founders)
pop_structure = pd.DataFrame(pop_structure)

pop_structure.columns = [f'PC{i}' for i in range(1, 11)]
pop_structure.to_csv(pop_structure_file)