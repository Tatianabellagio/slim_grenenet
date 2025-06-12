import numpy as np
import pandas as pd

ecotype_counts = snakemake.input['ecotype_counts'] 
pc_founders = snakemake.input['pc_founders'] 
pop_structure_file = snakemake.output['pop_structure'] 
mapenv = snakemake.params['mapenv'] 
env_file = snakemake.output['env_variable'] 


### ecotype freq normalized 
ecotype_counts = pd.read_csv(ecotype_counts)
## from ectoype counts to 
ecotype_counts = ecotype_counts.drop(columns = 'Unnamed: 0')
## only create pop strcuture for the population where nto all the indivuals died 
pop_alive = ecotype_counts.sum()[ecotype_counts.sum() != 0 ].index
## filter by the populations that are alive 
ecotype_counts = ecotype_counts[pop_alive]
ecotype_counts = ecotype_counts.set_index('ecotype')
## calculate ecotype freq
total = ecotype_counts.sum(axis=0)
ecotype_freq = ecotype_counts.div(total, axis=1)
ecotype_freq = ecotype_freq.fillna(0)

ecotype_freq = ecotype_freq.drop('other')


## create env variable  
env = ecotype_freq.columns.str.split('_subp').str[0].str.replace('optima', '').astype(int)
mapenv = eval(mapenv)
sites = env.map(mapenv)
repl = ecotype_freq.columns.str.split('_subp').str[1].astype(int)
sites_env = pd.DataFrame({'sites': sites, 'repl': repl, 'env': env})

# save sites ande nv varaibel 
sites_env.to_csv(env_file)


## not transfort it to multiply 
ecotype_freq = ecotype_freq.T

### matrix multiplciation between the  ecotype freq and the pcs, constrcut population structure correction 
pc_founders = pd.read_csv(pc_founders, sep = ' ', header=None)

pc_founders = pc_founders.drop(0,axis=1).rename(columns = {1: 'ecotype'})
pc_founders = pc_founders.set_index('ecotype')

pop_structure = np.dot(ecotype_freq, pc_founders)
pop_structure = pd.DataFrame(pop_structure)

pop_structure.columns = [f'PC{i}' for i in range(1, 21)]
pop_structure.to_csv(pop_structure_file)

