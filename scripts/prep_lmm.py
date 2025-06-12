import numpy as np
import pandas as pd

ecotype_counts = snakemake.input['ecotype_counts'] 
allele_freq_file = snakemake.input['allele_freq'] 
allele_freq_founder = snakemake.input['allele_freq_founder_offset'] 
pc_founders = snakemake.input['pc_founders'] 
mapenv = snakemake.params['mapenv'] 

env_file = snakemake.output['env_variable'] 
p_norm_file = snakemake.output['allele_freq_norm'] 
pop_structure_file = snakemake.output['pop_structure'] 

## import allele freq and founders allele freq 
allele_freq = pd.read_csv(allele_freq_file).drop(columns = 'Unnamed: 0')
print(allele_freq.head())
allele_freq_founder = pd.read_csv(allele_freq_founder)
## convert allele allele_freq_founder chrom pos to int so it can be merged with allele_freq
allele_freq_founder['chrom_pos'] = allele_freq_founder['chrom_pos'].astype(int)
## eliminate duplciates basically positions where the allele freq is the same in all 
#allele_freq = allele_freq.round(10)
#allele_freq = allele_freq.set_index('chrom_pos').drop_duplicates()
## convert back chrom pos to a column
#allele_freq = allele_freq.reset_index()
## all the values where there is no allele freq is because thea llele freq is basically 0 
allele_freq = allele_freq.fillna(0)
## first calculate delta p norm 
## imoport allele freq of the founder and normalize 
allele_freq = pd.merge(allele_freq,allele_freq_founder, on ='chrom_pos')
allele_freq = allele_freq.set_index('chrom_pos')
print(allele_freq.head())
p_norm = pd.DataFrame()
for col in allele_freq.columns:
    p_norm[col] = (allele_freq[col] - allele_freq['allele_freq_founder']) / allele_freq['deno_norm']

p_norm = p_norm.drop(['allele_freq_founder','deno_norm'],axis=1)

#p_norm = p_norm.round(10)
## eliminate rows with all the same valuesbecause the linear model will just nor tun 
p_norm = p_norm[p_norm.round(5).std(axis=1) > 0]
print(p_norm.head())


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

pc_founders = pc_founders[[1,2,3,4]].rename(columns = {1: 'ecotype'})
pc_founders = pc_founders.set_index('ecotype')

pop_structure = np.dot(ecotype_deltapn_t, pc_founders)
pop_structure = pd.DataFrame(pop_structure)
pop_structure.columns = ['PC1', 'PC2', 'PC3']


## create env variable  
env = p_norm.columns.str.split('_subp').str[0].str.replace('optima', '').astype(int)
mapenv = eval(mapenv)
sites = env.map(mapenv)
repl = p_norm.columns.str.split('_subp').str[1].astype(int)
sites_env = pd.DataFrame({'sites': sites, 'repl': repl, 'env': env})

# save sites ande nv varaibel 
sites_env.to_csv(env_file)
# save p norm 
p_norm.to_csv(p_norm_file)
## save pop structure
pop_structure.to_csv(pop_structure_file)