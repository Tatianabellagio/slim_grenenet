import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

## snakemake 

input_allele_freq = snakemake.input['allele_freq_table'] 
input_pop_size = snakemake.input['pool_sizes'] 
input_optima = snakemake.input['optimas'] 

output_geno = snakemake.output['gen'] 
output_env = snakemake.output['env'] 
ouput_n_comp = snakemake.output['num_components'] 

n_replicates = int(snakemake.params['n_replicates'])
## snakemake 

allele_freq = pd.read_csv(input_allele_freq)
allele_freq = allele_freq.drop('pos',axis=1)


####calculating mean ########## 
## calculating the mean across subp insde optima
allele_freq = allele_freq.groupby(allele_freq.columns.str.split('_').str[0], axis=1).mean()
##################

n_components = len(allele_freq.columns)

scaler = StandardScaler()
scaler.fit(allele_freq)
scaled = scaler.fit_transform(allele_freq)

scaled_allele_freq = pd.DataFrame(scaled, columns=allele_freq.columns)

# Perform PCA
pca = PCA(n_components=n_components)  # Set the desired number of components
pca.fit(scaled_allele_freq)

# Get the transformed data (projected onto the principal components)
transformed_data = pca.transform(scaled_allele_freq)

explain_var_ratio = pca.explained_variance_ratio_

cumulative_variance = np.cumsum(explain_var_ratio)

# Set the desired cumulative explained variance threshold
threshold = 0.96

# Find the number of components that exceed the threshold
num_components = np.sum(cumulative_variance <= threshold) + 1

with open(ouput_n_comp, 'w') as file:
    file.write(str(num_components))

allele_freq = allele_freq.astype(float)

allele_freq.to_csv(output_geno,header=None, index=False)

#env : i will use the optima

env = pd.read_csv(input_optima,header=None)

## filter the environments only do the ones that have subpopulations that made it 
#######
optima_in_allelec = [int(i.replace('optima', '')) for i in allele_freq.columns]
env = env[env.index.isin(optima_in_allelec)]
#######

##scale it
env[0] = (env[0] - np.mean(env[0])) / np.std(env[0])

## to account for replciates in each enviornment 
#env = [element for element in env[0] for _ in range(n_replicates)] ##comment since now is the mean

env[0] = env[0].round(4)
env.to_csv(output_env,header=None, sep = ' ', index=False)