import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

allele_freq = snakemake.input['allele_freq'] 
env_var_file = snakemake.input['env_var_lmm'] 
output_env = snakemake.output['env_var_lfmm'] 
ouput_n_comp = snakemake.output['num_components'] 

allele_freq = pd.read_csv(allele_freq)#.drop('Unnamed: 0',axis=1)
allele_freq = allele_freq.fillna(0)
allele_freq = allele_freq.drop('chrom_pos',axis=1)
n_components = len(allele_freq.columns)

scaler = StandardScaler()
scaler.fit(allele_freq)

scaled = scaler.fit_transform(allele_freq)

scaled_afn = pd.DataFrame(scaled, columns=allele_freq.columns)
# Perform PCA
pca = PCA(n_components=n_components)  # Set the desired number of components
pca.fit(scaled_afn)
# Get the transformed data (projected onto the principal components)
transformed_data = pca.transform(scaled_afn)
explain_var_ratio = pca.explained_variance_ratio_
cumulative_variance = np.cumsum(explain_var_ratio)
# Set the desired cumulative explained variance threshold
threshold = 0.96

# Find the number of components that exceed the threshold
num_components = np.sum(cumulative_variance <= threshold) + 1

with open(ouput_n_comp, 'w') as file:
    file.write(str(num_components))

env = pd.read_csv(env_var_file, index_col = [0])['env']
env.to_csv(output_env,header=None, sep = ' ', index=False)