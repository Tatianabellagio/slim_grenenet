
import numpy as np
import pandas as pd
import os
import yaml

pi =  str(snakemake.params['pi'])
print(pi)
beta = str(snakemake.params['beta'])
print(beta)
pool_sizes_output = snakemake.output['pool_sizes']
print(pool_sizes_output)
allele_counts_output = snakemake.output['allele_counts']
print(allele_counts_output)

arq = f'arq_pi{pi}_beta{beta}'
print(arq)

def find_results_files():
    ## collect all the files that are slim outputs 
    ## make sure the size is bigger than 580, because if not it might be corrupted 
    min_file_size = 580
    # Find files named 'slim_results' in subfolders
    file_paths = []
    for root, dirs, files in os.walk('results/'):
        for file in files:
            if 'slim_output.txt' in file:
                file_path = os.path.join(root, file)
                if os.path.getsize(file_path) > min_file_size:
                    file_paths.append(file_path)
    return file_paths

## create a dictionary where the architectures are the keys and the paths including
## and all the vcf files inside those architectures are the data
def create_dict_arq(arqs, file_paths):
    files_per_arq = {}
    for arq in arqs:
        vcfs_path = [file for file in file_paths if arq in file]
        files_per_arq[arq] = vcfs_path
    return files_per_arq


def calc_allele_c_pop_size (files_per_arq, arq):
    ## main fucntion tgat calculates the alle count and pop sizes for each arq 
    pop_sizes = pd.DataFrame()
    allele_counts = pd.DataFrame({'pos': []})

    for path in files_per_arq[arq]:
        start_index = path.find('optima')
        tag = path[start_index:].replace('_slim_output.txt', '').replace('/', '_')

        #load data
        data = pd.read_csv(path, skiprows=list(range(0,15)),header=None,low_memory=False)


        ## extract allele counts
        allele_counts1 = data[data[0] != 'pop_size'].copy()    
        allele_counts1 = allele_counts1.rename(columns={0:'pos',1:tag})
        if not allele_counts1.empty:
            allele_counts = allele_counts.merge(allele_counts1, how='outer', on='pos')

        ## extract the pop sizes
        pop_size = data[data[0] == 'pop_size'][[1]].copy()
        pop_size.columns = [tag]
        pop_sizes = pd.concat([pop_sizes, pop_size],axis=1)

    ## filna with 0 since it means no copies of alt allele in that pos
    allele_counts = allele_counts.fillna(0)
    return pop_sizes, allele_counts

file_paths = find_results_files()
files_per_arq = create_dict_arq([arq], file_paths)
pop_sizes, allele_counts = calc_allele_c_pop_size(files_per_arq, arq)

#output
pop_sizes.to_csv(pool_sizes_output, index=False)
allele_counts.to_csv(allele_counts_output, index=False)