
import numpy as np
import pandas as pd
import os
import yaml
import re

pi =  str(snakemake.params['pi'])
print(pi)
beta = str(snakemake.params['beta'])
print(beta)
allele_freq = str(snakemake.params['allele_freq'])
print(beta)

outputs_list = snakemake.output
og_allele_freq_file = snakemake.input['og_allele_freq']
og_allele_freq = pd.read_csv(og_allele_freq_file).drop('Unnamed: 0',axis=1)
arq = 'arq_' + allele_freq + '_' + pi + '_' + beta
print('rerun')

def find_results_files(arq):
    # Find files named 'slim_output.txt' in subfolders for a particualr arq 
    file_paths = []
    for root, dirs, files in os.walk('results/'):
        for file in files:
            if 'slim_output.txt' in file:
                if arq in root:
                    file_path = os.path.join(root, file)
                    file_paths.append(file_path)
    return file_paths

## function to create allele_counts, and phenotypes for each generation inside slim_ouput file 
def allele_counts_phenotypes_per_gen(file, path_subp, subp_name):
    with open(file, "r") as input_file:
        content = input_file.read()
    generations = re.split(r"(generation\d+)", content)[1:]
    for i in range(0, len(generations), 2):
        generation_number = generations[i]
        generation_content = generations[i + 1].strip()
        output_file_path = f"{path_subp}{subp_name}_{generation_number}.txt"

        with open(output_file_path, "w") as output_file:
            output_file.write(generation_content)

    for i in range(2,7):
        file_path = f"{path_subp}{subp_name}_generation{i}.txt"
        with open(file_path, "r") as input_file:
            content = input_file.read()
        phenotypes = content.split('allele_counts')[0].split('phenotypes\n')[1]
        allele_counts = content.split('allele_counts')[1]
        with open(f"{path_subp}{subp_name}_generation{i}_phenotypes.txt", "w") as output_file:
            output_file.write(phenotypes)
        with open(f"{path_subp}{subp_name}_generation{i}_allelecounts.txt", "w") as output_file:
            output_file.write(allele_counts)

def unique_ecotypes(file, path_subp, subp_name):
    ecot_counts = pd.DataFrame()
    for i in range(2,7):
        if os.path.getsize(f"{path_subp}{subp_name}_generation{i}_phenotypes.txt") <=1:
            with open(path_subp + subp_name + "_ecot_counts.csv", "w") as f:
                f.write('empty')
        else:
            phenotype = pd.read_csv(f"{path_subp}{subp_name}_generation{i}_phenotypes.txt",sep = ',',header=None).T

            ## extract pop size 
            #with open(f"{path_subp}{subp_name}_generation{i}_pop_size.txt", "w") as output_file:
            #    output_file.write(str(pop_size))

            unique_pheno = phenotype[0].value_counts().reset_index()
            unique_pheno.columns = ['pheno', 'unique_eco_gen' + str(i)]
            if ecot_counts.empty:
                ecot_counts = unique_pheno
            else:
                ecot_counts = ecot_counts.merge(unique_pheno, on= 'pheno', how='outer')
            ecot_counts.to_csv(path_subp + subp_name + '_ecot_counts.csv')

def find_allelecounts_files(arq):
    ## collect all the files that are slim outputs 
    ## make sure the size is bigger than 580, because if not it might be corrupted 
    # Find files named 'slim_results' in subfolders
    file_paths = []
    for root, dirs, files in os.walk('results/'):
        for file in files:
            if '_allelecounts.txt' in file:
                if arq in root:
                    file_path = os.path.join(root, file)
                    file_paths.append(file_path)
    return file_paths

## create a dictionary where the architectures are the keys and the paths including
## and all the vcf files inside those architectures are the data

def create_dict_generation(file_paths_arq):
    files_per_gen = {}
    generations = ['generation'+ str(i) for i in range(2,7)]
    for gen in generations:
        vcfs_path = [file for file in file_paths_arq if gen in file]
        files_per_gen[gen] = vcfs_path
    return files_per_gen

def calc_allele_c_pop_size (files_per_arq_per_gen, og_allele_freq):
    ## main fucntion tgat calculates the alle count and pop sizes for each arq 
    pop_sizes = pd.DataFrame()
    allele_counts = pd.DataFrame({'pos': []})

    for path in files_per_arq_per_gen:
        start_index = path.find('optima')
        tag = path[start_index:].replace('_allelecounts.txt', '').replace('/', '_')

        #load data
        data = pd.read_csv(path,header=None,low_memory=False)


        ## extract allele counts
        allele_counts1 = data[data[0] != 'pop_size'].copy()    
        allele_counts1 = allele_counts1.rename(columns={0:'pos',1:tag})
        ##### correct for slim position error #######################
        allele_counts1['pos'] = allele_counts1['pos'].astype(int) + 1
        ##########################################################
        if not allele_counts1.empty:
            allele_counts = allele_counts.merge(allele_counts1, how='outer', on='pos')

        ## extract the pop sizes
        pop_size = data[data[0] == 'pop_size'][[1]].copy()
        pop_size.columns = [tag]
        pop_sizes = pd.concat([pop_sizes, pop_size],axis=1)

    ## filna with 0 since it means no copies of alt allele in that pos
    allele_counts = allele_counts.fillna(0)
    
    ### calcualte allele freq 
    allele_freq_evol = allele_counts.drop('pos',axis=1).div(pop_sizes.iloc[0] * 2, axis='columns')    
    
    ## drop columns where all individuals died so that the lineal features deostn get consufed with individuals dying at the extremes
    allele_freq_evol = allele_freq_evol.dropna(how = 'all', axis=1)
    
    allele_freq_evol['pos'] = allele_counts['pos'].astype(int)
    merged = pd.merge(allele_freq_evol, og_allele_freq, on = 'pos', how = 'outer')
    merged = merged.fillna(0)
    og_freq = merged['a_freq']
    merged = merged.set_index('pos')
    og_allelefreq = merged['a_freq']
    merged = merged.drop('a_freq',axis=1)
    delta_p = merged.sub(og_allelefreq, axis=0)
    norm_coef = og_allelefreq*(1-og_allelefreq)
    delta_p_norm = delta_p.div(norm_coef,axis=0)
    
    return allele_freq_evol,delta_p_norm, pop_sizes

## get all teh result files
file_paths = find_results_files(arq)

for file in file_paths:
    ## get the path of the subp 
    path_subp = file.split('subp')[0]
    ## get the name of the subp 
    subp_name = file.split('/')[-1].split('_')[0]
    ## create txt of phenotypes and txtof allele counts for each generation 
    allele_counts_phenotypes_per_gen(file, path_subp, subp_name)
    ## create the file of unique ecotpyes per generation with the phenotype file 
    unique_ecotypes(file, path_subp, subp_name)

## get all teh allele count files
allelecounts_files = find_allelecounts_files(arq)

## separate them by generation 
arq_gen_dict = create_dict_generation(allelecounts_files)

## create the csv for allele counts and pop sizes 
for i in arq_gen_dict.keys():
    outpus_this_gen = [j for j in outputs_list if i in j]
    allele_counts,delta_p_norm, pop_sizes = calc_allele_c_pop_size(arq_gen_dict[i], og_allele_freq)
    #output
    allele_counts.to_csv(outpus_this_gen[0], index=False)
    delta_p_norm.to_csv(outpus_this_gen[1], index=False)
    pop_sizes.to_csv(outpus_this_gen[2] , index=False)