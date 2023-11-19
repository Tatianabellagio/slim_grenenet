import numpy as np
import pandas as pd
import allel
import os
from collections import defaultdict
import multiprocessing

## here all the vcfs from a particular combination of arq herit and selection and expanding all replicates and optimas
output_vcf_offset = snakemake.input['output_vcf_offset'] 
nonhet_pos = snakemake.input['nonhet_pos'] 
og_vcf_offset = snakemake.input['og_vcf_offset']
ecotypes_grenenet = snakemake.input['ecotypes_grenenet']
ecotype_counts = snakemake.output['ecotype_counts']

def filtering_pos (nonhet_pos, pos_new, geno_og, geno_new):
    pos_to_keep = np.intersect1d(nonhet_pos, pos_new)
    mask_pos_ogvcf = pd.Series(pos_og).isin(pos_to_keep)
    geno_og_rpos  = geno_og[mask_pos_ogvcf]
    mask_pos_newvcf = pd.Series(pos_new).isin(pos_to_keep)
    geno_new_rpos  = geno_new[mask_pos_newvcf]
    return geno_og_rpos, geno_new_rpos

def get_ecotype_geno_mapper(geno_og_rpos):
    geno_og_rpos = np.swapaxes(geno_og_rpos, 0, 1)
    ecotype_geno_mapper = {}
    for i,j in zip(geno_og_rpos, samples):
        geno = i.tobytes()
        ecotype_geno_mapper[geno] = j
    return ecotype_geno_mapper

def get_ecotype_counts(geno_new_rpos, pop_name, ecotype_geno_mapper):
    geno_new_rpos = np.swapaxes(geno_new_rpos, 0, 1)
    # Initialize a defaultdict to store the genotype counts
    ecotype_counts = defaultdict(int)
    # Count genotypes in geno_drift_resh
    for i in geno_new_rpos:
        sample = i.tobytes()
        ecotype = ecotype_geno_mapper.get(sample, 'other')
        ecotype_counts[ecotype] += 1
    name = pop_name
    ecotype_countsdf = pd.DataFrame(list(ecotype_counts.items()), columns=['ecotype', name])
    ecotype_countsdf['ecotype'] = ecotype_countsdf['ecotype'].str.split('_').str[0]
    return ecotype_countsdf

nonhet_pos = np.array(pd.read_csv(nonhet_pos))
vcf_og = allel.read_vcf(og_vcf_offset, fields=['calldata/GT', 'variants/POS', 'samples'])
geno_og = vcf_og['calldata/GT']
pos_og = vcf_og['variants/POS']
samples = vcf_og['samples']

#ecotypes_grenenet = pd.read_csv(ecotypes_grenenet, dtype=object)
#ecotypes_grenenet.columns= ['ecotype']
#ecotypes_grenenet = pd.concat([ecotypes_grenenet, pd.DataFrame(data = {'ecotype': ['other']}, index=[231])],axis=0)

def process_vcf(i):
    name = i.split('/')[-2] + '_' + i.split('/')[-1][0:5]
    if os.path.exists(i) and os.path.getsize(i) <= 1:
        pass
    elif os.path.exists(i) and os.path.getsize(i) > 1:
        ## import each of teh vcf 
        vcf_new = allel.read_vcf(i, fields = ['calldata/GT', 'variants/POS'])
        ##extract the posisions and the geno array
        pos_new = vcf_new['variants/POS']
        geno_new = vcf_new['calldata/GT']
        ## for each of them create the ecotype geno mapper, depending on the positions that made it 
        geno_og_rpos, geno_new_rpos = filtering_pos(nonhet_pos, pos_new, geno_og, geno_new)
        ecotype_geno_mapper = get_ecotype_geno_mapper(geno_og_rpos)
        ecotype_countsdf = get_ecotype_counts(geno_new_rpos, name, ecotype_geno_mapper)
        print(ecotype_countsdf)
        ## merge with previous 
        #ecotypes_grenenet = ecotypes_grenenet.merge(ecotype_countsdf, how='left', on ='ecotype')
    return ecotype_countsdf

if __name__ == "__main__":
    # Prepare arguments for each dataset generation
    tasks = [output_vcf_offset]

    # Number of CPU cores to use (adjust as needed)
    num_cores = 20

    # Create a pool of worker processes
    with multiprocessing.Pool(processes=num_cores) as pool:
        # Use the pool to parallelize the processing of tasks
        results = pool.map(generate_dataset, tasks)

print(results)
ecotypes_grenenet = pd.concat(results)

#ecotypes_grenenet = ecotypes_grenenet.merge(ecotype_countsdf, how='left', on ='ecotype')
ecotypes_grenenet.to_csv(ecotype_counts)
