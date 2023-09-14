import pandas as pd
import allel
import os
import numpy as np

output_vcf_fixpos = snakemake.input['output_vcf_fixpos'] 
pos_vcf_og = snakemake.input['pos_vcf_og'] 
output_allele_freq = snakemake.output['allele_freq'] 

pos_vcf_og = pd.read_csv(pos_vcf_og)

def extract_allele_freq(samples, geno_array, pos, name):
    total_alleles = len(samples) * 2 
    alt_freq = geno_array.sum(axis=2).sum(axis=1) / total_alleles
    alt_allele_freq = pd.DataFrame(data = {'pos': pos, name: alt_freq})
    return alt_allele_freq

for i in output_vcf_fixpos:
    print(i)
    name = i.split('/')[-2] + '_' + i.split('/')[-1][0:5]
    if os.path.exists(i) and os.path.getsize(i) <= 1:
        print('empty_vcf')
        pos_vcf_og[name] = np.nan
    elif os.path.exists(i) and os.path.getsize(i) > 1:
        vcf = allel.read_vcf(i)
        samples = vcf['samples']
        geno_array = vcf['calldata/GT']
        pos = vcf['variants/POS']
        alt_allele_freq = extract_allele_freq(samples,geno_array, pos, name)
        pos_vcf_og = pos_vcf_og.merge(alt_allele_freq, on ='pos', how='outer')
pos_vcf_og.to_csv(output_allele_freq)

print(output_allele_freq)
pos_vcf_og.to_csv(output_allele_freq)