import pandas as pd
import allel
import os
import numpy as np

output_vcf_fixpos = snakemake.input['output_vcf_fixpos'] 
pos_vcf_og = snakemake.input['pos_vcf_og'] 
output_allele_freq = snakemake.output['allele_freq'] 
output_allele_counts = snakemake.output['allele_counts'] 

pos_vcf_og_counts = pd.read_csv(pos_vcf_og)
pos_vcf_og_freq = pd.read_csv(pos_vcf_og)


def extract_allele_freq(samples, geno_array, pos, chrom, name):
    total_alleles = len(samples) * 2 
    alt_count = geno_array.sum(axis=2).sum(axis=1)
    alt_freq = alt_count / total_alleles

    chrom_pos = pd.Series(chrom.astype(str)) + '_' +  pd.Series(pos.astype(str))

    alt_allele_count = pd.DataFrame(data = {'chrom_pos': chrom_pos, name: alt_count})
    alt_allele_freq = pd.DataFrame(data = {'chrom_pos': chrom_pos, name: alt_freq})
    return alt_allele_count, alt_allele_freq

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
        chrom = vcf['variants/CHROM']
        alt_allele_count, alt_allele_freq = extract_allele_freq(samples,geno_array, pos, chrom, name)
        ## counts
        pos_vcf_og_counts = pos_vcf_og_counts.merge(alt_allele_count, on ='chrom_pos', how='outer')
        ## freq 
        pos_vcf_og_freq = pos_vcf_og_freq.merge(alt_allele_freq, on ='chrom_pos', how='outer')



print(len(pos_vcf_og))
pos_vcf_og_counts.to_csv(output_allele_counts)
pos_vcf_og_freq.to_csv(output_allele_freq)