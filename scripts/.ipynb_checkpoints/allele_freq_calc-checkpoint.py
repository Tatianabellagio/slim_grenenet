import pandas as pd
import allel

output_vcf_fixpos = snakemake.input['output_vcf_fixpos'] 
output_allele_freq = snakemake.output['allele_freq'] 

def extract_allele_freq(samples, geno_array, pos, name):
    total_alleles = len(samples) * 2 
    alt_freq = geno_array.sum(axis=2).sum(axis=1) / total_alleles
    alt_allele_freq = pd.DataFrame(data = {'pos': pos, name: alt_freq})
    return alt_allele_freq

for i in output_vcf_fixpos:
    name = i.split('/')[-2] + '_' + i.split('/')[-1][0:5]
    vcf = allel.read_vcf(i)
    samples = vcf['samples']
    geno_array = vcf['calldata/GT']
    pos = vcf['variants/POS']
    alt_allele_freq = extract_allele_freq(samples,geno_array, pos, name)
    pos_vcf_og = pos_vcf_og.merge(alt_allele_freq, on ='pos', how='outer')

print(output_allele_freq)
pos_vcf_og.to_csv(output_allele_freq)