import pandas as pd

fam_file = snakemake.input['fam_gwas'] 
pheno_hapfm = snakemake.output['pheno_hapfm'] 

pd.read_csv(fam_gwas,sep=' ')

geno_fam = pd.read_csv('geno.fam',sep=' ',header=None)[5]
geno_fam = geno_fam * 1000
geno_fam.to_csv(pheno_hapfm,header=None, index=None)
