import pandas as pd

fam_file = snakemake.input['fam_gwas'] 
pheno_hapfm = snakemake.output['pheno_hapfm'] 

geno_fam = pd.read_csv(fam_file,sep=' ',header=None)[5]
geno_fam = geno_fam * 1000
geno_fam.to_csv(pheno_hapfm,header=None, index=None)
