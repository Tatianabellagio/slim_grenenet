vcf_file="${snakemake_input[vcf_file]}"
maf="${snakemake_input[maf]}"
corr_ld="${snakemake_input[corr_ld]}"


## first compress 
plink --vcf "{$vcf_file}" --make-bed --out data/chr5_grenenet 2> "${snakemake_log[0]}"

#ok im gonna filter by maf before

plink --bfile data/chr5_grenenet --maf "${maf}" --make-bed --out data/chr5_grenenet_filteredmaf 2>> "${snakemake_log[0]}"
#1135068 variants removed due to minor allele threshold(s)
#kept 567106


# ok so how to do ld pruning 
# this si a comman to kinda check 
plink --bfile data/chr5_grenenet_filteredmaf --indep-pairwise 50 5 "${corr_ld}" --out data/chr5_grenenet_filteredmaf 2>> "${snakemake_log[0]}"

