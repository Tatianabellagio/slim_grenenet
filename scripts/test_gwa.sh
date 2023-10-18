## Check if all individuals have the phenotype value as 0
if awk -F ' ' '{print $6}' "../data/greneNet_final_v1.1.recode.fam"| grep -vqE '^0(\.0+)?$'; then
#if awk -F ' ' '{print $6}' "../data/greneNet_final_v1.1.recode.fam" | grep -qE -v '^0(\.0+)?$'; then
    ## Run GWAS with Gemma gemma -bfile geno -lmm -k kinship.cXX.txt -o results -maf 0.00001
    echo run_gemma
else
    ## Create a file indicating all individuals died
    echo dont run gemma 
fi
 