
echo awk 'BEGIN { rows=0; cols=0; } { rows++; cols=NF; } END { print rows, cols; }' ../data/greneNet_final_v1.1.recode.fam

#if awk -F ' ' 'NF==6 {exit 1}' "../data/greneNet_final_v1.1.recode.fam"; then
if awk -F ' ' 'NF!=6 {flag=1; exit} END {if (flag) exit 1}' "../data/greneNet_final_v1.1.recode.fam"; then

    ## Run GWAS with Gemma gemma -bfile geno -lmm -k kinship.cXX.txt -o results -maf 0.00001
    echo run_gemma
else
    ## Create a file indicating all individuals died
    echo dont run gemma 
fi



#if awk -F ' ' '{print $6}' "geno.fam"| grep -vqE '^0(\.0+)?$'; then
    