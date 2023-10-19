fam_file="${snakemake_input[fam_file]}"
bed_file="${snakemake_input[bed_file]}"
bim_file="${snakemake_input[bim_file]}"
kinship="${snakemake_input[kinship]}"
wd="${snakemake_param[wd]}"

output="${snakemake_output[output_gwas]}"
## extract the directory in which the fam file is 
fam_dir=$(dirname "$fam_file")

#move to that directoryd 
cd "$fam_dir"
echo "$fam_dir"
## create the hardlinks to bed bim and kinship matrix
echo $wd
echo "../../../../../../$bed_file"
ln -sf "../../../../../../$kinship" kinship.cXX.txt
ln -sf "../../../../../../$bim_file" geno.bim
ln -sf "../../../../../../$bed_file" geno.bed

## adn now that all the files are in the same folder, run gemma 
awk 'BEGIN { rows=0; cols=0; } { rows++; cols=NF; } END { print rows, cols; }' geno.fam


#if awk -F ' ' 'NF==6 {exit 1}' "../data/greneNet_final_v1.1.recode.fam"; then
if awk -F ' ' 'NF!=6 {flag=1; exit} END {if (flag) exit 1}' "geno.fam"; then
    ## Run GWAS with Gemma gemma -bfile geno -lmm -k kinship.cXX.txt -o results -maf 0.00001
    gemma -bfile geno -lmm -k kinship.cXX.txt -o results -maf 0.00001
else
    ## Create a file indicating all individuals died
    mkdir output
    echo "All individuals died" > "output/results.assoc.txt"
fi
