fam_file="${snakemake_input[fam_file]}"
bed_file="${snakemake_input[bed_file]}"
bim_file="${snakemake_input[bim_file]}"
kinship="${snakemake_input[kinship]}"

output="${snakemake_output[output_gwas]}"
## extract the directory in which the fam file is 
fam_dir=$(dirname "$fam_file")

#move to that directory 
cd "$fam_dir"
echo "$fam_dir"
## create the hardlinks to bed bim and kinship matrix

echo "../../../../../../$bed_file"
ln -sf "../../../../../../$kinship" kinship.cXX.txt
ln -sf "../../../../../../$bim_file" geno.bim
ln -sf "../../../../../../$bed_file" geno.bed

## adn now that all the files are in the same folder, run gemma 

## Check if all individuals have the phenotype value as 0
if awk -F ' ' '{print $6}' "geno.fam"| grep -vqE '^0(\.0+)?$'; then
    ## Run GWAS with Gemma gemma -bfile geno -lmm -k kinship.cXX.txt -o results -maf 0.00001
    gemma -bfile geno -lmm -k kinship.cXX.txt -o results -maf 0.00001
else
    ## Create a file indicating all individuals died
    mkdir output
    echo "All individuals died" > "output/results.assoc.txt"
fi
 