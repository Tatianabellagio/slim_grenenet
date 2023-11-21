base_directory=$(pwd)
pheno="${snakemake_input[pheno_hapfm]}"
haplotypeDM="${snakemake_input[haplotypeDM]}"
covariates="${snakemake_input[covariates]}"


## extract the directory in which the pheno file is 
pheno_dir=$(dirname "$pheno")
#move to that directoryd 
cd "$pheno_dir"
echo "$pheno_dir"
## create the hardlinks to bed bim and kinship matrix
echo "${base_directory}/${haplotypeDM}"

ln -sf "${base_directory}/${haplotypeDM}" haplotypeDM
ln -sf "${base_directory}/${covariates}" covariates

# Check if all values in the 6th column are 0
if awk -F ' ' '{print $1}' "geno.fam" | grep -qE '^[^0]|0[^.].*$'; then
    python3 ${base_directory}/scripts/hapfm/HapFM_mapping.py -i haplotypeDM -y pheno.txt -o cov3pheno1000 -c covariates
else
    if [ -d "output" ]; then
        # If the output directory exists, just create the file
        echo "All individuals died" > "cov3pheno1000_block_pip.txt"
    else
        # If the output directory does not exist, create it first, then create the file
        mkdir output
        echo "All individuals died" > "cov3pheno1000_block_pip.txt"
    fi
fi

