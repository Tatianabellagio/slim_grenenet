fam_file="${snakemake_input[fam_file]}"
bed_file="${snakemake_input[bed_file]}"
bim_file="${snakemake_input[bim_file]}"
kinship="${snakemake_input[kinship]}"
base_directory=$(pwd)
output="${snakemake_output[output_gwas]}"

## extract the directory in which the fam file is 
fam_dir=$(dirname "$fam_file")
#move to that directoryd 
cd "$fam_dir"
echo "$fam_dir"
## create the hardlinks to bed bim and kinship matrix
echo "${base_directory}/${bed_file}"

ln -sf "${base_directory}/${kinship}" kinship.cXX.txt
ln -sf "${base_directory}/${bim_file}" geno.bim
ln -sf "${base_directory}/${bed_file}" geno.bed

# Count the number of columns in the file
num_columns=$(awk -F ' ' 'NR==1 {print NF; exit}' "geno.fam")

# Check if the number of columns is not 6
if [ "$num_columns" -ne 6 ]; then
    echo "fam file wrong dimensions"
    exit 1
else
    # Check if all values in the 6th column are 0
    if awk -F ' ' '{print $6}' "geno.fam" | grep -qE '^[^0]|0[^.].*$'; then
        gemma -bfile geno -lmm -k kinship.cXX.txt -o results_nmaf
    else
        mkdir output
        echo "All individuals died" > "output/results_nmaf.assoc.txt"
    fi
fi

