fasta="${snakemake_input[fasta]}"
vcf="${snakemake_input[vcf]}"
optima_slim="${snakemake_input[optima_slim]}"
output_folder=$(echo "${snakemake_output[0]}" | cut -d'/' -f1-2)
optima="${snakemake_params[optima]}"
initial_pop="${snakemake_params[initial_pop]}"

echo "fasta: $fasta" > "${snakemake_log[0]}"
echo "vcf: $vcf" >> "${snakemake_log[0]}"
echo "optima_slim: $optima_slim" >> "${snakemake_log[0]}"
echo "output_folder: $output_folder" >> "${snakemake_log[0]}"
echo "optima: $optima" >> "${snakemake_log[0]}"
echo "initial_pop: $initial_pop" >> "${snakemake_log[0]}"


slim \
    -d "ref_fasta='$fasta'" \
    -d "main_vcf='$vcf'" \
    -d "optima_file='$optima_slim'" \
    -d "output_folder='$output_folder'" \
    -d "optima='$optima'" \
    -d "initial_pop='$initial_pop'" \
    scripts/arabidopsis_evolve.slim >> "${snakemake_log[0]}" 2>> "${snakemake_log[0]}"

    #> output.txt 2> error.txt


