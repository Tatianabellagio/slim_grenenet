
allele_counts_table="${snakemake_input[allele_counts_table]}"
prune_in="${snakemake_input[snps_to_keep]}"
allele_counts_table_filtered="${snakemake_output[allele_counts_table_filtered]}"

awk -F ',' 'BEGIN{OFS=","} NR==FNR{values[substr($1, 3)]; next} {gsub(/^5_/, "", $1); if ($1 in values) print}' "${prune_in}" "${allele_counts_table}" > "${allele_counts_table_filtered}" 2> "${snakemake_log[0]}"
