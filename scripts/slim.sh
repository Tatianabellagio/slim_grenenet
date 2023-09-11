tree_seq_causalloci="${snakemake_input[tree_seq_causalloci]}"
optima_values="${snakemake_input[optima_values]}"
variance_values="${snakemake_input[variance_values]}"
optima_index="${snakemake_params[optima_index]}"
selection="${snakemake_params[selection]}"
output_file="${snakemake_output[0]}"

# Map 'selection' to its numeric value using a case statement
case "$selection" in
  'strongsel')
    variance_index=1
    ;;
  'moderatesel')
    variance_index=2
    ;;
  'lowsel')
    variance_index=3
    ;;
  *)
    echo "Invalid selection"
    exit 1
    ;;
esac

echo "tree_seq_causalloci: $tree_seq_causalloci" #>> "${snakemake_log[0]}"
echo "optima_values: $optima_values" #>> "${snakemake_log[0]}"
echo "variance_values: $variance_values" #>> "${snakemake_log[0]}"
echo "selection: $selection" #>> "${snakemake_log[0]}"
echo "variance_index: $variance_index"
echo "optima_index: $optima_index" #>> "${snakemake_log[0]}"
echo ${snakemake_output[0]}
echo "${snakemake_output[0]}"

mkdir -p "$output_folder/optima_index$optima_index"
echo "$output_folder/optima_index$optima_index"

slim \
    -d "tree='$tree_seq_causalloci'" \
    -d "optima_index='$optima_index'" \
    -d "optima_file='$optima_values'" \
    -d "variance_index='$variance_index'" \
    -d "variance_values='$variance_values'" \
    -d "output_file='$output_file'" \
    scripts/arabidopsis_evolve_treeseq.slim 
    


