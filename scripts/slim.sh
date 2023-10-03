tree_seq_causalloci="${snakemake_input[tree_seq_causalloci]}"
optima="${snakemake_params[optima]}"
selection="${snakemake_params[selection]}"
heritability_state="${snakemake_params[heritability]}"
h2="${snakemake_params["$heritability_state"]}"
output_tree_gen4="${snakemake_output[output_tree_gen4]}"
output_tree_gen10="${snakemake_output[output_tree_gen10]}"
output_pop_size="${snakemake_output[output_pop_size]}"
output_va="${snakemake_output[output_va]}"
output_vpheno="${snakemake_output[output_vpheno]}"
output_mfitness="${snakemake_output[output_mfitness]}"
output_vfitness="${snakemake_output[output_vfitness]}"

# Map 'selection' to its numeric value using a case statement
case "$selection" in
  'strongsel')
    variance=0.1
    ;;
  'moderatesel')
    variance=1
    ;;
  'lowsel')
    variance=3
    ;;
  *)
    echo "Invalid selection"
    exit 1
    ;;
esac

echo $optima

slim \
    -d "tree='$tree_seq_causalloci'" \
    -d "h2='$h2'" \
    -d "optima='$optima'" \
    -d "variance='$variance'" \
    -d "output_tree_gen4='$output_tree_gen4'" \
    -d "output_tree_gen10='$output_tree_gen10'" \
    -d "output_pop_size='$output_pop_size'" \
    -d "output_va='$output_va'" \
    -d "output_vpheno='$output_vpheno'" \
    -d "output_mfitness='$output_mfitness'" \
    -d "output_vfitness='$output_vfitness'" \
    scripts/arabidopsis_evolve_treeseq.slim 


