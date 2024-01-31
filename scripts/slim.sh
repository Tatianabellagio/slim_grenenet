#!/usr/bin/env bash
echo hola
tree_seq_causalloci="${snakemake_input[tree_seq_causalloci]}"
echo hola2
echo ${snakemake_params[lelo]}
echo ${snakemake_params[optima]}

optima=${snakemake_params[optima]}
selection=${snakemake_params[selection]}
heritability_state=${snakemake_params[heritability]}
output_tree_gen4="${snakemake_output[output_tree_gen4]}"
#output_tree_gen10="$6"
output_pop_size_early="${snakemake_output[output_pop_size_early]}"
output_pop_size_late="${snakemake_output[output_pop_size_late]}"
output_va="${snakemake_output[output_va]}"
output_mfitness="${snakemake_output[output_mfitness]}"
output_vfitness="${snakemake_output[output_vfitness]}"
output_mpheno="${snakemake_output[output_mpheno]}"
output_vpheno="${snakemake_output[output_vpheno]}"

# Map 'selection' to its numeric value using a case statement
case "$selection" in
  'exstrongsel')
    variance=0.001
    ;;
  'estrongsel')
    variance=0.01
    ;;
  'vstrongsel')
    variance=0.05
    ;;
  'strongsel')
    variance=0.1
    ;;
  'strongmod')
    variance=0.5
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

case "$heritability_state" in
  '1')
    h2=0.1
    ;;
  '2')
    h2=0.3
    ;;
  '3')
    h2=0.5
    ;;
  '4')
    h2=0.7
    ;;
  '5')
    h2=0.9
    ;;
  *)
    echo "Invalid selection"
    exit 1
    ;;
esac

echo $selection
echo $variance

slim \
    -d "tree='$tree_seq_causalloci'" \
    -d "h2='$h2'" \
    -d "optima='$optima'" \
    -d "variance='$variance'" \
    -d "output_tree_gen4='$output_tree_gen4'" \
    -d "output_tree_gen10='$output_tree_gen10'" \
    -d "output_pop_size_early='$output_pop_size_early'" \
    -d "output_pop_size_late='$output_pop_size_late'" \
    -d "output_va='$output_va'" \
    -d "output_mfitness='$output_mfitness'" \
    -d "output_vfitness='$output_vfitness'" \
    -d "output_mpheno='$output_mpheno'" \
    -d "output_vpheno='$output_vpheno'" \
    scripts/arabidopsis_evolve_treeseq.slim 


