#!/bin/bash
echo 'hola'

tree_seq_causalloci="$1"
echo $tree_seq_causalloci
optima="$2"
selection="$3"
heritability_state="$4"
output_tree_gen4="$5"
output_tree_gen10="$6"
output_pop_size="$7"
output_va="$8"
output_vpheno="$9"
output_mfitness="$10"
output_vfitness="$11"

echo $output_mfitness
echo $output_vfitness

echo $heritability_state
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

echo 'hola2'

case "$heritability_state" in
  'lowh')
    h2=0.1
    ;;
  'mediumh')
    h2=0.5
    ;;
  'highh')
    h2=0.9
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


