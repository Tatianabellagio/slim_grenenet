#!/bin/bash
echo 'hola'

tree_seq_causalloci="$1"
optima="$2"
selection="$3"
heritability_state="$4"
output_tree_gen4="$5"
output_tree_gen10="$6"
output_pop_size_early="$7"
output_pop_size_late="$8"
output_va="$9"
output_vpheno="${10}"
output_mfitness="${11}"
output_vfitness="${12}"

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
    -d "output_vpheno='$output_vpheno'" \
    -d "output_mfitness='$output_mfitness'" \
    -d "output_vfitness='$output_vfitness'" \
    scripts/arabidopsis_evolve_treeseq.slim 


