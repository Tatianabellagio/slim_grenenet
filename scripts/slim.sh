#!/bin/bash
echo 'hola'

tree_seq_causalloci="$1"
optima="$2"
selection="$3"
heritability_state="$4"
output_tree_gen4="$5"
#output_tree_gen10="$6"
output_pop_size_early="$6"
output_pop_size_late="$7"
output_va="$8"
output_vpheno="${9}"
output_mfitness="${10}"
output_vfitness="${11}"
output_mean_pheno="${12}"
output_sd_pheno="${13}"
output_st_phenom="${14}"
output_st_phenov="${15}"

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
    -d "output_vpheno='$output_vpheno'" \
    -d "output_mfitness='$output_mfitness'" \
    -d "output_vfitness='$output_vfitness'" \
    -d "output_mean_pheno='$output_mean_pheno'" \
    -d "output_sd_pheno='$output_sd_pheno'" \
    -d "output_st_phenom='$output_st_phenom'" \
    -d "output_st_phenov='$output_st_phenov'" \
    scripts/arabidopsis_evolve_treeseq.slim 


