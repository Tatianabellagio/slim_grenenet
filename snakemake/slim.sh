#!/opt/homebrew/bin/bash

# this is an auxiliary script to call slim, since all the parameters have to be included and there are some troubles with double quotes 

while getopts f:v:s:o:p:i: flag
do
    case "${flag}" in
        f) fasta=${OPTARG};;
        v) vcf=${OPTARG};;
        s) optima_slim=${OPTARG};;
        o) output_folder=$(echo "${OPTARG}" | cut -d'/' -f1);;
        p) optima=${OPTARG};;
        i) initial_pop=${OPTARG};;
    esac
done
echo "fasta: $fasta";
echo "vcf: $vcf";
echo "optima_slim: $optima_slim";
echo "output_folder: $output_folder";
echo "optima: $optima";
echo "initial_pop: $initial_pop";

slim -d "ref_fasta='$fasta'" -d "main_vcf='$vcf'" -d "optima_file='$optima_slim'" -d "optima='$optima'" -d "initial_pop='$initial_pop'" -d "output_folder='$output_folder'" arabidopsis_evolve.slim

