#!/bin/bash

output_file="qtl2_contrib1_final2.vcf.gz"
file="qtl2_contrib.vcf.gz" 

command="bcftools merge --threads 5 --force-samples -Oz"

for i in {1..64}; do
    command+=" $file"
done

command+=" -o $output_file"

eval "$command"


# command=""
# for i in `seq 1 48` ; do

#     command="${command} file.vcf"

# done

# vcftools -other -options $command

# bgzip qtl2_contrib.vcf
# echo awk '{if ($1 == "#CHROM"){print NF-9; exit}}' qtl2_contrib.vcf.gz ## now only 225
# bcftools tabix -f qtl2_contrib.vcf.gz

# bcftools merge --threads 10 --force-samples -O z qtl2_contrib.vcf.gz qtl2_contrib.vcf.gz qtl2_contrib.vcf.gz qtl2_contrib.vcf.gz -o qtl2_contrib1.vcf.gz

