##### phase 1 #######################################
## filter ch5 and grenenet ecotypes from the big file 

bcftools view 1001_genomes_snps_missing0.8_merged_imputed_biallelic_named.vcf.gz --regions 5 -o chr5.vcf.gz
bcftools view -S ecotypes_grenenet.txt chr5.vcf.gz -o chr5_grenenet.vcf

## for testing only sample 0.01
bcftools view chr5_grenenet.vcf | vcfrandomsample -r 0.01 > chr5_grenenet_subset.vcf

##### phase 2 #######################################
######################################################
####### run python script ############################
######################################################


##### phase 3 #######################################
############## annotate vcf with sc

## first
bgzip selection_coef_chr5_subset.bed 
bgzip chr5_grenenet_subset.vcf 
bcftools tabix selection_coef_chr5_subset.bed.gz 
bcftools tabix -f chr5_grenenet_subset.vcf.gz 

## and then annotate
bcftools annotate \
  -a selection_coef_chr5_subset.bed.gz \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  chr5_grenenet_subset.vcf.gz \
  -o chr5_grenenet_subset_ann.vcf


##### phase 4 #######################################
###### amplify number of ecotypes #########################################################
bcftools query -l chr5_grenenet_subset_ann.vcf | wc -l ## now only 225

bgzip chr5_grenenet_subset_ann.vcf
bcftools tabix -f chr5_grenenet_subset_ann.vcf.gz 

bcftools merge --threads 4 --force-samples chr5_grenenet_subset_ann.vcf.gz chr5_grenenet_subset_ann.vcf.gz chr5_grenenet_subset_ann.vcf.gz chr5_grenenet_subset_ann.vcf.gz -o chr5_grenenet_subset_ann1.vcf.gz
bcftools tabix -f chr5_grenenet_subset_ann1.vcf.gz 

#bcftools query -l chr5_grenenet_subset_ann1.vcf.gz | wc -l ## 900

bcftools merge --threads 4 --force-samples chr5_grenenet_subset_ann1.vcf.gz chr5_grenenet_subset_ann1.vcf.gz chr5_grenenet_subset_ann1.vcf.gz chr5_grenenet_subset_ann1.vcf.gz -o chr5_grenenet_subset_ann2.vcf.gz
bcftools tabix -f chr5_grenenet_subset_ann2.vcf.gz 

#bcftools query -l chr5_grenenet_subset_ann2.vcf.gz | wc -l ## 3600

bcftools merge --threads 4 --force-samples chr5_grenenet_subset_ann2.vcf.gz chr5_grenenet_subset_ann2.vcf.gz chr5_grenenet_subset_ann2.vcf.gz chr5_grenenet_subset_ann2.vcf.gz -o chr5_grenenet_subset_ann3.vcf.gz
bcftools tabix -f chr5_grenenet_subset_ann3.vcf.gz 

#bcftools query -l chr5_grenenet_subset_ann3.vcf.gz | wc -l ## 14400

bcftools merge --threads 4 --force-samples chr5_grenenet_subset_ann3.vcf.gz chr5_grenenet_subset_ann3.vcf.gz -o chr5_grenenet_subset_ann4.vcf.gz
bcftools tabix -f chr5_grenenet_subset_ann4.vcf.gz 

#bcftools query -l chr5_grenenet_subset_ann4.vcf.gz | wc -l ## 28800
# so this will give 2400 ind per subpop 

##### phase 5 #######################################
## run slim code 