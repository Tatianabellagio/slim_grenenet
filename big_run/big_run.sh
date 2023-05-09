#!/bin/bash

## load modules on cluster
#module load Conda/3.7
module load BCFtools/1.10.2
module load SLiM/4.0.1
module load HTSlib/1.10.2 

###activate env on cluster
#conda activate pipeline_slim
source /carnegie/binaries/centos7/conda/3.7/bin/activate pipeline_slim

####### run the pipeline ############################
######################################################

##in case it is compressed from other runs 
bgzip -d chr5_grenenet.vcf.gz 

# python pi_beta_optima_fasta.py pi beta
echo PHASE 1 based on beta and pi calculate optimum phenotypes based on selected snps and their effect sizes 
python pi_beta_calc.py 0.0001 5

##### phase 3 #######################################
############## annotate vcf with sc 
echo PHASE 2 annotate vcf with effect sizes 
## first
bgzip -f selection_coef_chr5.bed 
bgzip -f chr5_grenenet.vcf 
bcftools tabix -f selection_coef_chr5.bed.gz 
bcftools tabix -f chr5_grenenet.vcf.gz 
# -f force in case this is past the first run of the pipeline and these files already exist 

## and then annotate
bcftools annotate \
  -a selection_coef_chr5.bed.gz \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  chr5_grenenet.vcf.gz \
  -o chr5_grenenet_ann.vcf

#decompress so it can be usen by the python when run this script again 
bgzip -d chr5_grenenet.vcf.gz 

echo PHASE 4 amplify number of ecotypes to simulate starting seedmix
##### phase 4 #######################################
###### amplify number of ecotypes #########################################################
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' chr5_grenenet_ann.vcf ## now only 225
bgzip -f chr5_grenenet_ann.vcf
bcftools tabix -f chr5_grenenet_ann.vcf.gz 

bcftools merge --threads 4 --force-samples -O z chr5_grenenet_ann.vcf.gz chr5_grenenet_ann.vcf.gz chr5_grenenet_ann.vcf.gz chr5_grenenet_ann.vcf.gz -o chr5_grenenet_ann1.vcf.gz
bcftools tabix -f chr5_grenenet_ann1.vcf.gz 
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' chr5_grenenet_ann1.vcf  ## 900

bcftools merge --threads 4 --force-samples -O z chr5_grenenet_ann1.vcf.gz chr5_grenenet_ann1.vcf.gz chr5_grenenet_ann1.vcf.gz chr5_grenenet_ann1.vcf.gz -o chr5_grenenet_ann2.vcf.gz
bcftools tabix -f chr5_grenenet_ann2.vcf.gz 
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' chr5_grenenet_ann2.vcf
  ## 3600

bcftools merge --threads 4 --force-samples chr5_grenenet_ann2.vcf.gz chr5_grenenet_ann2.vcf.gz chr5_grenenet_ann2.vcf.gz -o chr5_grenenet_ann3.vcf
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' chr5_grenenet_ann3.vcf
   ## 10800
## 
awk '{if ($1 == "#CHROM"){print NF-9; exit}}' chr5_grenenet_ann3.vcf




## phase 4
#echo PHASE 4 run slim simulation  
./bash_slim.sh

##### i think i will have to do ld pruning for slim to run faster
