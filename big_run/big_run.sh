####### run python script ############################
######################################################
## input: chr5_grenenet.vcf

### pi is 0.0001 since length_chr5 = 1702174 * 0.0001 = 170
# python pi_beta_optima_fasta.py pi beta 

python pi_beta_optima_fasta.py 0.0001 5

##### phase 3 #######################################
############## annotate vcf with sc

## first
bgzip selection_coef_chr5.bed 
bgzip chr5_grenenet.vcf 
bcftools tabix selection_coef_chr5.bed.gz 
bcftools tabix -f chr5_grenenet.vcf.gz 

## and then annotate
bcftools annotate \
  -a selection_coef_chr5.bed.gz \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  chr5_grenenet.vcf.gz \
  -o chr5_grenenet_ann.vcf

##### phase 4 #######################################
###### amplify number of ecotypes #########################################################
bcftools query -l chr5_grenenet_ann.vcf | wc -l ## now only 225

bgzip chr5_grenenet_ann.vcf
bcftools tabix -f chr5_grenenet_ann.vcf.gz 

bcftools merge --threads 4 --force-samples chr5_grenenet_ann.vcf.gz chr5_grenenet_ann.vcf.gz chr5_grenenet_ann.vcf.gz chr5_grenenet_ann.vcf.gz -o chr5_grenenet_ann1.vcf.gz
bcftools tabix -f chr5_grenenet_ann1.vcf.gz 

#bcftools query -l chr5_grenenet_ann1.vcf.gz | wc -l ## 900

bcftools merge --threads 4 --force-samples chr5_grenenet_ann1.vcf.gz chr5_grenenet_ann1.vcf.gz chr5_grenenet_ann1.vcf.gz chr5_grenenet_ann1.vcf.gz -o chr5_grenenet_ann2.vcf.gz
bcftools tabix -f chr5_grenenet_ann2.vcf.gz 

#bcftools query -l chr5_grenenet_ann2.vcf.gz | wc -l ## 3600

bcftools merge --threads 4 --force-samples chr5_grenenet_ann2.vcf.gz chr5_grenenet_ann2.vcf.gz chr5_grenenet_ann2.vcf.gz chr5_grenenet_ann2.vcf.gz -o chr5_grenenet_ann3.vcf.gz
bcftools tabix -f chr5_grenenet_ann3.vcf.gz 

#bcftools query -l chr5_grenenet_ann3.vcf.gz | wc -l ## 14400

bcftools merge --threads 4 --force-samples chr5_grenenet_ann3.vcf.gz chr5_grenenet_ann3.vcf.gz -o chr5_grenenet_ann4.vcf.gz
bcftools tabix -f chr5_grenenet_ann4.vcf.gz 

#for slim use, decompress
gunzip chr5_grenenet_ann4.vcf.gz 

#bcftools query -l chr5_grenenet_ann4.vcf.gz | wc -l ## 28800
# so this will give 2400 ind per subpop 

##### phase 5 #######################################
### run slim code thorugh all the optima ##########
#add permisison to run 
#chmod +x bash_slim.sh

./bash_slim.sh

