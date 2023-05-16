#ok for calculating ld with plink first we need to add ids to the pos for pruning porposes 


## dont need this part because it already has a name 
gzip subp.vcf
bcftools tabix subp.vcf.gz
bcftools annotate --set-id +'%CHROM\_%POS' subp3.vcf.gz > subp3.vcf

#second we need to transform my vcf to the plink format 

plink --vcf ../chr5_grenenet.vcf --make-bed --out chr5_grenenet

#now to calculate pairwise ld in between all my snps 

plink --bfile subp3 --r2 --matrix --out subp3_ld

#ok im gonna filter by maf before

plink --bfile chr5_grenenet --maf 0.01 --make-bed --out chr5_grenenet_filteredmaf


#now to calculate pairwise ld in between all my snps 

plink --bfile subp3_filteredmaf --r2 --matrix --out subp3_filteredmaf_ld

# ok so how to do ld pruning 
# this si a comman to kinda check 
plink --bfile chr5_grenenet_filteredmaf --indep-pairwise 50 5 0.6 --out chr5_grenenet_filteredmaf


awk 'END {print NR}' chr5_grenenet_filteredmaf.bim


# for doing ld pruning you usually need to specify 
# a sliding window size (50in this case) , and a step size (5 in this case)
# the sliding window size woudl be the distance between 2 basepairs
## being compared and the step size would be how many steps the
## sliding window will make each time 

# a threshold of 0.2 means that SNPs with a pairwise LD correlation 
# greater than or equal to 0.2 will be pruned

## this last comman will generate 2 ouputs
subp0_filteredmafandld.prune.in
subp0_filteredmafandld.prune.out

## ok this was giving me empty results because the id value is . 

### first we have to set an id by concatenating chr and pos, will do it at the begining 

## now this comand that actualyl makes the pruning aka filtering is 

plink --bfile subp3_filteredmaf --exclude subp3_filteredmafandld.prune.out --make-bed --out subp3_filteredmafandld

# and now use subp3_filteredmafandld to calcualte a matrix and cehck output 
plink --bfile subp3_filteredmafandld --r2 --matrix --out subp3_filteredmafandld





