March 8th 

I am gonna try and download from the 1001 genome proejct the vcf files for the ecotypes we used
in grennet. I am gonna use only chr 1 from 0 to 30427671 since that is the length 
that it seems to have basen on the fasta file that i got from subseting the fas from all the chromosomes
from tair 

Chr1:0-30427671

Ok im testing if i can do this with the api from the 1001 genomes 
I ran this:

https://tools.1001genomes.org/api/v1/vcfsubset/strains/9998,9999/regions/Chr1:1000..1010,Chr2:2000..2010/type/fullgenome/format/vcf

and got a vcf file: 
9998,9999_1_1000-1010,2_2000-2010_fullgenome.vcf

So i decided to use the request module from python to interact with the 1001g project api
I wrote a script called download_vcf_1001g.py were i take all the ecotypes from the grenenet project THAT CAME FROM THE 1001G (becuase
some of them are from other sources so i could not retrieve them), i got an error when i tried getting the whole chromosome 1 

something that did not work is the api thing because it only allows me to download a small number, so I am going to use bcf tools (clasica y confiable)

in the cluster 

/home/tbellagio/safedata/1001g/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz

this is the file with all the ecotypes and all the chrom for vcf 

this command here should be filtering the samples in chr2 and by the ones stated in the tsv file 
bcftools view -r 2 -S ecotypes_grenenet.tsv /home/tbellagio/safedata/1001g/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz -Oz -o grenenet_ecotypes.vcf.gz
thi

bcftools view -Oz -o output.vcf.gz /home/tbellagio/safedata/1001g/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff.gz

so, basically i have to many sites with no ALT or .|. in the ecotype, this might mean monomorphic site but i am not sure
i will use the vcf file lucas has been using for grenenet that is in caltech, i should ask for calthec access 


reminder!!
when you transfer a vcf.gz file always remember to trasfer it together with a tbi file
that is basically the idnex of a compress file


March 14 

______
Classic command in bcftools 

number of positions: 
bcftools view -H input.vcf.gz | wc -l  

name of samples: 
bcftools query -l input.vcf

name of samples save it in a txt 
bcftools query -l input.vcf -o sample_names.txt 

filter samples 
bcftools view -S samplelist.txt input.vcf -o output.vcf

## filter regions or positions
bcftools view input.vcf.gz --regions 1 -o output.vcf
bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf > output.vcf
bcftools view -T positions.txt input.vcf > output.vcf


sometimes bcftools is annoying with compressing decompressing 
bgzip file.vcf ## to compress it 
gunzip file.vcf.gz ## to decompress it 

sometimes problems with index happen in compressed files, to create a new tbi file (index)
tabix -p vcf chr2.vcf.gz 

______


## So basically the file in the cluster named 
1001_genomes_snps_missing0.8_merged_imputed_biallelic_named.vcf.gz
## has all the ecotypes from the 1001 genomes, so I am going to filter by our ecotypes 

## Mar 15th ---

## Last time i realized how chr position is fucked up and i did the filter with bcftools wrong. Now i will fo it again 


# ok so: right way to do it is 1) first filter by regions and then 2) filter by samples 
bcftools view 1001_genomes_snps_missing0.8_merged_imputed_biallelic_named.vcf.gz --regions 1 -o chr1.vcf.gz
bcftools view -S ecotypes_grenenet.txt chr1.vcf.gz -o chr1_grenenet_ecotypes.vcf

## tip: try using the -o from vcf tools rather than the -> since it is safer 

## So, now my important file is chr1_grenenet_ecotypes.vcf


____________________
##March 23 

## I tried to do nucleotide based simulations but the vcf file might be too big, so i might ened to trim it 

## So, the chr1_grenenet_ecotypes.vcf file has 1969055 positions 

## from the 1969055 i will only keep 1969055 * 0.01 

## So based on the internet the best toold for doing this is vcflib
## that contains a function called vcfrandomsample

bcftools view chr1_grenenet_ecotypes.vcf | vcfrandomsample -r 0.01 > chr1_grenenet_ecotypes_subset.vcf
bgzip chr1_grenenet_ecotypes_subset.vcf
## once a filtered the vcf file no w have only 19681 variants 

## now i have to think about how to assign 
## different selection coefficient to different variant in different environments 

## now im trying to subset a set of regions of our subset chr1 vcf file into 4 parts: 4920 


## ok this might be weird, but im gonna create 4 files including the position for the 4 parts of my chromosome


### something is bad here but var2 and qtl2 were overlapping 
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | head -n 4920 > pos_qtl1.txt
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | tail -n +4921 | head -n 4920 > pos_inv1.txt
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | tail -n +9841 | head -n 4920 > pos_qtl2.txt
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | tail -n 4920 > pos_inv2.txt



bcftools view -r chr1:1000-2000,chr2:3000-4000 input.vcf > output.vcf


bcftools view -T pos_qtl1.txt chr1_grenenet_ecotypes_subset.vcf.gz -o qtl1.vcf
bcftools view -T pos_inv1.txt chr1_grenenet_ecotypes_subset.vcf.gz -o inv1.vcf
bcftools view -T pos_qtl2.txt chr1_grenenet_ecotypes_subset.vcf.gz -o qtl2.vcf


num_vars=$(bcftools view -H qtl1.vcf | wc -l)
shuf -i 1-$num_vars | split -d -l $((num_vars/2)) - vcf_indices_

## to now take 0.8 and 0.2 of the variants in the vcf file 
## we need to take 0.2 and 0.8 of the posiitons in the txt files: 


# Count the total number of lines in the input file
total=$(wc -l < pos_qtl1.txt)  ## 4920

# Calculate the number of lines for each file
n1=$(echo "scale=0; $total * 0.2" | bc) 
echo "scale=0; 4920 * 0.2" | bc  ## 984.0

n2=$(echo "scale=0; $total * 0.8" | bc)
echo "scale=0; 4920 * 0.8" | bc  ## 3936

# Shuffle the lines randomly
shuf pos_qtl1.txt > shuffled_pos_qtl1.txt
shuf pos_qtl2.txt > shuffled_pos_qtl2.txt
# Extract the first file
head -n 984 shuffled_pos_qtl1.txt > qtl1_contrib.txt
head -n 984 shuffled_pos_qtl2.txt > qtl2_contrib.txt

# Extract the second file
tail -n 3936 shuffled_pos_qtl1.txt > qtl1_neutral.txt
tail -n 3936 shuffled_pos_qtl2.txt > qtl2_neutral.txt


##and now i can filter those based on shuffled positions 

bcftools view -T qtl1_contrib.txt qtl1.vcf -o qtl1_contrib.vcf
bcftools view -T qtl1_neutral.txt qtl1.vcf -o qtl1_neutral.vcf


bcftools view -T qtl2_contrib.txt qtl2.vcf -o qtl2_contrib.vcf
bcftools view -T qtl2_neutral.txt qtl2.vcf -o qtl2_neutral.vcf


## now i changed my mind and im only doing one qtl 


bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | head -n 4920 > pos_inv1.txt
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | tail -n +4921 | head -n 4920 > pos_qtl1.txt
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_ecotypes_subset.vcf.gz | grep -v '^#' | tail -n +9841 > pos_inv2.txt

bcftools view -T pos_inv1.txt chr1_grenenet_ecotypes_subset.vcf.gz -o inv1.vcf
bcftools view -T pos_qtl1.txt chr1_grenenet_ecotypes_subset.vcf.gz -o qtl1.vcf
bcftools view -T pos_inv2.txt chr1_grenenet_ecotypes_subset.vcf.gz -o inv2.vcf


## to now take 0.8 and 0.2 of the variants in the vcf file 
## we need to take 0.2 and 0.8 of the posiitons in the txt files: 

# Count the total number of lines in the input file
## 4920

# Calculate the number of lines for each file
n1=$(echo "scale=0; $total * 0.2" | bc) 
echo "scale=0; 4920 * 0.2" | bc  ## 492

n2=$(echo "scale=0; $total * 0.8" | bc)
echo "scale=0; 4920 * 0.8" | bc  ## 4428

# Shuffle the lines randomly
shuf pos_qtl1.txt > shuffled_pos_qtl1.txt
# Extract the first file
head -n 492 shuffled_pos_qtl1.txt > qtl1_contrib.txt
tail -n 4428 shuffled_pos_qtl1.txt > qtl1_neutral.txt

##and now i can filter those based on shuffled positions 

bcftools view -T qtl1_contrib.txt qtl1.vcf -o qtl1_contrib.vcf
bcftools view -T qtl1_neutral.txt qtl1.vcf -o qtl1_neutral.vcf

