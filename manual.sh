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

#number of positions: 
bcftools view -H input.vcf.gz | wc -l  

#number of samples
bcftools query -l inv2.vcf | wc -l

#name of samples: 
bcftools query -l input.vcf

#name of samples save it in a txt 
bcftools query -l input.vcf -o sample_names.txt 

#filter samples 
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


____________________

March 27th 
# ok, simulations are workign now 
# and I was able to ouput vcf files after 7 generations, where selection started in generation 4 

# now what to doy with the vcf files?
# well, the whole idea is to simulate the grenenet experiment and for that, we have to simulate the fact that there was certain depth in the reads which was not total 

# for that we are going to use a tool called 
art_illumina

# I should solve this later but for now, 
export PATH=$PATH:/usr/local/bin/art_bin_MountRainier
art_illumina
## to make it run 

# ok I just realized that als these softwares are to simulate single individual sequencing, so what about pool seq?
# but meixi had a great idea of using a normal software and pretend all the individuals in a vcf file are actually different chromosomes
# so the softwares should accept multiple fasta files as inputs 

## so first step would be to convert my vcf file into multiple fasta files 

_____ march 28

fisrt think about how would you convert vcf file to a fasta file 

So, first of all you should know that for each chromosomes in a diploid individual there is tipically
2 fasta files
since each fasta file corresponds to 1 molecule of dna 
fasta files do not include information about hetercygocity is a plain chain of nucleotides 

But, becuase we are supposidly workign with only homocgygot, we are not supposed to see any heterocygocity 

ok i managed to write a python code to create fasta files from vcf files 

now using art-illumina to generate reads from fasta files: 


## this si art illumina basic functionality:
art_illumina -ss HS25 -i <input.fasta> -l <read_length> -f <coverage> -o <output_prefix>

art_illumina -ss HS25 -i genome.fasta -l 100 -f 30 -o reads
# ART will generate two output files with the extensions ".1.fq" and ".2.fq" for paired-end reads

ok so i have no idea for each run which platform we used but i found a pdf inside the fodler of the reads
from novogene , so i will just run it with the example 

art_illumina -ss HS25 -i i0.fasta -l 100 -f 30 -p -m 200 -s 10 -o reads


##ok now i have the reads, what i need to do is run the whole grenenet pipeline 
also i should exlore non hw equilibrium since im ending with too many individuals which is not quite real



###


____ march 29th 

i have some problems in slim, because i want to run a non wright fisher model to be able to 
not have a fixed population size, so plants can die based on their fitness, 
but the nonwrite fisher does not allow me to have addsubpop

I can think fo 2 options, the first one is to directly add many repetitions of individuals in the vcf
but that seems like a lot 

the second one is to:
1. create the pop of 224 individuals with the mut but keeping the mutation effect on 0 
2. make it grow exponentially until a certain value
3. create the subpopulations 
4. make the individuals migrate from the big popualtion to the subpop 


____________________
march 30 

one fo the commment i got was that is better if i get the popualtion 0 directly from a vcf file 
since i will be adding randomness everytime i run the simulation if i do the expoential growth things


## so i will take my original vcf file 
chr1_grenenet_ecotypes_subset ## 225 samples
## and i will basically merge it with itself 4 times, so i will 4 times the number of samples 
## hint: concat is used for vertically pasting, merge is used for horizontal pasting 

## first i need the compressed format
bcftools merge --threads 4 --force-samples chr1_grenenet_ecotypes_subset.vcf.gz chr1_grenenet_ecotypes_subset.vcf.gz chr1_grenenet_ecotypes_subset.vcf.gz chr1_grenenet_ecotypes_subset.vcf.gz -o chr1_grenenet_multiple_ecotypes_subset.vcf.gz

## output now is 
chr1_grenenet_multiple_ecotypes_subset.vcf.gz

## checking 
bcftools query -l chr1_grenenet_ecotypes_subset.vcf.gz | wc -l
bcftools query -l chr1_grenenet_multiple_ecotypes_subset.vcf.gz | wc -l
## so now the starting population is 900


# now i have to regenerate the qtl 

## so now there are not going to be 3 types of segments in the genome, 
## there will rather be only one where neutral and constributing variants are taken from 

# Calculate the number of lines for each file
19625 ## is the total nubmer of positions 

196 # 19625*0.05   contributing 
19429 # neutral   

## get a list of all the positions:
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_multiple_ecotypes_subset.vcf.gz | grep -v '^#' > pos.txt

# Shuffle the lines randomly
shuf pos.txt > shuffled_pos.txt
# Extract the first file
head -n 196 shuffled_pos.txt > qtl1_contrib.txt
tail -n 19429 shuffled_pos.txt > qtl1_neutral.txt

##and now i can filter those based on shuffled positions 

bcftools view -T qtl1_contrib.txt chr1_grenenet_multiple_ecotypes_subset.vcf.gz -o qtl1_contrib.vcf
bcftools view -T qtl1_neutral.txt chr1_grenenet_multiple_ecotypes_subset.vcf.gz -o qtl1_neutral.vcf

#check 
bcftools view -H qtl1_neutral.vcf | wc -l
bcftools view -H qtl1_contrib.vcf | wc -l



____________________march 31

one other thing we discussed was how to define optimas 



_____ April 3
So, slim does something weird where it basically: if the position of the snp in the vcf file is 4, 
then in slim its gonna be 4, so the actual position -1. 
I think this si related to slim saving positions starting at 0 or something 

## so this is how to retrieve the effects and position of all the mutations 
ut =  sim.mutationsOfType(m2);
effects = ut.selectionCoeff; # so this will give the effects
cat(paste(ut.position, sep="\n")); # this will give the positions 


_______ april4 
ok, so right now, im getting selection coefficients for the contributing snps in the qtl from a normal,
but every time i run the script again, i get new selction coefficeints which is not replicable,
so i should give the selection coefficents with the vcf fiel so they are fixed. 

this is called annotation of vcf file, bascially adding a column to the info file 

So i created a test vcf file called test_addinfo.vcf (with only 8 positions) where im gonna try to annotate it in the INFO column with the slectio ncoeffcient

so 

## first
bgzip selection_coef_test.bed 
bgzip test_addinfo.vcf 
bcftools tabix selection_coef_test.bed.gz 
bcftools tabix test_addinfo.vcf.gz 

## and then annotate
bcftools annotate \
  -a selection_coef_test.bed.gz \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  test_addinfo.vcf.gz \
  -o test_addinfo_.vcf

## the code above works 

## so now i can pass the selection coefficients AS I WANT 

## so now i should generate the selection coefficients with python, based on teh nubmer fo variants, and then pass them to the vcf files 
## and then consume the vcf files

