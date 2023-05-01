## im trying to do some forward-time simulations from vcf files from the 1001 genome project 
## for this im researching bcf tools to extract info from a a vcf.gz file which is basically 
## a compress vcf file. The original one from the 1001 genomes 
## is like 60 gb and the same one is in the cluster in:

/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz

## so i will use 

bcftools query -l 1001genomes_snp-short-indel_only_ACGTN.vcf.gz > sample_ids.txt
/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz

bcftools query -l /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz  > sample_ids_1001g.txt

## to get all the sample ids from the 1001 genomes 
## and then take a subsample of that with 

bcftools view -S list_of_samples.txt 1001genomes_snp-short-indel_only_ACGTN.vcf.gz > subset.vcf

## actually: 
bcftools view -S subsampled_ids_1001g.txt /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz > subset_1001g.vcf


## to try on run some simulations 

## also for this i had to load the bcf module
module load BCFtools/1.10.2

## once i got the samples_ids.txt file i nees to subsample it since i will try to run some mock 
## simulations and i dont want al the genomes 

## i will use bash script shuf for this, where you can indicate the percentage of lines you need 

shuf -n 10% sample_ids_1001g.txt > subsampled_ids_1001g.txt

## bueno 10 % parece ser muy poco since there are 1135 files actually 
## got the numbre of lines with:
wc -l sample_ids_1001g.txt


## also for playing with the simulations in an ordered way I created a new conda environemtn 
## one issue that i run into is actualyl being able to see my new environemtn from jupyer notebook
## i solved this by: 

 conda install ipykernel ## inside env
 python -m ipykernel install --user --name simulations --display-name "simulations"

 ## cosas que descubre el caminante al caminar 
 ## parece ser que los vcf files pueden tener 1 sample o multiple samples 

fileformat=VCFv4.3
FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1

## ese es el formato y de ahi se pueden ir agregando sample2 sample3 sample4 etc. 


### ok GRACIAS GPT 
bcftools view -s Sample1,Sample2 input.vcf | bcftools view -c 10 -Ov > output.vcf
## se supone que con este script puedo seleccion ciertos samples y la cantidad de snps 
bcftools view -s 430,9683 1001g_vcf_comp | bcftools view -c 10 -Ov > subs_samplesandsnps.vcf


## also i created a symlink beacuse the path is really long and annoying
## so i can access the datafile directly 
ln -s /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz 1001g_vcf_comp.vcf.gz


bcftools view -s 6909 /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz > col0.vcf



## end
im trying to run this 2 no succesfully 
 ### for getting the vcf file for col and subsample snps
 bcftools view -s 6909 /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz > col0.vcf
 
 ### for creating a little vcf file fo 2 samples and 10 snps 
 bcftools view -s 430,9683 /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz | bcftools view -c 10 -Ov > subs_samplesandsnps.vcf


#col0.vcf : vcf file of col0 accession to check it out 
#subs_samplesandsnps.vcf vcf of 2 samples to check out 

### SOMETHIGN TO BE AWARE OF!!!!!
## WHEN I WAS https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_chromosome_files
# in the tair page i was trying to filter the mithocondrial chromosome thinking it was chr1 

## ChrM~~  avoid



### ok now im trying to follow the manual from slim on how to generate simulations based on a 'real' population,
## since the inpu is going to be a table containing snps and insdividual from real data (from vcf files from the 1001 genome)

## the manual said that i need 2 files: the vcf file of all the indivduals and a fasta file from the reference genome 
## i will assume the reference is the one from tair 

## im just gona work with one part of the genome so i needed to subset that part 

# first based on the file a downloaded from tair im gona subset chromosome 1 
gunzip -c filename.fasta.gz | awk '/^>Chr1 /{p=1}p' > chr1.fasta
gunzip -c TAIR10_chr_all.fas | awk '/^>Chr1 /{p=1}p' > chr1.fasta

TAIR10_chr_all.fas
## you can do that by: 
## create a file called

# subset.bed 
Chr1    9001000    10001000 #this is a tab separeted file
## the nsomething that it seems that you need to do to filter a fasta file is create a file with extension .fai

samtools faidx genome.fasta.gz
## samtols can do that with the coman faidx 

#once you ahve that 
bedtools getfasta -fi genome.fasta.gz -bed subset.bed -fo subset.fasta

bedtools getfasta -fi chr1.fasta -bed subset.bed -fo reference.fasta

## now we got the subset.fasta that has 
#>Chr1:9001000-10001000 as header 


so i got my reference genome from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release%2FTAIR10_chromosome_files
and i will get my vcf files from 30 accessions from https://tools.1001genomes.org/vcfsubset/#select_loci


## i need to elarn how to run things in slim 



## getting position in a vcf file
grep -v '^#' variants_30ecotypes.vcf | awk '{print $2}' | sort -n | awk 'BEGIN{min=1000000000}{if($1<min)min=$1}END{print "min position:",min}' ; awk '{print $2}' variants_30ecotypes.vcf | sort -nr | awk 'BEGIN{max=0}{if($1>max)max=$1}END{print "max position:",max}'

##length of seq in fasta file 
awk '/^>/ {if (seqlen){print seqlen}; seqlen=0;next; } { seqlen += length($0)} END {print seqlen}' input.fasta

## chequeando esto me di cuenta que mi fasta file 
## y mi vcf file tienen diferente length 
# jeje 
1000000

1000001

## parece se rque todo estopaso porque al general el subset de fasta hay qe fijarse bien si el 
## subsetting se hace inlcuyendo the upper end o no 



#la verdad que no tengo idea que carajo de problema
#est apasadno en slim pero me sigeu tirando el 
#error de que tal posicion esta out fo range 

ERROR (Genome_Class::ExecuteMethod_readFromVCF): VCF file POS value 9000999 out of range.
#por lo que voy a eliminarla y ver que pasado

#como eliminar una dada posicion con vcf tools
bcftools view -e 'POS == 9000999' variants_30ecotypes.vcf -o variants_30ecotypes.vcf


### SOLUCION
## basicamente para poder correr esta simulacion en lguar de poner un dado segmento del fasta file y el mismo 
## segmento del vcf file, la forma de hacer es pasar el fasta file COMPLETO del cromosoma y el vcf file 
## de los pedazos que haya

## esto tiene sentido biologico, ya que la tasa de recomebinacion etc va a ser a nivel cromosoma
## y no solo en los segmentos que le paso 

## una vez pasado el fasta compelto de chromosome 1 y el vcf de una parte de varias accessions tuve otros errores 

## fasta
## slim solo acepta ciertas variantes en los fasta file que son basicamente A,C,G,T y no otras, algo 
## que se usa en genomica 

"M": ["A"],
"R": ["A"],
"W": ["A"],
"S": ["C"],
"Y": ["C"],
"K": ["G"],
"V": ["A"],
"H": ["A"],
"D": ["A"],
"B": ["C"],
"N": ["A"]

## si osi hay que remplazarlos para que slim los lea 

sed 's/Y/C/g' TAIR10_chr1.fas > TAIR10_chr1_edited.fas

# esto se puede hacer con sed uno por uno o con un python script para hacerlos todos a la vez 
## podria armar uno alguna vez

## VCF
## luego, otra cosa que slim no acepta es que en la columna REF y ALT existan otras cosas que no sean los 4 nucleotidos 
## a veces en estas columnas en ref y alt existen mas de un tipo de nucleotidos, calculo que indicando que cualqueira
## de los dos podria ser 

## para eso lo que yo hice fue filtrar solo columnas donde el REF y ALT tenian un unico valor 

bcftools view -i 'REF=="A" || REF=="C" || REF=="G" || REF=="T"' input.vcf > output.vcf ## lo mismo para el ALT 
bcftools view -i 'ALT=="A" || ALT=="C" || ALT=="G" || ALT=="T"' test1.vcf > test.vcf

# en realidad la forma correcta de hacerlo es tal vez quedarse con el primer elemento 

## otra cosa que no acepta slim partes del genoma basicamente sin coverage, que serian las celdas sonde under samples 
## veo esto:    ./.  , por loque por ahora elimine todas las aprtes del genoma donde no habia full coverage

bcftools view -v snps -e 'GT="0/0" || GT="./."' input.vcf > output.vcf


bcftools view -v snps -e 'GT="./."' test.vcf > subset.vcf

## la forma correcta de ahcer esto seria en realidad reemplazar esto por el REF allele?

## y leugo de todas estasc orreciones por fin si corrio la simulacion 


## also another script that i have been using a lot is: 
bedtools getfasta -fi TAIR10_chr_all.fas -bed subset.bed -fo test.fasta
## this scripts allows you filter a fasta/fas file based on a subset.bed file to get another fasta file 

### BE CAREFUL ###
# be careful with the different tools you use to trim fasta and vcf files 

## for example: 

# bedtools getfasta 
$ cat test.fa
>chr1
AAAAAAAACCCCCCCCCCCCCGCTACTGGGGGGGGGGGGGGGGGG

$ cat test.bed
chr1 5 10  

$ bedtools getfasta -fi test.fa -bed test.bed
>chr1:5-10
AAACC  ## as you can see it is not inclusing hte position 5 if you start counting on 1 
       ## it starrts counting on 0 

## on the other hand when downloading vcf files from the 1001 genome projects the lower and upper end are included 
Chr1:9001000..9002000

## also i used this python script a lot to compare the length of the fasta and vcf file 

import pysam

# Open FASTA and VCF files
fasta_file = pysam.FastaFile("referencevcf3.fasta")
vcf_file = pysam.VariantFile("test.vcf")

# Count number of positions in FASTA and VCF files
num_positions_fasta = sum([len(fasta_file[contig]) for contig in fasta_file.references])
num_positions_vcf = sum([1 for record in vcf_file])

print(num_positions_fasta)
print(num_positions_vcf)

## was not useful at the end since they did not have to have the same length, but might be useful for other things 


for (mut in sim.mutations[1]) {
    catn("Mutation ID: " + mut.id);
    catn("Location: " + mut.position);
    catn("Effect: " + mut.selectionCoeff);
}
'"'
#####

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


#sometimes bcftools is annoying with compressing decompressing 
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


## now i will make it 1800
bcftools merge --threads 4 --force-samples chr1_grenenet_multiple_ecotypes_subset.vcf.gz chr1_grenenet_multiple_ecotypes_subset.vcf.gz -o chr1_grenenet_multiple_ecotypes_subset2.vcf.gz
bcftools query -l chr1_grenenet_multiple_ecotypes_subset2.vcf.gz | wc -l


# now i have to regenerate the qtl 

## so now there are not going to be 3 types of segments in the genome, 
## there will rather be only one where neutral and constributing variants are taken from 

# Calculate the number of lines for each file
19625 ## is the total nubmer of positions 

196 # 19625*0.05   contributing 
19429 # neutral   

## get a list of all the positions:
bcftools query -f '%CHROM\t%POS\n' chr1_grenenet_multiple_ecotypes_subset2.vcf.gz | grep -v '^#' > pos.txt

# Shuffle the lines randomly
shuf pos.txt > shuffled_pos.txt
# Extract the first file
head -n 196 shuffled_pos.txt > qtl1_contrib.txt
tail -n 19429 shuffled_pos.txt > qtl1_neutral.txt

##and now i can filter those based on shuffled positions 

bcftools view -T qtl1_contrib.txt chr1_grenenet_multiple_ecotypes_subset2.vcf.gz -o qtl1_contrib.vcf
bcftools view -T qtl1_neutral.txt chr1_grenenet_multiple_ecotypes_subset2.vcf.gz -o qtl1_neutral.vcf

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
bgzip selection_coef_qtl1_contrib.bed 
bgzip qtl1_contrib.vcf 
bcftools tabix selection_coef_qtl1_contrib.bed.gz 
bcftools tabix -f qtl1_contrib.vcf.gz 

## and then annotate
bcftools annotate \
  -a selection_coef_qtl1_contrib.bed.gz \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  qtl1_contrib.vcf.gz \
  -o qtl1_contrib_sc.vcf

## the code above works 

## so now i can pass the selection coefficients AS I WANT 

## so now i should generate the selection coefficients with python, based on teh nubmer fo variants, and then pass them to the vcf files 
## and then consume the vcf files



_____ april 11
so, now i have my optima defined by a gradient of temperature:

so i took random ecotypes from a gradient of temperature they are local to 
and then calculated the optima based on this ecotypes (their prs)

now i ahve to decide how to scale up this population 

but first i need to solve the number of plants this can be detemriend by the cesus data 




i will increase the number of starting populations


## first i need the compressed format
bcftools merge --threads 4 --force-samples chr1_grenenet_ecotypes_subset.vcf.gz chr1_grenenet_ecotypes_subset.vcf.gz chr1_grenenet_ecotypes_subset.vcf.gz chr1_grenenet_ecotypes_subset.vcf.gz -o chr1_grenenet_multiple_ecotypes_subset.vcf.gz

## output now is 
chr1_grenenet_multiple_ecotypes_subset.vcf.gz

## checking 
bcftools query -l chr1_grenenet_ecotypes_subset.vcf.gz | wc -l
bcftools query -l chr1_grenenet_multiple_ecotypes_subset.vcf.gz | wc -l


figure out where mutation effect, fitness scaling etc are acting, poutput each of them at every stage of the simulation 

____ april 12th 
ok i finally figured out the bug in my slim code. basically, it is important to pay attention to what is the difference bettwen wf and nonwf models, 
and what each callback does 

----- now i have to think about rescaling and what i need from the ouput 

so lets start with bayenv 
i know that for bayenv i need allele freq 


##3 baypass run 

## so it seems from my previous tries, that i only need 2 files: 
geno_table	pool_size_file

g_baypass -npop 28  \
          -gfile geno_table_grenenet_fil  \
          -poolsizefile pool_size_grenenet_fil  \
          -d0yij 2.8  \  
          -outprefix filtrun.BP  \
          -npilot 50

art_illumina -ss HS25 -i vcf_slim/fasta_slim/subp0.fasta -l 100 -f 30 -p -m 200 -s 10 -o vcf_slim/fasta_slim/reads_slim/subp0.fasta_reads

art_illumina -ss HS25 \
             -i vcf_slim/fasta_slim/subp1.fasta \
             -l 150 \
             -f 7 \
             -p \
             -m 800 \
             -s 100 \
             -o vcf_slim/fasta_slim/reads_slim/subp1.fasta



ok, i have bypass in my computer and in the cluster
and i have bayenv in my computer

to run baypass remember that the comman is actually g_baypass

also remember that to run things on the cluster: 

### jobs in the cluster

# Memex/calc interactive session:
# create a tmux session 
tmux new -s session_name
#run this in a tmux session, so that the interactive session stays alive:
srun --partition DPB --time 0 --mem 20G --cpus-per-task 4 --nodes 1 --ntasks 1 --pty bash -i


### creating ther eads for baypass 
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

art_illumina -ss HS25 -i i0.fasta -l 150 -f 7 -p -m 300 -s 10 -o reads

art_illumina -ss HS25 # right now im assuming HS25 - HiSeq 2500
             -i i0.fasta #input 
             -l 150 #length of sequences
             -f 7 # depth 
             -p # paired reads
             -m 800 #mean length 
             -s 100 # 
             -o reads

ok so fragment length nowadays is about 700 or 800 hundred with a standars dev of 100 

#-m   --mflen
# the mean size of DNA/RNA fragments for paired-end simulations

__ 
so talking with meixi and xing they recommended goign with the default: either novaseq or hiseq and usually the read length is 150 and the coverage for our 
experiment is in between 5 to 10, and it would be a good idea to model 2 extreme scenarios where the coverage is very low or the vocerage is very high 


i created a sh file 
simulate_reads.sh

#add permisison to run 
chmod +x simulate_reads.sh

and to run the sh 
./simulate_reads.sh


## one thing that i just obbserved and waskinda suprised abotu is that the vcf files generated by slim are in the magnitude of 5 th snps, when the og vcf had like 19 th variants
## but this makes sense, because there was a lot of ecotype sorting, so noe the variants are more, and sites that are not variant are not included in the vcf file 


## ok lets run baypass

### jobs in the cluster

# Memex/calc interactive session:
# create a tmux session 
tmux new -s session_name
#run this in a tmux session, so that the interactive session stays alive:
srun --partition DPB --time 0 --mem 20G --cpus-per-task 16 --nodes 1 --ntasks 1 --pty bash -i

--cpus-per-task 16 since baypass can use 16 threads

g_baypass -npop 16  \
          -gfile geno_table_slim  \
          -poolsizefile pool_size_slim  \
          -d0yij 2.8  \
          -outprefix filtrun.BP  \
          -npilot 50 \
          -nthreads 16



          geno_table_slim

          pool_size_slim



###### installing slim with conda in the cluster

conda install -c conda-forge slim 

conda create -n slim 
conda activate slim 
conda install -c conda-forge slim


###### if you are in need of getting soemthing from the cluster: 
scp tbellagio@calc.dpb.carnegiescience.edu:safedata/ath_evo/grenephase1/data/worldclim_ecotypesdata.csv ~/Desktop/

## my slim script has 3 parameters optima, replicates, gen_to_run
slim -d optima=0 -d replicates=12 -d gen_to_run=3 arabidopsis_evolve.slim



conda install 


## conda environemnt problems 
conda create --name simulations_slim python=3.8.9
conda activate simulations_slim
conda create --name simulations_slim --clone simulations
conda env remove --name my_old_env
