#!/usr/bin/bash:wq

# this is for running as local
# #!/opt/homebrew/bin/bash
echo $BASH_VERSION

module load BCFtools/1.10.2
module load HTSlib/1.10.2

## this didnt work 
#echo "${snakemake_input[bed_sc]}"
#echo "${snakemake_input[bed_sc]}.gz"
#echo "${snakemake_input[vcf_file]}.gz"

vcf_file="$1"
bed_sc="$2"
output="$3"
threads_ann="$4"
folder="${output%%/*}"

#vcf_file_gz="${snakemake_input[vcf_file]}.gz"
vcf_file_gz="$vcf_file.gz"
echo $vcf_file_gz

#bed_sc_gz="${snakemake_input[bed_sc]}.gz"
bed_sc_gz="$bed_sc.gz"
echo $bed_sc_gz

annotated1="$folder/annotated1.vcf.gz"
annotated2="$folder/annotated2.vcf.gz"
annotated3="$folder/annotated3.vcf.gz"

## first compress 
cat "$bed_sc" | bgzip > ${bed_sc_gz}
cat "$vcf_file" | bgzip > ${vcf_file_gz}


bcftools tabix -f "${bed_sc_gz}"
bcftools tabix -f "${vcf_file_gz}"
# -f force in case this is past the first run of the pipeline and these files already exist 

## and then annotate

bcftools annotate \
  -O z \
  -a "${bed_sc_gz}" \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  "${vcf_file_gz}" \
  -o "${annotated1}"


echo PHASE 4 amplify number of ecotypes to simulate starting seedmix
##### phase 4 #######################################
###### amplify number of ecotypes #########################################################
echo awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "${annotated1}" ## now only 225
bcftools tabix -f "${annotated1}"



bcftools merge --threads "${threads_ann}" --force-samples -O z "${annotated1}" "${annotated1}" "${annotated1}" "${annotated1}" -o "${annotated2}"
bcftools tabix -f "${annotated2}"
echo awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "${annotated2}"  ## 900


echo first amplification done 

bcftools merge --threads "${threads_ann}" --force-samples -O z "${annotated2}" "${annotated2}" "${annotated2}" "${annotated2}" -o "${annotated3}"
bcftools tabix -f "${annotated3}"
echo awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "${annotated3}"
  ## 3600

echo second amplification done 

bcftools merge --threads "${threads_ann}" --force-samples "${annotated3}" "${annotated3}" "${annotated3}" -o "${output}"
echo awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "${snakemake_output}"
   ## 10800


echo all amplifications done 
