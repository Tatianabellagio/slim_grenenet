vcf_file="${snakemake_input[vcf_file]}"
bed_sc="${snakemake_input[bed_sc]}"
threads_ann="${snakemake[threads]}"
path="${snakemake_output}"
folder=$(dirname "$path")

echo "${folder}" > "${snakemake_log[0]}"

vcf_file_gz="${vcf_file}.gz"
#bed_sc_gz="$bed_sc.gz"

annotated1="${folder}/annotated1.vcf.gz"

## first compress 
cat "$bed_sc" | bgzip > ${bed_sc_gz} 2>> "${snakemake_log[0]}"
#cat "$vcf_file" | bgzip > ${vcf_file_gz} 2>> "${snakemake_log[0]}"


bcftools tabix -f "${bed_sc_gz}" 2>> "${snakemake_log[0]}"
#bcftools tabix -f "${vcf_file_gz}" 2>> "${snakemake_log[0]}"
# -f force in case this is past the first run of the pipeline and these files already exist 

## and then annotate

#echo done1 > "${snakemake_log[0]}"

bcftools annotate \
  -O z \
  -a "${bed_sc_gz}" \
  -c CHROM,FROM,TO,S \
  -h <(echo '##INFO=<ID=S,Number=.,Type=Float,Description="Selection Coefficient">') \
  "${vcf_file_gz}" \
  -o "${annotated1}" 2>> "${snakemake_log[0]}"


#echo PHASE 4 amplify number of ecotypes to simulate starting seedmix 2>> "${snakemake_log[0]}"
##### phase 4 #######################################
###### amplify number of ecotypes #########################################################
#echo awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "${annotated1}" ## now only 225
bcftools tabix -f "${annotated1}"

echo done2 >> "${snakemake_log[0]}"


output_file="${snakemake_output}"
file="${annotated1}"

echo "${output_file}" >> "${snakemake_log[0]}"
echo "${file}" >> "${snakemake_log[0]}"

command="bcftools merge --threads "${threads_ann}" --force-samples"

for i in {1..9}; do
    command+=" $file"
done

command+=" -o $output_file"

echo "$command" >> "${snakemake_log[0]}"
eval "$command" 2>> "${snakemake_log[0]}"

awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "${output_file}" >> "${snakemake_log[0]}"

