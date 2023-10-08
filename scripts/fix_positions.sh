#!/bin/bash
output_vcf="${snakemake_input[output_vcf]}"
output_vcf_fixpos="${snakemake_output[output_vcf_fixpos]}"

awk -F'\t' 'BEGIN{OFS="\t"}
    $1 ~ /^#/ {print}
    $1 !~ /^#/ {
        pos=$2;
        if (pos > 30427671 && pos <= 50125959) {
            pos -= 30427671;
            $1="2"
        }
        else if (pos > 50125959 && pos <= 73585789) {
            pos -= 50125960;
            $1="3"
        }
        else if (pos > 73585789 && pos <= 92170845) {
            pos -= 73585790;
            $1="4"
        }
        else if (pos > 92170845 && pos <= 119146348) {
            pos -= 92170846;
            $1="5"
        }
        $2=pos;
        print
    }' $output_vcf > $output_vcf_fixpos

