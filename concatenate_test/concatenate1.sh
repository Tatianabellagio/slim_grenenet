
bcftools merge --threads 5 --force-samples -O z qtl2_contrib.vcf.gz qtl2_contrib.vcf.gz qtl2_contrib.vcf.gz qtl2_contrib.vcf.gz -o qtl2_contrib1.vcf.gz
bcftools tabix -f qtl2_contrib1.vcf.gz

bcftools merge --threads 5 --force-samples -O z qtl2_contrib1.vcf.gz qtl2_contrib1.vcf.gz qtl2_contrib1.vcf.gz qtl2_contrib1.vcf.gz -o qtl2_contrib_final1.vcf.gz
bcftools tabix -f qtl2_contrib_final1.vcf.gz
  
## 3600

bcftools merge --threads 5 --force-samples -O z qtl2_contrib_final1.vcf.gz qtl2_contrib_final1.vcf.gz qtl2_contrib_final1.vcf.gz qtl2_contrib_final1.vcf.gz -o qtl2_contrib_final10.vcf.gz
bcftools tabix -f qtl2_contrib_final10.vcf.gz
#14400