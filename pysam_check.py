import pysam

# Open FASTA and VCF files
fasta_file = pysam.FastaFile("genome.fasta")
vcf_file = pysam.VariantFile("variants.vcf")

# Count number of positions in FASTA and VCF files
num_positions_fasta = sum([len(fasta_file[contig]) for contig in fasta_file.references])
num_positions_vcf = sum([1 for record in vcf_file])

# Close files
fasta_file.close()
vcf_file.close()

# Compare number of positions
if num_positions_fasta == num_positions_vcf:
    print("The FASTA and VCF files have the same number of positions.")
else:
    print("The FASTA and VCF files do not have the same number of positions.")
