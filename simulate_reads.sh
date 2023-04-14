#!/bin/bash
# Basic for loop
names='subp1.fasta subp3.fasta subp15.fasta subp11.fasta subp7.fasta subp5.fasta subp13.fasta subp2.fasta subp14.fasta subp0.fasta subp4.fasta subp12.fasta subp10.fasta subp6.fasta subp8.fasta subp9.fasta'
for name in $names
do
art_illumina -ss HS25 \
             -i vcf_slim/fasta_slim/$name \
             -l 150 \
             -f 7 \
             -p \
             -m 800 \
             -s 100 \
             -o vcf_slim/fasta_slim/reads_slim/$name

done
echo Reads generated