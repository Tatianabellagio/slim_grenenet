import os 
### generation of reference fasta from tair fasta with all chr

## params for gen of fasta from vcf 
or_fasta_file = 'chr5.fasta'
slim_fasta_file = 'slim_' + or_fasta_file

# Execute this block of code only if the file does not exist aka if we dont have the fasta file in slim format for the chr
if not os.path.exists(slim_fasta_file):
    ## save the first line (chr name) and all the other lines (seq) separatedly 
    with open(or_fasta_file, 'r') as file:
        chro = ""
        seq = ""
        for i, line in enumerate(file):
            if i == 0:
                chro += line.strip()
            elif i != 0:
                seq += line.strip()
    ## replace unknown variant with acgt so slim can read it 
    replacement_dict = {"M": "A", "R": "A", "W": "A", "S": "C", "Y": "C", "K": "G", "V": "A", "H": "A", "D": "A", "B": "C", "N": "A"}

    for key in replacement_dict:
        seq = seq.replace(key, replacement_dict[key])

    with open(slim_fasta_file, "w") as f:
        f.write(f'{chro}\n')
        f.write(seq + '\n')
        # Write the string to the file