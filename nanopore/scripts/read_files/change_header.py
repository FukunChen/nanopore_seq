from Bio import SeqIO
from Bio.Seq import Seq

input_file = rf"D:\nanopore\data\fasta\idealDNA.fasta"
output_file = rf"D:\nanopore\data\fasta\idealDNA_new.fasta"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('>'):
            # 去除标题行中的空格
            cleaned_line = line.replace(' ', '')
            outfile.write(cleaned_line)
        else:
            outfile.write(line)