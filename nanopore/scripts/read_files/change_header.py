from Bio import SeqIO
from Bio.Seq import Seq

input_file = rf"D:\nanopore\data\fasta\idealDNA_repIU2N.fasta"
output_file = rf"D:\nanopore\data\fasta\idealDNA_repIU2N_repeat_ch.fasta"

prefix_seq = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for i, record in enumerate(SeqIO.parse(infile, "fasta")):
        record.id = f"read{i}"
        record.name = f"read{i}"
        record.description = f"read{i}"
        
        original_seq = str(record.seq)
        repeated_seq = original_seq * 10
        new_seq = Seq(prefix_seq) + repeated_seq
        record.seq = new_seq
        SeqIO.write(record, outfile, "fasta")