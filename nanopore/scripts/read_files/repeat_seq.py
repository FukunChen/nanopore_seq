from Bio import SeqIO
from Bio.Seq import Seq

input_file = rf"D:\nanopore\data\fasta\idealDNA_repIU2N.fasta"
output_file = rf"D:\nanopore\data\fasta\idealDNA_repIU2N_repeat_read0.fasta"

prefix_seq = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"
count = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        header = "read{count}"
        count + 1
        original_seq = str(record.seq)
        repeated_seq = original_seq * 10
        new_seq = Seq(prefix_seq) + repeated_seq
        record.seq = new_seq
        SeqIO.write(record, outfile, "fasta")