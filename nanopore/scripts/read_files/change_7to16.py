from Bio import SeqIO
from Bio.Seq import Seq

input_file = rf"D:\nanopore\data\fasta\idealDNA.fasta"
output_file = rf"D:\nanopore\data\fasta\idealDNA_rep7to16.fasta"
records = []

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for i, record in enumerate(SeqIO.parse(infile, "fasta")):
        seq = str(record.seq)
        if record.id.endswith("fwd"):
            new_seq = seq[:6] + 'N' * 10 + seq[16:]
            record.seq = Seq(new_seq)
        records.append(record)
    SeqIO.write(records, outfile, "fasta")