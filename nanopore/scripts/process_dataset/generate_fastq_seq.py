import os
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_reads_from_fastq_gz_dir(fastq_dir):
    read_dict = {}
    for filename in os.listdir(fastq_dir):
        if filename.endswith(".fastq.gz") or filename.endswith(".fq.gz"):
            filepath = os.path.join(fastq_dir, filename)
            #print(f"Processing: {filepath}")
            with gzip.open(filepath, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    read_id = record.id.split()[0]
                    if read_id not in read_dict: 
                        read_dict[read_id] = str(record.seq)
    return read_dict

fastq_folder = r"D:\nanopore\data\pod5\nanopore\basecalling\pass"
read_seq_dict = extract_reads_from_fastq_gz_dir(fastq_folder)

output_fasta = "read_sequences_from_fastq.fasta"
records = [SeqRecord(Seq(seq), id=read_id, description="") for read_id, seq in read_seq_dict.items()]

with open(output_fasta, "w") as f:
    SeqIO.write(records, f, "fasta")

print(f"Saved {len(records)} unique reads to {output_fasta}")