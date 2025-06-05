import gzip
import os
from Bio import SeqIO

input_folder = rf"D:\nanopore\data\basecalling\pass"      
output_folder = rf"D:\nanopore\data\basecalling\bascalling_txt"   

os.makedirs(output_folder, exist_ok=True)

count = 0
for filename in os.listdir(input_folder):
    if filename.endswith(".fastq.gz"):
        count + 1
        input_path = os.path.join(input_folder, filename)
        output_filename = filename.replace(".fastq.gz", ".txt")
        output_path = os.path.join(output_folder, output_filename)
        
        print(f"Processing: {filename}")

        with gzip.open(input_path, "rt") as handle, open(output_path, "w") as outfile:
            for record in SeqIO.parse(handle, "fastq"):
                header = record.description
                sequence = str(record.seq)
                quality_ascii = record.letter_annotations["phred_quality"]

                outfile.write(f"@{header}\n")
                outfile.write(f"{sequence}\n")
                outfile.write(f"+\n")
                outfile.write(f"{quality_ascii}\n")

print(f"total file = {count}")