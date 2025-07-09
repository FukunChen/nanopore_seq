from Bio import SeqIO
import os
import argparse

def main(arg):
    input_folder = args.input_folder
    output_file = args.output_file

    all_records = []

    for filename in os.listdir(input_folder):
        if filename.endswith(".fasta"):
            filepath = os.path.join(input_folder, filename)
            print(f"Loading: {filename}")
            records = list(SeqIO.parse(filepath, "fasta"))
            all_records.extend(records)


    SeqIO.write(all_records, output_file, "fasta")
    print(f"Num: {len(all_records)} -> {output_file}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine all FASTA files in a folder into a single FASTA file.")
    parser.add_argument("--input_folder", type=str, required=True, help="Path to folder containing .fasta files")
    parser.add_argument("--output_file", type=str, required=True, help="Path to save the combined .fasta file")
    
    args = parser.parse_args()
    main(args)