import pysam
import pod5
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import argparse
from Bio import SeqIO

read_seq_dict = {}
with open("read_sequences_from_fastq.fasta") as f_fasta:
    for record in SeqIO.parse(f_fasta, "fasta"):
        read_seq_dict[record.id] = str(record.seq)

def main():
    parser = argparse.ArgumentParser(description="Extract U/I reads from BAM and save as NPZ.")
    parser.add_argument("--bam_path", required=True, help="Path to the BAM file.")
    parser.add_argument("--bed_path", required=True, help="Path to the BED file containing U/I annotations.")
    parser.add_argument("--pod5_path", required=True, help="Path to the POD5 file.")
    args = parser.parse_args()
    
    bam_path = args.bam_path
    bed_path = args.bed_path
    pod5_path = args.pod5_path
    output_u_npz = "calls_u.npz"
    output_i_npz = "calls_i.npz"

    bam_fh = pysam.AlignmentFile(bam_path, "rb")
    read_damage_type = defaultdict(set)
    bed_read_ids = set()
    bed_reads = 0

    print("Read DNA damage message from bed file...")
    read_damage_positions = defaultdict(list)
    with open(bed_path) as f_bed:
        for line in f_bed:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue

            ctg, start, end, damage_type = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if damage_type not in {"U", "I"}:
                continue
            read_damage_positions[ctg].append((start, end, damage_type))
            bed_read_ids.add(ctg) 
            bed_reads += 1

            # for ctg in bam_fh:
            #     #if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            #     read_damage_type[ctg].add(damage_type)

    bam_fh.close()

    print(f"Find {bed_reads} reads have U/I damage")

    bam_fh = pysam.AlignmentFile(bam_path, "rb")
    pod5_reader = pod5.Reader(pod5_path)

    data_u = {"signal": [], "seq": [], "label": [], "move": []}
    data_i = {"signal": [], "seq": [], "label": [], "move": []}
    processed = set()

    print("Extract reads from BAM file...")
    for read in bam_fh.fetch(until_eof=True):
        read_id = read.query_name

        if read_id not in bed_read_ids:
            continue

        try:
            pod5_read = next(pod5_reader.reads(selection=[read_id]))
        except RuntimeError:
            print("pod5 error")
            continue

        try:
            move_table = read.get_tag("mv")
        except KeyError:
            print("No move table")
            continue

        if read.query_sequence is None:
            print("find seq in fastq")
            read_id = read.query_name
            sequence = read_seq_dict.get(read_id)
            if sequence is None:
                continue
            else:
                read.query_sequence = sequence 
        else:
            sequence = read.query_sequence
        #sequence = read.query_sequence
        signal = pod5_read.signal

        label_u = [0] * len(sequence)
        label_i = [0] * len(sequence)
        for start, end, dtype in read_damage_positions.get(read_id, []):
            if dtype == "U":
                for i in range(start, min(max(start + 1, end), len(sequence))):
                    label_u[i] = 1
            elif dtype == "I":
                for i in range(start, min(max(start + 1, end), len(sequence))):
                    label_i[i] = 1

        if any(label_u):
            data_u["seq"].append(list(sequence))
            data_u["signal"].append(signal)
            data_u["label"].append(label_u)
            data_u["move"].append(move_table)

        if any(label_i):
            data_i["seq"].append(list(sequence))
            data_i["signal"].append(signal)
            data_i["label"].append(label_i)
            data_i["move"].append(move_table)

        processed.add(read_id)

    print("Saving npz files...")

    np.savez(output_u_npz,
            signal=np.array(data_u["signal"], dtype=object),
            seq=np.array(data_u["seq"], dtype=object),
            label=np.array(data_u["label"], dtype=object),
            move=np.array(data_u["move"], dtype=object))

    np.savez(output_i_npz,
            signal=np.array(data_i["signal"], dtype=object),
            seq=np.array(data_i["seq"], dtype=object),
            label=np.array(data_i["label"], dtype=object),
            move=np.array(data_i["move"], dtype=object))

    print(f" calls_u.npz: {len(data_u['seq'])} U reads")
    print(f" calls_i.npz: {len(data_i['seq'])} I reads")
    print("All finished!")

if __name__ == "__main__":
    main()