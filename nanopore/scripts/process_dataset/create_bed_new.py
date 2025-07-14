import os
from Bio import SeqIO
import pysam
import json
import argparse
from collections import defaultdict


def create_bed(ref_fa_path, bam_path, output_dir, filename):
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(ref_fa_path, "fasta")}

    coverage_dict = defaultdict(list)
    bed_lines = []
    found_reads = 0

    for aln in bam.fetch(until_eof=True):       
        if aln.is_unmapped or aln.query_sequence is None:
            continue

        if aln.is_secondary or aln.is_supplementary:
            continue
        
        #breakpoint()
        read_id = aln.query_name
        ref_name = aln.reference_name
        ref_seq = ref_dict.get(ref_name)
        if ref_seq is None:
            continue

        ref_idx = aln.reference_start
        #ref_idx = 0
        #print({ref_idx})
        ref_idx_list = []  # Store the aligned sequence index corresponding to the genome coordinate

        for base in ref_seq:
                ref_idx_list.append(ref_idx) ###from the starting position(the read) of reference store the following bases
                ref_idx += 1
        
        #breakpoint()
        for i, base in enumerate(ref_seq):
            if base not in ['U', 'I']:
                continue
            
            #breakpoint()
            pos_start = max(0, i)
            pos_end = min(len(ref_seq)-1, i + 1)
            context_ref  = base
            #print(f"DNA damage is: {context_ref} and it's basecalled as : and its position is: {pos_end} ")
        
            if pos_start is not None and pos_end is not None:
                bed_lines.append(f"{ref_name}\t{pos_start}\t{pos_end+1}")
            
            # if x_pos_start is not None and x_pos_end is not None:
            #     bed_lines.append(f"{ref_name}\t{x_pos_start}\t{x_pos_end+1}")

            # Save details to JSON
            entry = {
                "read_id": read_id,
                f"{base}_context": context_ref,
                f"{base}_pos": [pos_start, pos_end],
            }

            coverage_dict[f"{ref_name}:{ref_idx_list[i]+1}"].append(entry)
            
        found_reads += 1

    bam.close()

    # Output JSON and BED files
    output_json = os.path.join(output_dir, filename)
    output_bed = os.path.join(output_dir, filename)

    with open(output_json, "w") as f:
        json.dump(coverage_dict, f, indent=2)

    with open(output_bed, "w") as f:
        for line in bed_lines:
            f.write(line + "\n")

    print(f"\nDetected {found_reads} mapped reads with U/I.")
    print(f"Output saved: {output_json}, {output_bed}")

def main(args):
    output_dir = args.output_dir
    ref_fa_path = args.ref_fa_path
    bam_path = args.bam_path
    filename = args.filename

    create_bed(ref_fa_path, bam_path, output_dir, filename)
    print("\nAll finished!!!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BAM and reference to generate BED output")

    parser.add_argument("--output_dir", type=str, required=True, help="Directory to store output BED")
    parser.add_argument("--ref_fa_path", type=str, required=True, help="Path to reference FASTA file")
    parser.add_argument("--bam_path", type=str, required=True, help="Path to input BAM file")
    parser.add_argument("--filename", type=str, required=True, help="Output BED file name")

    args = parser.parse_args()
    main(args)
