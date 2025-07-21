import os
from Bio import SeqIO
import pysam
import json
from collections import defaultdict
import argparse

read_seq_dict = {}
with open("read_sequences_from_fastq.fasta") as f_fasta:
    for record in SeqIO.parse(f_fasta, "fasta"):
        read_seq_dict[record.id] = str(record.seq)
        #print("Get read from fastq")

def reconstruct_alignment(read, ref_seq, output_path="alignment.json"):
    aligned_ref = []
    aligned_read = []

    ref_pos = read.reference_start ###reference.start is the position where read starts on reference
    #print(f"reference start is {ref_pos}")
    read_pos = 0
    #breakpoint()
    
    aligned_ref = []
    aligned_read = []

    for (cigar_type, length) in read.cigartuples:
        if cigar_type == 0:  # match / mismatch (M)
            # Both sides advance 1
            aligned_ref.extend(ref_seq[ref_pos : ref_pos + length])
            aligned_read.extend(read.query_sequence[read_pos : read_pos + length])
            ref_pos += length
            read_pos += length

        elif cigar_type == 1:  # insertion to reference (I)
            # If ref is empty, use '-'
            aligned_ref.extend(['-'] * length)
            aligned_read.extend(read.query_sequence[read_pos : read_pos + length])
            read_pos += length

        elif cigar_type == 2:  # deletion from reference (D)
            # read is blank, use '-'
            aligned_ref.extend(ref_seq[ref_pos : ref_pos + length])
            aligned_read.extend(['-'] * length)
            ref_pos += length
        elif cigar_type == 4:  # soft clip (S)
            aligned_ref.extend(['-'] * length)
            aligned_read.extend(read.query_sequence[read_pos : read_pos + length])
            read_pos += length

        elif cigar_type == 5:  # hard clip (H)
            # Do not consume ref, do not consume read, skip
            continue

        else:
            print(f"Unhandled CIGAR type: {cigar_type}")

    return aligned_ref, aligned_read



def create_bed(ref_fa_path, bam_path, output_dir, filename):
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(ref_fa_path, "fasta")}

    coverage_dict = defaultdict(list)
    bed_lines = []
    found_reads = 0
    no_reference = 0

    for aln in bam.fetch(until_eof=True):
        aligned_ref = []
        aligned_read = []

        if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
        #if aln.is_unmapped:
              continue

        if aln.query_sequence is None:
            #print("find sequence in fastq")
            read_id = aln.query_name
            sequence = read_seq_dict.get(read_id)
            if sequence is None:
                continue
            else:
                aln.query_sequence = sequence 
            #continue
        else:
            sequence = aln.query_sequence
        
        read_id = aln.query_name
        ref_name = aln.reference_name
        if ref_name is None:
            no_reference += 1
        ref_seq = ref_dict.get(ref_name)
        if ref_seq is None:
            continue


        aligned_ref, aligned_read = reconstruct_alignment(aln, ref_seq)

        ref_idx = aln.reference_start
        ref_idx_list = []  # Store the aligned sequence index corresponding to the genome coordinate

        for base in aligned_ref:
                ref_idx_list.append(ref_idx) #
                ref_idx += 1
                
        for i, base in enumerate(aligned_ref):
            if base not in ['U', 'I']:
                continue
            
            #print("find I/U processing...")
            i_u = sum(1 for b in aligned_read[0:i+1] if b != '-') ##Convert the position information in reference to the position information in read
            if i_u >= len(sequence):
                continue
            pos_start = i_u
            pos_end = i_u + 1
            if pos_end >= len(sequence):
                pos_end = i_u
            context_ref  = base
            context_read = aligned_read[i]
            #print(f"DNA damage is: {context_ref} and it's basecalled as : {context_read} and its position is: {i_u} ")
            
            # Find the position after the 15th reference base
            # ref_count = 0
            # j = i
            # while j < len(aligned_ref) and ref_count < 15:
            #     if aligned_ref[j] != '-':
            #         ref_count += 1
            #     j += 1 

            # if ref_count == 15 and j < len(aligned_ref):
            #     x = sum(1 for b in aligned_read[0:j+1] if b != '-')
            #     x_pos_start = x
            #     x_pos_end = x + 1
            #     x_context = aligned_ref[j]
            #     #print(f"control base for {context_ref} is : {x_context} and its position is : {x}")
            #     if x_pos_end or x_pos_start  >= len(aln.query_sequence):
            #         x_pos_start = None
            #         x_pos_end = None
                    
                    
            # aligned_ref_str = ''.join(aligned_ref)
            # aligned_read_str = ''.join(aligned_read)

            # print(f"Aligned ref:  {aligned_ref_str}")
            # print(f"Aligned read: {aligned_read_str}")

            if pos_start is not None and pos_end is not None:
                bed_lines.append(f"{ref_name}\t{pos_start}\t{pos_end}\t{base}")
            
            # if x_pos_start is not None and x_pos_end is not None:
            #     bed_lines.append(f"{ref_name}\t{x_pos_start}\t{x_pos_end}\t{base}")

            entry = {
                "read_id": read_id,
                f"{base}_context": context_ref,
                f"{base}_context__read": context_read,
                f"{base}_pos": [pos_start, pos_end],
                # f"{base}_X_context": x_context,
                # f"{base}_X_pos": [x_pos_start, x_pos_end]
            }
            # if context_ref == "I":
            #     print(f"Saving to coverage_dict[{ref_name}:{ref_idx_list[i]+1}]:")
            #     for k, v in entry.items():
            #         print(f"  {k}: {v}")
                
                #breakpoint()
            coverage_dict[f"{ref_name}:{ref_idx_list[i]+1}"].append(entry)
            found_reads += 1

    bam.close()

    output_json = os.path.join(output_dir, filename + ".json")
    output_bed = os.path.join(output_dir, filename + ".bed")

    with open(output_json, "w") as f:
        json.dump(coverage_dict, f, indent=2)

    with open(output_bed, "w") as f:
        for line in bed_lines:
            f.write(line + "\n")

    print(f"\n There are {no_reference} reads don't have reference in mapped reads")
    print(f"\nDetected {found_reads} mapped reads with U/I.")
    print(f"Output saved: {output_json}, {output_bed}")

def main():
    parser = argparse.ArgumentParser(description="Extract U/I positions and write to BED and JSON.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the output BED and JSON files.")
    parser.add_argument("--ref_fa_path", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--bam_path", required=True, help="Path to the input BAM file.")
    parser.add_argument("--filename", required=True, help="Base name (no extension) for output files.")

    args = parser.parse_args()

    create_bed(args.ref_fa_path, args.bam_path, args.output_dir, args.filename)
    print("\nAll finished!!!")

if __name__ == "__main__":
    main()