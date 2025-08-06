import os
from Bio import SeqIO
import pysam
import json
import argparse
from collections import defaultdict
from parasail_function import (GAP_OPEN_PENALTY,
                               GAP_EXTEND_PENALTY,
                               LOCAL_ALIGN_FUNCTION,
                               MATRIX
                               )
import pysam
import argparse
from Bio import SeqIO
from collections import defaultdict
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

CHUNK_SIZE = 5000
################################################################
# Use parasail to find the corresponding u position in the read sequence and make an index
###############################################################
read_seq_dict = {}
with open("read_sequences_from_fastq.fasta") as f_fasta:
    for record in SeqIO.parse(f_fasta, "fasta"):
        read_seq_dict[record.id] = str(record.seq)
        #print("Get read from fastq")

# def create_bed(ref_fa_path, bam_path, output_dir, filename):
#     bam = pysam.AlignmentFile(bam_path, "rb")
#     ref_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(ref_fa_path, "fasta")}

#     coverage_dict = defaultdict(list)
#     bed_lines = []
#     found_reads = 0

#     for aln in bam.fetch(until_eof=True):       
#         if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
#             continue

#         if aln.query_sequence is None:
#             #print("find sequence in fastq")
#             read_id = aln.query_name
#             sequence = read_seq_dict.get(read_id)
#             if sequence is None:
#                 no_reference += 1
#                 continue
#             else:
#                 aln.query_sequence = sequence 
#             #continue
#         else:
#             sequence = aln.query_sequence
        
#         #breakpoint()
#         read_id = aln.query_name
#         ref_name = aln.reference_name
#         ref_seq = ref_dict.get(ref_name)
#         if ref_seq is None:
#             continue
        
#         result = LOCAL_ALIGN_FUNCTION(sequence, ref_seq, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, MATRIX)
#         ref = result.ref
#         query = result.query

#         print(f"the length of ref is {len(ref)}, the length of query is {len(query)}")
#         #breakpoint()
#         for i, base in enumerate(ref):
#             if base not in ['U', 'I']:
#                 continue
            
#             #breakpoint()
#             pos_start = sum(1 for q in query[:i] if q != '-')
#             if pos_start > 200:
#                 continue
#             pos_end = pos_start + 1
            
#             #print(f"DNA damage is: {context_ref} and it's basecalled as : and its position is: {pos_end} ")
        
#             if pos_start is not None and pos_end is not None:
#                 bed_lines.append(f"{read_id}\t{pos_start}\t{pos_end+1}")
            
#         found_reads += 1

#     bam.close()

#     # Output JSON and BED files
#     output_bed = os.path.join(output_dir, filename)

#     with open(output_bed, "w") as f:
#         for line in bed_lines:
#             f.write(line + "\n")

#     print(f"\nDetected {found_reads} mapped reads with U/I.")
#     print(f"Output saved: {output_bed}")

# def main(args):
#     output_dir = args.output_dir
#     ref_fa_path = args.ref_fa_path
#     bam_path = args.bam_path
#     filename = args.filename

#     create_bed(ref_fa_path, bam_path, output_dir, filename)
#     print("\nAll finished!!!")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Process BAM and reference to generate BED output")

#     parser.add_argument("--output_dir", type=str, required=True, help="Directory to store output BED")
#     parser.add_argument("--ref_fa_path", type=str, required=True, help="Path to reference FASTA file")
#     parser.add_argument("--bam_path", type=str, required=True, help="Path to input BAM file")
#     parser.add_argument("--filename", type=str, required=True, help="Output BED file name")

#     args = parser.parse_args()
#     main(args)
    
def process_chunk(chunk_reads, ref_dict):
    """处理 BAM chunk，返回 BED 结果"""
    bed_lines = []
    found_reads = 0

    for aln in chunk_reads:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue

        sequence = aln.query_sequence
        if sequence is None:
            read_id = aln.query_name
            sequence = read_seq_dict.get(read_id)
            if sequence is None:
                continue

        read_id = aln.query_name
        ref_name = aln.reference_name
        ref_seq = ref_dict.get(ref_name)
        if ref_seq is None:
            continue

        # 本地比对
        result = LOCAL_ALIGN_FUNCTION(sequence, ref_seq, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, MATRIX)
        ref = result.ref
        query = result.query

        for i, base in enumerate(ref):
            if base not in ['U', 'I']:
                continue

            pos_start = sum(1 for q in query[:i] if q != '-')
            if pos_start > 200:
                continue
            pos_end = pos_start + 1
            bed_lines.append(f"{read_id}\t{pos_start}\t{pos_end+1}")

        found_reads += 1

    return bed_lines, found_reads


def create_bed(ref_fa_path, bam_path, output_dir, filename):
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_dict = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(ref_fa_path, "fasta")}

    print("Loading BAM alignments into memory...")
    all_alignments = list(bam.fetch(until_eof=True))
    bam.close()

    total_reads = len(all_alignments)
    print(f"Total reads: {total_reads}")

    chunks = [all_alignments[i:i + CHUNK_SIZE] for i in range(0, total_reads, CHUNK_SIZE)]

    bed_lines = []
    found_reads_total = 0

    with Pool(processes=cpu_count()) as pool:
        for bed_result, found_reads in tqdm(pool.imap_unordered(
            lambda chunk: process_chunk(chunk, ref_dict), chunks),
            total=len(chunks),
            desc="Processing reads"
        ):
            bed_lines.extend(bed_result)
            found_reads_total += found_reads

    output_bed = os.path.join(output_dir, filename)
    with open(output_bed, "w") as f:
        for line in bed_lines:
            f.write(line + "\n")

    print(f"\nDetected {found_reads_total} mapped reads with U/I.")
    print(f"Output saved: {output_bed}")


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
