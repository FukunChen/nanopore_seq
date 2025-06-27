import parasail
import numpy as np
import csv
import json
from Bio import SeqIO
import os
import gzip
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

#input setting
sub_matrix = parasail.Matrix(r"D:\nanopore\scripts\alignment\didu_matrix_withIU.txt")
#sub_matrix = parasail.Matrix(r"D:\nanopore\scrips\alignment\didu_matrix.txt")
fastq_folder = rf"D:\nanopore\data\pod5\nanopore\basecalling\pass" 
fasta_folder = rf"D:\nanopore\data\fasta\idealDNA.fasta"
#fasta_folder = rf"D:\nanopore\data\fasta\idealDNA_rep7to16.fasta"
output_folder = rf"D:\nanopore\data\aligned_reference"


        
#check matrix
valid_chars = sub_matrix.matrix
print(f"valid_chars: {valid_chars}")

#read files
def read_from_fastq(fastq_path):
    reads = []
    headers = []
    with gzip.open(fastq_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            headers.append(record.id)
            reads.append(str(record.seq))
                    
    return reads,headers

def read_from_fasta(input_file):
    refs = []
    with open(input_file, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"): 
            refs.append(record)
    return refs


#process read name to avaliable format
def sanitize_filename(name):
    return re.sub(r'[\\/*?:"<>|= @]', '_', name)

#alignment for one read
def align_and_mask_with_output(read, oligos, read_id):
    max_iter = (len(read) // 46) + 2 # 
    min_remaining=4  #why 15???
    original_read = list(read)
    read = list(read)
    mask = [False] * len(read)
    results = []
    aligned_seq = ['-'] * len(read)
    iteration = 0
    selected_oligos = oligos
    header_matched = False
    tail_matched = False 
        
    for iteration in range(max_iter):
        current_read = ''.join([read[i] if not mask[i] else 'X' for i in range(len(read))])
        #print(f"current read is {current_read}")
        if current_read.count('X') > len(read) - min_remaining:
            break

        best_score = -float('inf')
        best_oligo = None
        best_result = None        
        
        #divide oligos
        fwd_oligos = [o for o in oligos if o.description.endswith("fwd")]
        rev_oligos = [o for o in oligos if o.description.endswith("rev")]
        
        #if header or tail is detected, remove them in the oligos
        selected_oligos = [
        o for o in oligos
        if not (header_matched and o.description == "header")
        and not (tail_matched and o.description == "tail")
        ]

        for oligo in selected_oligos:        
            in_oligo = str(oligo.seq)
            #print(f"oligo seq: {in_oligo}")
            result = parasail.sg_qx_trace_striped_16(current_read, in_oligo, 10, 5, sub_matrix)
            ref = result.ref
            query = result.query
            #current_score = compute_partial_score(ref,query) 
            current_score = result.score
            
            #testing traceback
            tb = result.traceback
            if tb is None:
                print("traceback is None")
            else:
                try:
                    print("ref:", tb.ref)
                    print("query:", tb.query)
                except ValueError:
                    print("traceback exists, but internal .ref or .query is NULL (C-level)")
                    break
                    
            #check count results
            
            if current_score > best_score:
                best_score = current_score
                best_oligo = oligo
                best_result = result
    
        if best_result is None:
            return None, []
        if best_score <= 0:
            break
        
        #check header
        if best_oligo is not None:
            if best_oligo.description == "header":
                header_matched = True
            elif best_oligo.description == "tail":
                tail_matched = True
        #select continue oligos
        if iteration == 0 and best_oligo is not None:
            if best_oligo.description.endswith("fwd"):
                selected_oligos = fwd_oligos
            elif best_oligo.description.endswith("rev"):
                selected_oligos = rev_oligos
                    
        # use traceback to compute matched positions
        if best_result is None:
            return None, []
        query_aln = best_result.traceback.query
        ref_aln = best_result.traceback.ref
        read_pos = 0
        mask_indices = []

        for q_char, r_char in zip(query_aln, ref_aln):
            if q_char != '-':
                if r_char != '-' and r_char != 'N':
                    mask_indices.append(read_pos)
                    aligned_seq[read_pos] = r_char
                elif r_char == 'N':
                    mask_indices.append(read_pos)
                    aligned_seq[read_pos] = q_char
                read_pos += 1  # only advance when not a gap in query

        # mark matched positions
        for i in mask_indices:
            if 0 <= i < len(read):
                mask[i] = True
                read[i] = "X"

        #save result for each iteration
        results.append({
            "read_id": read_id,
            "iteration": iteration + 1,
            "oligo id":str(best_oligo.id),
            "oligo": str(best_oligo.seq),
            "score": best_result.score,
            "query_aln": query_aln,
            "ref_aln": ref_aln
        })
        
        #in each iteration, update final_reference
        final_seq = ''.join(aligned_seq)
        
        #iteration count
        #iteration = iteration + 1
    
    #save reference result after all iterations
    final_record = SeqRecord(
        Seq(final_seq),
        id=f"{sanitize_filename(read_id)}",
        description=""
    )

    return final_record, results

def main():
    oligos = []
    querys = []
    headers = []
    oligos = read_from_fasta(fasta_folder)
    
    #process only one read in one loop
    for filename in os.listdir(fastq_folder):
        if filename.endswith(".fastq.gz"):
            fastq_path = os.path.join(fastq_folder, filename)
        querys,headers = read_from_fastq(fastq_path)
        
        all_results = []
        fasta_records = []

        for query, header in zip(querys, headers):
            if len(query) > 400:
                continue
            final_record, results = align_and_mask_with_output(query, oligos, header)
            if final_record is None:
                continue
            fasta_records.append(final_record)
            all_results.extend(results)

       
        safe_name = sanitize_filename(filename.replace(".fastq.gz", ""))
        shared_fasta_path = os.path.join(output_folder, f"{safe_name}.fasta")
        shared_json_path = os.path.join(output_folder, f"{safe_name}.json")

        with open(shared_fasta_path, "w") as fasta_out:
            SeqIO.write(fasta_records, fasta_out, "fasta")

        with open(shared_json_path, "w") as jf:
            json.dump(all_results, jf, indent=4)

if __name__ == '__main__':
    main()