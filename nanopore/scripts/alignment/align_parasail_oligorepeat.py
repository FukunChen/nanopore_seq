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
import math
import argparse
import matplotlib.pyplot as plt
from parasail_function import (GAP_OPEN_PENALTY,
                               GAP_EXTEND_PENALTY,
                               LOCAL_ALIGN_FUNCTION,
                               MATRIX
                               )

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

def find_match_ratio(arr):
    matches_nu = 0
    matches = 0
    aligned_len = 0
    
    match_symbols = np.where(np.isin(arr[1], ['|', ':', '.']))[0]
    if len(match_symbols) == 0:
        return 0
    arr = arr[:, np.min(match_symbols):np.max(match_symbols)+1]
    matches = np.sum(np.isin(arr[1], ['|']))
    matches_nu = np.sum(np.isin(arr[1], [':']))
    aligned_len = arr.shape[1] - matches_nu
    match_ratio = matches / aligned_len if aligned_len > 0 else 0
    return match_ratio

#find header/tail adapter in the read
def find_adapter(header_adapter, tail_adapter, query):
    header = 0
    tail = 0

    if header_adapter:
        result_header = LOCAL_ALIGN_FUNCTION(query, header_adapter[0], GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, MATRIX)
        tb = result_header.traceback

        q_aln = tb.query
        r_aln = tb.ref
        comp = tb.comp

        # Trimming arr
        arr = np.array([list(q_aln), list(comp), list(r_aln)])
        match_ratio = find_match_ratio(arr)
        
        match_indices = np.where(np.isin(arr[1], ['|', ':', '.']))[0]
        if len(match_indices) == 0:
            return None, query, ""
        aln_end = np.max(match_indices) + 1
        query_aligned_segment = arr[0, :aln_end]  
        non_gap_bases_count = np.sum(query_aligned_segment != '-')  

        if non_gap_bases_count <= 40 and match_ratio > 0.5:  # the lenth of header adapter is 36. is it shoould be 0.9 or 0.5 ? but the posiotion is here
            header = 1

    if tail_adapter:
        result_tail = LOCAL_ALIGN_FUNCTION(query, tail_adapter[0], GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, MATRIX)
        tb = result_tail.traceback
        q_aln = tb.query
        r_aln = tb.ref
        comp = tb.comp

        arr = np.array([list(q_aln), list(comp), list(r_aln)])
        match_ratio = find_match_ratio(arr)

        match_indices = np.where(np.isin(arr[1], ['|', ':', '.']))[0]
        if len(match_indices) == 0:
            return None, query, ""
        aln_start = np.min(match_indices)
        query_aligned_segment = arr[0, :aln_start]
        non_gap_bases_count = np.sum(query_aligned_segment != '-')  
        distance_to_tail = len(query) - non_gap_bases_count

        if distance_to_tail <= 30 and match_ratio > 0.5: # the lenth of tail is 27 
            tail = 1

    return header, tail


#process read name to avaliable format
def sanitize_filename(name):
    return re.sub(r'[\\/*?:"<>|= @]', '_', name)

#alignment for one read
def find_oligo(read, oligos):
    #read = str(read.seq)
    best_score = 0
    best_oligo = None
    best_result = None
    #pre-process the oligos
    # if oligo.id endwith '_fwd', change [7,16] to N
    # if oligo.id endwith '_rev', change [19,28] to N

    for oligo in oligos:  
        original_oligo = oligo
        
        #processed the oligo
        oligo_seq = list(oligo.seq)  
        if oligo.id.endswith('_fwd'):
            for i in range(7, 17): 
                if i < len(oligo_seq):
                    oligo_seq[i] = 'N'
        elif oligo.id.endswith('_rev'):
            for i in range(19, 29):
                if i < len(oligo_seq):
                    oligo_seq[i] = 'N'
        oligo_seq = ''.join(oligo_seq)   
        
        in_oligo = str(oligo_seq)
        #print(f"oligo seq: {in_oligo}")
        result = LOCAL_ALIGN_FUNCTION(read, in_oligo, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, MATRIX)
        ref = result.ref
        query = result.query

        
        current_score = result.score
        
        #testing traceback
        tb = result.traceback
        if tb is None:
            print("traceback is None")
        else:
            try:
                _ = tb.ref
                _ = tb.query
            except ValueError:
                print("traceback exists, but internal .ref or .query is NULL (C-level)")
                break
        
        q_aln = tb.query
        r_aln = tb.ref
        comp = tb.comp

        arr = np.array([list(q_aln), list(comp), list(r_aln)])
        #check count results
        mr = find_match_ratio(arr)
        # matched_part = np.where(np.isin(list(comp),['|'])) [0]
        # current_score = matched_part.shape[0]
        #print(f"current score: {current_score}")
        if current_score > best_score and mr == 1:  ### only when the oligo entirely matched with the read, it can be identified as the Best Oligo
            best_score = current_score
            best_oligo = original_oligo
            best_result = result
    
    if best_result is None or best_score <= 0:
        return None
    

    return best_oligo


def main(args):
    ### set input
    fastq_folder = args.fastq_folder
    fasta_folder = args.fasta_file
    output_folder = args.output_folder
    
    oligos = []
    querys = []
    headers = []
    oligos = read_from_fasta(fasta_folder)
    
    #if there is header or tail in the read or not
    header_yn = 0
    tail_yn = 0
    header_adapter = ['TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTT']  ##end repire add T at 5'
    tail_adapter = ['GCAATACGTAACTGAACGAAGTACAGGA'] ## end repire add A at 3'
    ref_len = 0
    
    #check match rate for all of the files
    match_ratios = []
    
    #process only one read in one loop
    for filename in os.listdir(fastq_folder):
        #parameter for summary
        total_reads = 0
        matched_reads = 0
        reads_mr_eight = 0
        reads_mr_seven = 0
        if filename.endswith(".fastq.gz"):
            fastq_path = os.path.join(fastq_folder, filename)
        querys,headers = read_from_fastq(fastq_path)
        
        all_results = []
        fasta_records = []
        json_results = []

        for query, header in zip(querys, headers):
            #initial
            ref_len = 0
            total_reads = total_reads + 1
            #Filter out reads that are too long
            if len(query) > 500:
                continue
            header_yn, tail_yn = find_adapter(header_adapter, tail_adapter, query)
            best_oligo = find_oligo(query, oligos)
            if best_oligo is None:
                continue
            
            if header_yn == 1 or tail_yn == 1:
                ref = ''
                if header_yn == 1:
                    ref += header_adapter[0]
                    ref_len = len(ref)
                if tail_yn == 1:
                    ref_len = len(ref) + 27 - 10
                
                oligo_len = len(best_oligo.seq)
                needed_len = len(query) - ref_len + 10
                repeat_times = math.ceil(needed_len / oligo_len)   
                ref += str(best_oligo.seq) * repeat_times
                
                ref = ref[:-10] #cut the hang-out part
                if tail_yn == 1:
                    ref += tail_adapter[0]

            # without adapter
            else:
                repeat_times = math.ceil(len(query) / len(best_oligo.seq))
                ref = str(best_oligo.seq) * repeat_times
            # if len(ref) > len(query):
            #     print('ref is longer than query')
            
            result = LOCAL_ALIGN_FUNCTION(query, ref, GAP_OPEN_PENALTY, GAP_EXTEND_PENALTY, MATRIX)
            tb = result.traceback
            q_aln = tb.query
            r_aln = tb.ref
            comp = tb.comp

            arr = np.array([list(q_aln), list(comp), list(r_aln)])
            
            #check match rate
            match_r = find_match_ratio(arr)
            match_ratios.append(match_r)
            matched_reads = matched_reads + 1

            #print matched results
            # print("\nAlignment Result:")
            # print("Query : ", q_aln)
            # print("Comp  : ", comp)
            # print("Ref   : ", r_aln)
            
            #save reference: json and fasta
            json_results.append({
                        "read_id": header,
                        "query_aln": q_aln,
                        "comp": comp,
                        "ref_aln": r_aln
                    })
            corrected_ref = list(r_aln)
            for i in range(len(corrected_ref)):
                if corrected_ref[i] == 'N' and q_aln[i] != '-':
                    corrected_ref[i] = q_aln[i]
            corrected_ref_str = ''.join(corrected_ref)
            fasta_records = SeqRecord(
                Seq(corrected_ref_str),
                id=f"{sanitize_filename(header)}",
                description=""
            )
            all_results.append(fasta_records)

                    

        print(f'--------------------For {filename} ----------------------') 
        print(f'Total reads (per fastq) is: {total_reads}')
        print(f'Matched_reads with best oligos: {matched_reads}')

        # output results
        safe_name = sanitize_filename(filename.replace(".fastq.gz", ""))
        shared_fasta_path = os.path.join(output_folder, f"{safe_name}.fasta")
        shared_json_path = os.path.join(output_folder, f"{safe_name}.json")

        with open(shared_fasta_path, "w") as fasta_out:
            SeqIO.write(all_results, fasta_out, "fasta")

        with open(shared_json_path, "w") as jf:
             json.dump(json_results, jf, indent=4)
             
    # if match_ratios:
    #     sorted_ratios = np.sort(match_ratios)
    #     cdf = np.arange(1, len(sorted_ratios)+1) / len(sorted_ratios)

    #     plt.figure(figsize=(8, 6))
    #     plt.plot(sorted_ratios, cdf, marker='.', linestyle='-')
    #     plt.xlabel("Match Ratio")
    #     plt.ylabel("CDF")
    #     plt.title("CDF of Match Ratio")
    #     plt.grid(True)
    #     plt.tight_layout()
    #     plt.show()
    # else:
    #     print("No match ratios to plot.") 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run alignment using parasail and custom matrix")
    
    #parser.add_argument("--matrix", type=str, required=True, help="Path to substitution matrix file (.txt)")
    parser.add_argument("--fastq_folder", type=str, required=True, help="Path to folder containing FASTQ files")
    parser.add_argument("--fasta_file", type=str, required=True, help="Path to reference FASTA file")
    parser.add_argument("--output_folder", type=str, required=True, help="Folder to store output alignments")

    args = parser.parse_args()
    main(args)
