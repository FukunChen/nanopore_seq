import numpy as np
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import gzip
import re
import math
import argparse
import matplotlib.pyplot as plt
import concurrent.futures
import time
from tqdm import tqdm  # 用于进度条
from parasail_function import (GAP_OPEN_PENALTY,
                               GAP_EXTEND_PENALTY,
                               LOCAL_ALIGN_FUNCTION,
                               MATRIX
                               )

def safe_str(x):
    return ''.join(x) if isinstance(x, list) else x

def parse_args():
    parser = argparse.ArgumentParser(description="Process FASTQ files with oligo matching.")
    parser.add_argument('--fastq_folder', required=True, help='Path to folder containing .fastq.gz files')
    parser.add_argument('--fasta_file', required=True, help='Path to FASTA file containing oligos')
    parser.add_argument('--output_folder', required=True, help='Path to output folder')
    return parser.parse_args()

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
    new_ref = []
    
    match_symbols = np.where(np.isin(arr[1], ['|', ':', '.']))[0]
    if len(match_symbols) == 0:
        return 0
    arr = arr[:, np.min(match_symbols):np.max(match_symbols)+1]
    matches = np.sum(np.isin(arr[1], ['|']))
    matches_nu = np.sum(np.isin(arr[1], [':']))
    aligned_len = arr.shape[1] - matches_nu
    match_ratio = matches / aligned_len if aligned_len > 0 else 0
    #print(f'check arr[2]: {arr[2]}')
    for i in range(len(arr[2])):
        if arr[2, i] == 'N':
            # print(f'arr[2]: {arr[2,i]}')
            # print(f'arr[0]: {arr[0,i]}')
            arr[2, i] = arr[0, i]
    
    new_ref = ''.join(arr[2])
    return match_ratio, new_ref

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
    new_oligo = []

    for oligo in oligos:  
        original_oligo = oligo
        
        #processed the oligo
        oligo_seq = list(oligo.seq)  
        for i in range(10, 20): 
            if i < len(oligo_seq) and oligo_seq[i] != 'U'and oligo_seq[i] != 'I':
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
        mr, new_oligo = find_match_ratio(arr)
        if current_score > best_score and mr > 0.9:  ### only when the oligo entirely matched with the read, it can be identified as the Best Oligo
            comp_str = ''.join(comp)
            start = comp_str.find('|')
            end = comp_str.rfind('|') + 1  # 包含最后一个匹配位

            query_list = list(q_aln)
            new_query = ''.join(query_list[:start]) + new_oligo + ''.join(query_list[end:])

            best_oligo = new_query
            best_score = current_score
            # print(f"original seq is: {q_aln}")
            # print(f"query start is: {query_list[:start]}")
            # print(f"oligo segment is: {new_oligo}")
            # print(f"query end is: {query_list[end:]}")
            # print(f"new query is: {new_query}")
            best_result = result
    
    if best_result is None or best_score <= 0:
        return None
    
    return best_oligo

def process_fastq_file(filename, fastq_folder, oligos, output_folder):
    if not filename.endswith(".fastq.gz"):
        return

    safe_name = sanitize_filename(filename.replace(".fastq.gz", ""))
    shared_fasta_path = os.path.join(output_folder, f"{safe_name}.fasta")
    shared_json_path = os.path.join(output_folder, f"{safe_name}.json")

    # Skip if output file already exists
    if os.path.exists(shared_fasta_path):
        print(f"Skipping {filename} (output already exists)")
        return
    
    fastq_path = os.path.join(fastq_folder, filename)
    querys, headers = read_from_fastq(fastq_path)
    
    # header_adapter = ['TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTT']  ##end repire add T at 5'
    # tail_adapter = ['GCAATACGTAACTGAACGAAGTACAGGA'] ## end repire add A at 3'
    # m13_seq = ['TTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATACGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAA']
    
    total_reads = 0
    matched_reads = 0
    all_results = []
    fasta_results = []
    results = []
    count =0

    for query, header in zip(querys, headers):
        #initial
        ref_len = 0
        reads_too_long = 0
        total_reads = total_reads + 1
        if len(query) < 1000:
            #count += 1
            fasta_results = find_oligo(query, oligos)
            #print(fasta_results)
            #breakpoint()
            if fasta_results is None:
                continue
            else:
                matched_reads += 1

                results = ''.join(fasta_results)
                
                
                fasta_records = SeqRecord(
                    Seq(results),
                    id=f"{sanitize_filename(header)}",
                    description=""
                )
                all_results.append(fasta_records)

    #print(f" {filename}: {matched_reads}/{count} matched (matched read/ 200-400bp)")
    with open(shared_fasta_path, "w") as fasta_out:
        SeqIO.write(all_results, fasta_out, "fasta")

def main():
    ### set input
    args = parse_args()
    fastq_folder = args.fastq_folder
    fasta_file = args.fasta_file
    output_folder = args.output_folder
    
    # fastq_folder = r"D:\nanopore\data\pod5\dIdU_MIX\basecalling\pass"
    # fasta_file = r"D:\nanopore\data\fasta\oligos_PCR.fasta"
    # output_folder = r"D:\nanopore\data\pod5\dIdU_MIX\basecalling\pass"
    
    os.makedirs(output_folder, exist_ok=True)

    oligos = read_from_fasta(fasta_file)
    fastq_files = [f for f in os.listdir(fastq_folder) if f.endswith(".fastq.gz")]

    start_time = time.time() 

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        future_to_file = {}
        for f in fastq_files:
            output_path = os.path.join(output_folder, os.path.splitext(os.path.splitext(f)[0])[0] + ".fasta")
            if os.path.exists(output_path):
                print(f" Skipping {f} (already processed)")
                continue
            future = executor.submit(process_fastq_file, f, fastq_folder, oligos, output_folder)
            future_to_file[future] = f

        for future in tqdm(concurrent.futures.as_completed(future_to_file), total=len(future_to_file), desc="Processing"):
            f = future_to_file[future]
            try:
                future.result()
                print(f" Finished processing: {f}")
            except Exception as e:
                print(f" Error processing {f}: {e}")

    end_time = time.time()  
    print(f"\n Total processing time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
    
#     oligos = []
#     querys = []
#     oligos = read_from_fasta(fasta_folder)

#     fastq_files = [f for f in os.listdir(fastq_folder) if f.endswith(".fastq.gz")]

#     with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
#         futures = [
#             executor.submit(process_fastq_file, f, fastq_folder, oligos, output_folder)
#             for f in fastq_files
#         ]
#         for f in concurrent.futures.as_completed(futures):
#             try:
#                 f.result()
#             except Exception as e:
#                 print(f"Error: {e}")
                
# if __name__ == "__main__":
#     main()