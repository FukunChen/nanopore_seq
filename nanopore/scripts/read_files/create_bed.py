import os
from Bio import SeqIO
import pysam
import json
from collections import defaultdict

# 设置路径
fasta_dir = "/mnt/d/nanopore/data/aligned_reference"
bam_dir = "/mnt/d/nanopore/data/pod5/nanopore/basecalling/pass/alignment"
output_dir = "/mnt/d/nanopore/data/bed"
os.makedirs(output_dir, exist_ok=True)

def create_bed(ref_fa_path, bam_path, filename):
    found_reads = 0
    # 解析 reference 中的 U/I 位置
    ref_to_ui_positions = dict()
    for record in SeqIO.parse(ref_fa_path, "fasta"):
        seq = str(record.seq).upper()
        ref_id = record.id
        ui_positions = [(i, base) for i, base in enumerate(seq) if base in ["U", "I"]]
        if ui_positions:
            ref_to_ui_positions[ref_id] = ui_positions

    print(f"Loaded {len(ref_to_ui_positions)} reference sequences with U/I.")

    # 设置输出路径
    output_json = os.path.join(output_dir, filename.replace(".bam", ".json"))
    output_bed = os.path.join(output_dir, filename.replace(".bam", ".bed"))

    # 打开 BAM 文件并遍历 primary mapped reads
    print(f"\nProcessing {filename}...")
    bam = pysam.AlignmentFile(bam_path, "rb")
    coverage_dict = defaultdict(list)
    bed_lines = []

    for aln in bam.fetch(until_eof=True):
        if aln.is_unmapped:
        #if not (aln.flag & 0x100):
            #print("f")
            continue

        read_id = aln.query_name
        ref_name = aln.reference_name
        ref_start = aln.reference_start
        
        print(f"ref_name should be:",ref_name )

        if ref_name not in ref_to_ui_positions:
            print("no ui in the ref")
            continue  # 当前 reference 中没有 U/I
        
        found_reads += 1 
        
        for rel_pos, base in ref_to_ui_positions[ref_name]:
            genome_pos = ref_start + rel_pos ###这个位置可能有错误，总感觉应该是减法
            coverage_dict[f"{ref_name}:{genome_pos+1}"].append(read_id)
            bed_lines.append(f"{ref_name}\t{genome_pos}\t{genome_pos+1}\t{base}")
        
        #breakpoint()
    
    print(f"\nDetected {found_reads} primary mapped reads with reference containing U/I.")

    bam.close()


    # 输出 JSON 和 BED 文件
    with open(output_json, "w") as f:
        json.dump(coverage_dict, f, indent=2)

    with open(output_bed, "w") as f:
        for line in bed_lines:
            f.write(line + "\n")

    print(f"Output saved: {output_json}, {output_bed}")

def main():
    #for filename in os.listdir(fasta_dir):
        #if not filename.endswith(".fasta"):
        #    continue

        #base_name = filename.replace(".fasta", ".bam")
        base_name = "aln_bed"
        #ref_fa_path = os.path.join(fasta_dir, filename)
        ref_fa_path = r"/mnt/d/nanopore/data/combined_reference.fasta"
        bam_path = r"/mnt/d/nanopore/data/bed/aln.bam"
        #bam_path = os.path.join(bam_dir, base_name)

        # if not os.path.exists(bam_path):
        #     print(f"BAM file not found for {filename}, skipping.")
        #     continue

        create_bed(ref_fa_path, bam_path, base_name)
        #breakpoint()

        print("\nAll finished!!!")

if __name__ == "__main__":
    main()