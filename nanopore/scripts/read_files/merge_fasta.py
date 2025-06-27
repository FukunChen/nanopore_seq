from Bio import SeqIO
import os

# 输入和输出路径
input_folder = r"D:\nanopore\data\aligned_reference"  # ← 修改为你的 fasta 文件夹路径
output_file = r"D:\nanopore\data\aligned_reference\combined_reference.fasta"  # ← 输出合并后的 fasta

# 初始化一个列表收集所有记录
all_records = []

# 遍历文件夹中的所有 .fasta 文件
for filename in os.listdir(input_folder):
    if filename.endswith(".fasta"):
        filepath = os.path.join(input_folder, filename)
        print(f"正在读取: {filename}")
        records = list(SeqIO.parse(filepath, "fasta"))
        all_records.extend(records)

# 写入合并的 fasta 文件
SeqIO.write(all_records, output_file, "fasta")
print(f"共写入 {len(all_records)} 条序列到 {output_file}")