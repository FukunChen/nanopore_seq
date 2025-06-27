import os
from Bio import SeqIO

# 设置输入输出文件夹路径
fasta_folder = r"D:\nanopore\data\aligned_reference"
output_folder = r"D:\nanopore\data\aligned_reference\filted_ref"

# 创建输出文件夹（如果不存在）
os.makedirs(output_folder, exist_ok=True)

# 筛选条件
MAX_GAP_COUNT = 6
MAX_GAP_RATIO = 0.05

# 遍历所有 fasta 文件
for filename in os.listdir(fasta_folder):
    if filename.endswith(".fasta"):
        input_path = os.path.join(fasta_folder, filename)
        output_filename = filename.replace(".fasta", ".filtered.fasta")
        output_path = os.path.join(output_folder, output_filename)

        filtered_records = []
        total_records = 0
        with open(input_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                total_records += 1
                seq = str(record.seq)
                total_len = len(seq)
                gap_count = seq.count('-')
                gap_ratio = gap_count / total_len if total_len > 0 else 0

                if gap_count <= MAX_GAP_COUNT and gap_ratio <= MAX_GAP_RATIO:
                    filtered_records.append(record)

        print(f"{filename}: 原始 read 数量 = {total_records}")
        
        # 写入筛选后的 reads
        if filtered_records:
            SeqIO.write(filtered_records, output_path, "fasta")
            print(f"{output_filename}: 筛选出 {len(filtered_records)} 条 read")
        else:
            print(f"{output_filename}: 无符合条件的 read，未保存")