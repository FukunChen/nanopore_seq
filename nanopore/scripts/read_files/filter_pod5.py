import os
from Bio import SeqIO
from pod5.reader import Reader
from pod5.writer import Writer

n = 1

filtered_fasta_dir = r"D:\nanopore\data\aligned_reference\filted_ref"
input_pod5_file = rf"D:\nanopore\data\pod5\nanopore\FBC86537_skip_a61744f1_40a78090_{n}.pod5" # 假设你整合了全部 read 的 POD5
output_pod5_dir = rf"D:\nanopore\data\filtered_pod5_{n}"

os.makedirs(output_pod5_dir, exist_ok=True)

# 遍历所有 filtered fasta 文件
for fasta_file in os.listdir(filtered_fasta_dir):
    if not fasta_file.endswith(".filtered.fasta"):
        continue

    fasta_path = os.path.join(filtered_fasta_dir, fasta_file)

    # 提取 read_ids
    read_ids_raw = [record.id for record in SeqIO.parse(fasta_path, "fasta")]
    read_ids = {rid.strip().lower() for rid in read_ids_raw}
    #read_ids = ["a61744f1-cde4-4b61-bc16-2d903bf0f709"]

    if not read_ids:
        print(f"{fasta_file}: 没有 read 被保留，跳过")
        continue

    output_pod5_path = os.path.join(output_pod5_dir, fasta_file.replace(".filtered.fasta", ".filtered.pod5"))

    # 从原始 pod5 中筛选
    with Reader(input_pod5_file) as reader, Writer(output_pod5_path) as writer:
        count = 0
        for read in reader.reads():
            #print(f"read id of pod5 is {read.read_id}")
            if str(read.read_id).strip().lower() in read_ids:
                writer.add_read(read.to_read())
                count += 1

    print(f"{fasta_file}: 提取 {count} 条 read 到 {output_pod5_path}")