from Bio import SeqIO
import matplotlib.pyplot as plt

fastq_path = r"D:\nanopore\esox\demo\dev\fastq\demo.fastq"  # 替换为你的实际路径

read_lengths = []

with open(fastq_path, "r") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        read_lengths.append(len(record.seq))

# 绘图
plt.figure(figsize=(10, 5))
plt.hist(read_lengths, bins=50, color='skyblue', edgecolor='black')

# 添加竖线标记 250 bp
plt.axvline(x=250, color='red', linestyle='--', linewidth=2, label='250 bp')
plt.text(250 + 5, plt.ylim()[1]*0.9, '250 bp', color='red')

plt.xlabel("Read Length (bp)")
plt.ylabel("Frequency")
plt.title("Distribution of Read Lengths in FASTQ File")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()