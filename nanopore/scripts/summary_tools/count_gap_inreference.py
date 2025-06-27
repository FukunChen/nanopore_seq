from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
fasta_path = r"D:\nanopore\data\aligned_reference\fastq_runid_a61744f1-cde4-4b61-bc16-2d903bf0f709_0_0.fasta"


gap_counts = []
gap_ratios = []

with open(fasta_path, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq = str(record.seq)
        total_len = len(seq)
        gap_count = seq.count('-')
        gap_ratio = gap_count / total_len if total_len > 0 else 0
        gap_counts.append(gap_count)
        gap_ratios.append(gap_ratio)

# 计算 CDF
gap_counts_sorted = np.sort(gap_counts)
cdf_counts = np.arange(1, len(gap_counts_sorted)+1) / len(gap_counts_sorted)

gap_ratios_sorted = np.sort(gap_ratios)
cdf_ratios = np.arange(1, len(gap_ratios_sorted)+1) / len(gap_ratios_sorted)

# 绘图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Gap 数量 CDF
ax1.plot(gap_counts_sorted, cdf_counts, marker='.', linestyle='-')
ax1.set_xlabel("Gap Count per Read")
ax1.set_ylabel("CDF")
ax1.set_title("Gap Count CDF")
ax1.grid(True)

# Gap 占比 CDF
ax2.plot(gap_ratios_sorted, cdf_ratios, marker='.', linestyle='-')
ax2.set_xlabel("Gap Ratio (gap count / total length)")
ax2.set_ylabel("CDF")
ax2.set_title("Gap Ratio CDF")
ax2.grid(True)

plt.tight_layout()
plt.show()