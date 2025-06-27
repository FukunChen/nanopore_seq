import pod5
import matplotlib.pyplot as plt

# 设置 pod5 文件路径和 read_id
pod5_path = r"D:\nanopore\data\pod5\nanopore\FBC86537_skip_a61744f1_40a78090_0.pod5"
read_id = "0ad7abae-38c5-423f-a3a9-4efff6753a8a"

# 打开 pod5 文件并读取目标 read 的信号
with pod5.Reader(pod5_path) as reader:
    for read in reader.reads(selection=[read_id]):
        signal = read.signal  # numpy array
        print(f"Signal length: {len(signal)}")

        # 绘图
        plt.figure(figsize=(25, 4))
        plt.plot(range(len(signal)), signal, linewidth = 0.5)
        plt.xlim(0, len(signal))
        plt.title(f"Raw Signal for Read {read_id}")
        plt.xlabel("Sample Index")
        plt.ylabel("Current Signal (ADC units)")
        plt.tight_layout()
        plt.show()