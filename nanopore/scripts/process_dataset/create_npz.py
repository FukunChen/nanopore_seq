import pysam
import pod5
import numpy as np
from tqdm import tqdm
from collections import defaultdict
import argparse


def main():
    parser = argparse.ArgumentParser(description="Extract U/I reads from BAM and save as NPZ.")
    parser.add_argument("--bam_path", required=True, help="Path to the BAM file.")
    parser.add_argument("--bed_path", required=True, help="Path to the BED file containing U/I annotations.")
    parser.add_argument("--pod5_path", required=True, help="Path to the POD5 file.")
    args = parser.parse_args()
    
    bam_path = args.bam_path
    bed_path = args.bed_path
    pod5_path = args.pod5_path
    output_u_npz = "calls_u.npz"
    output_i_npz = "calls_i.npz"

    # === 打开 BAM 读取器，收集 read_id 对应的损伤类型（U/I） ===
    bam_fh = pysam.AlignmentFile(bam_path, "rb")
    read_damage_type = defaultdict(set)

    print("读取 BED 中的损伤信息...")
    with open(bed_path) as f_bed:
        for line in f_bed:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue

            ctg, start, end, damage_type = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if damage_type not in {"U", "I"}:
                continue

            for read in bam_fh.fetch(ctg, start, end):
                if not read.is_secondary and not read.is_supplementary:
                    read_damage_type[read.query_name].add(damage_type)

    bam_fh.close()

    print(f"总计检测到 {len(read_damage_type)} 个包含 U/I 的 reads。")

    # === 重新打开 BAM 和 POD5，准备提取数据 ===
    bam_fh = pysam.AlignmentFile(bam_path, "rb")
    pod5_reader = pod5.Reader(pod5_path)

    data_u = {"signal": [], "seq": [], "label": [], "move": []}
    data_i = {"signal": [], "seq": [], "label": [], "move": []}
    processed = set()

    print("正在提取 BAM 中对应的 primary reads...")
    for read in tqdm(bam_fh.fetch(until_eof=True), desc="提取中"):
        read_id = read.query_name

        if read_id in processed:
            continue
        if read_id not in read_damage_type:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.query_sequence is None:
            continue

        try:
            pod5_read = next(pod5_reader.reads(selection=[read_id]))
        except RuntimeError:
            continue

        try:
            move_table = read.get_tag("mv")
        except KeyError:
            continue

        sequence = read.query_sequence
        signal = pod5_read.signal

        if "U" in read_damage_type[read_id]:
            label_u = [1 if base == "U" else 0 for base in sequence]
            data_u["seq"].append(list(sequence))
            data_u["signal"].append(signal)
            data_u["label"].append(label_u)
            data_u["move"].append(move_table)

        if "I" in read_damage_type[read_id]:
            label_i = [1 if base == "I" else 0 for base in sequence]
            data_i["seq"].append(list(sequence))
            data_i["signal"].append(signal)
            data_i["label"].append(label_i)
            data_i["move"].append(move_table)

        processed.add(read_id)

    # === 保存为 .npz 文件 ===
    print("正在保存 npz 文件...")

    np.savez(output_u_npz,
            signal=np.array(data_u["signal"], dtype=object),
            seq=np.array(data_u["seq"], dtype=object),
            label=np.array(data_u["label"], dtype=object),
            move=np.array(data_u["move"], dtype=object))

    np.savez(output_i_npz,
            signal=np.array(data_i["signal"], dtype=object),
            seq=np.array(data_i["seq"], dtype=object),
            label=np.array(data_i["label"], dtype=object),
            move=np.array(data_i["move"], dtype=object))

    print(f"✅ calls_u.npz: 包含 {len(data_u['seq'])} 条 U read")
    print(f"✅ calls_i.npz: 包含 {len(data_i['seq'])} 条 I read")
    print("🎉 全部完成！")

if __name__ == "__main__":
    main()