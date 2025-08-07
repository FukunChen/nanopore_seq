import pod5
import pysam
import numpy as np
import json
from tqdm import tqdm
from remora import io

# 输入路径
pod5_path = "/mnt/d/nanopore/data/pod5/dIdU_MIX/merged_didumix.pod5"
bam_path = "/mnt/d/nanopore/calls_didumix.bam/didumix_onlydamage.bam"

# 打开文件
pod5_dr = pod5.DatasetReader(pod5_path)
bam_fh = io.ReadIndexedBam(bam_path)

# 获取所有 BAM read_id
bam_read_ids = list(bam_fh.read_ids)
total_reads = len(bam_read_ids)
print(f"📌 BAM contains {total_reads} reads")

# 每 5% 一批（共 20 批）
batch_size = total_reads // 20

# 遍历每个 batch
for batch_idx in range(20):
    start = batch_idx * batch_size
    end = (batch_idx + 1) * batch_size if batch_idx < 19 else total_reads
    batch_read_ids = bam_read_ids[start:end]

    output_jsonl = f"ref_views_didumix_batch_{batch_idx+1}.jsonl"
    print(f"\n🚀 Processing batch {batch_idx+1}/20: {len(batch_read_ids)} reads → {output_jsonl}")

    with open(output_jsonl, "w") as out_f:
        for read_id in tqdm(batch_read_ids,
                            total=len(batch_read_ids),
                            desc=f"[Batch {batch_idx+1:02d}/20]",
                            ncols=100):
            # 获取 BAM read
            bam_read = bam_fh.get_first_alignment(read_id)
            if bam_read is None:
                continue

            # 获取 pod5 read
            try:
                pod5_read = pod5_dr.get_read(read_id)
            except KeyError:
                continue

            try:
                # 构造 io_read
                io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)
                if io_read.ref_reg is None:
                    continue

                # 调整 reference window 区间
                region = io_read.ref_reg.adjust(start_adjust=0, end_adjust=200 - io_read.ref_reg.len)
                ref_view = io_read.extract_ref_reg(region)

                # 构造结果字典
                ref_view_data = {
                    "raw_signal": io_read.dacs.tolist(),
                    "query_to_signal": io_read.query_to_signal.tolist(),
                    "read_id": ref_view.read_id,
                    "sequence": ref_view.seq,
                    "signal": ref_view.norm_signal.tolist(),
                    "seq_to_sig_map": ref_view.seq_to_sig_map.tolist(),
                    "ref_region": {
                        "contig": ref_view.ref_reg.ctg,
                        "strand": ref_view.ref_reg.strand,
                        "start": ref_view.ref_reg.start,
                        "end": ref_view.ref_reg.end,
                    },
                    "sig_start": int(ref_view.sig_start),
                }

                # 写入一行 JSON
                out_f.write(json.dumps(ref_view_data) + "\n")

            except Exception:
                continue  # 出错的 read 跳过

    print(f"✅ Finished batch {batch_idx+1}, saved to {output_jsonl}")