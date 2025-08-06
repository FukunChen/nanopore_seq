import pod5
import matplotlib.pyplot as plt
import pysam
import numpy as np
from pathlib import Path
import plotnine as p9
import pandas as pd
import matplotlib.patches as patches
import json
from tqdm import tqdm 

import remora
from remora import io

pod5_path = r"/mnt/d/nanopore/data/pod5/dIdU_MIX/merged_didumix.pod5"
bam_path = r"/mnt/d/nanopore/calls_didumix.bam/didumix_onlydamage.bam"


def compute_base_space_sig_coords(seq_to_sig_map):
    """Compute coordinates for signal points, interpolating signal assigned to
    each base linearly through the span of each covered base.
    """
    return np.interp(
        np.arange(seq_to_sig_map[-1] - seq_to_sig_map[0]),
        seq_to_sig_map,
        np.arange(seq_to_sig_map.size),
    )

output_json = "ref_views_didumix.json"

pod5_dr = pod5.DatasetReader(pod5_path)
bam_fh = io.ReadIndexedBam(bam_path)

ref_views = []

# 先提取 BAM 中所有 read_id
bam_read_ids = list(bam_fh.read_ids)
print(f"📌 BAM contains {len(bam_read_ids)} reads")

for read_id in tqdm(bam_read_ids, total=len(bam_read_ids), desc="Processing reads"):
    bam_read = bam_fh.get_first_alignment(read_id)
    if bam_read is None:
        continue  # unmapped

    try:
        pod5_read = pod5_dr.get_read(read_id)
    except KeyError:
        # 这个 read_id 不在 pod5 中
        continue

    io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)  ###io_read 是我的数据
    #print(io_read)

    if io_read.ref_reg is None:
        continue  # unmapped read

    region = io_read.ref_reg.adjust(
        start_adjust=0,
        end_adjust=200 - io_read.ref_reg.len
    ) ###只是返回了一个Refregion格式的region

    ref_view = io_read.extract_ref_reg(region) 

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
    ref_views.append(ref_view_data)
    #breakpoint()

# 保存 JSON
with open(output_json, "w") as f:
    json.dump(ref_views, f, indent=2)

print(f"✅ Saved {len(ref_views)} ref_views to {output_json}")