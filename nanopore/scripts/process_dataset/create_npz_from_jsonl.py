import sys
import json
import numpy as np
from pathlib import Path
from tqdm import tqdm

MAX_LEN = 200  # length of the position vector

def load_bed_single(bed_path: Path):
    """
    BED with 4 tab-delimited columns: seq_id, start, end, type (I or U).
    Assumption: each seq_id appears at most once; only one site and one type.
    Returns: dict[seq_id] = (pos_1based, type)
    """
    idx = {}
    with bed_path.open("r", encoding="utf-8") as f:
        for ln, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                sys.stderr.write(f"[WARN] BED line {ln} has fewer than 4 columns. Skipped: {line}\n")
                continue
            seq_id, start_s, _end_s, typ = parts[0], parts[1], parts[2], parts[3]
            try:
                pos_1based = int(start_s)
            except ValueError:
                sys.stderr.write(f"[WARN] BED line {ln} start is not an integer. Skipped: {line}\n")
                continue
            if typ not in {"I", "U"}:
                sys.stderr.write(f"[WARN] Unexpected type '{typ}' at BED line {ln}. Skipped.\n")
                continue
            if seq_id in idx:
                # Keep the first; warn once.
                sys.stderr.write(f"[WARN] Duplicate seq_id {seq_id} in BED. Keeping the first record.\n")
                continue
            idx[seq_id] = (pos_1based, typ)
    return idx

def build_pos_vector_single(p1: int):
    """
    Map a single 1-based position into a 0/1 vector (length MAX_LEN).
    Return None if out of range (>MAX_LEN or <1).
    """
    if not (1 <= p1 <= MAX_LEN):
        return None
    vec = np.zeros(MAX_LEN, dtype=np.uint8)
    vec[p1 - 1] = 1
    return vec

def write_npz_chunk(out_prefix: Path, chunk_idx: int, records_chunk, compress=True) -> Path:
    """
    Write one NPZ chunk with an object array 'records'.
    """
    out_path = out_prefix.parent / f"{out_prefix.name}.part{chunk_idx:04d}.npz"
    arr = np.array(records_chunk, dtype=object)
    if compress:
        np.savez_compressed(out_path, records=arr)
    else:
        np.savez(out_path, records=arr)
    return out_path

def process_one_jsonl(jsonl_path: Path, bed_index: dict, chunk_size: int = 100_000, compress: bool = True):
    """
    Stream a single JSONL, add 'damage' if contig in BED and pos<=MAX_LEN,
    and write multiple NPZ chunks to avoid OOM.
    Output files: <jsonl_basename>.part0000.npz, part0001.npz, ...
    """
    out_prefix = jsonl_path.with_suffix("")  # remove .jsonl, keep name for prefix
    records_chunk = []
    chunk_idx = 0
    written = []

    # Optional: pre-count lines for better tqdm total; skip if file is huge to save IO.
    try:
        with jsonl_path.open("r", encoding="utf-8") as f:
            total_lines = sum(1 for _ in f)
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to count lines in {jsonl_path}: {e}\n")
        total_lines = None

    with jsonl_path.open("r", encoding="utf-8") as f:
        for ln, line in enumerate(tqdm(f, total=total_lines, desc=f"Processing {jsonl_path.name}"), 1):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError as e:
                sys.stderr.write(f"[WARN] JSON parse error in {jsonl_path} at line {ln}: {e}\n")
                continue

            # Attach damage if we can map the contig and position is within MAX_LEN
            try:
                contig = rec["ref_region"]["contig"]
            except Exception:
                # Keep record as-is if missing the expected field
                records_chunk.append(rec)
            else:
                entry = bed_index.get(contig)
                if entry is not None:
                    pos1, typ = entry
                    vec = build_pos_vector_single(pos1)
                    if vec is not None:
                        # Only one mark and one type per record (I or U)
                        rec["damage"] = {"pos": vec.tolist(), "type": typ}
                records_chunk.append(rec)

            # Flush a chunk to NPZ
            if len(records_chunk) >= chunk_size:
                out_path = write_npz_chunk(out_prefix, chunk_idx, records_chunk, compress=compress)
                written.append(out_path)
                records_chunk.clear()
                chunk_idx += 1

    # Flush remaining
    if records_chunk:
        out_path = write_npz_chunk(out_prefix, chunk_idx, records_chunk, compress=compress)
        written.append(out_path)

    sys.stderr.write(f"[INFO] {jsonl_path.name}: wrote {len(written)} NPZ chunks.\n")
    for p in written:
        sys.stderr.write(f"  - {p}\n")

def main():
    if len(sys.argv) < 3:
        sys.stderr.write(
            f"Usage: {sys.argv[0]} <input_dir_with_jsonl> <input.bed> [chunk_size]\n"
            f"Note : For each *.jsonl, outputs <basename>.part0000.npz, part0001.npz, ... in the same directory.\n"
            f"       NPZ contains an object array 'records' (each is the full JSON record + optional 'damage').\n"
        )
        sys.exit(1)

    input_dir = Path(sys.argv[1])
    bed_path  = Path(sys.argv[2])
    chunk_size = int(sys.argv[3]) if len(sys.argv) >= 4 else 100_000  # adjust to control memory

    if not input_dir.exists() or not input_dir.is_dir():
        sys.stderr.write(f"[ERROR] Input dir not found or not a directory: {input_dir}\n")
        sys.exit(1)
    if not bed_path.exists():
        sys.stderr.write(f"[ERROR] BED file not found: {bed_path}\n")
        sys.exit(1)

    bed_index = load_bed_single(bed_path)
    sys.stderr.write(f"[INFO] Loaded BED with {len(bed_index)} entries.\n")

    jsonl_files = sorted(input_dir.glob("*_clean.jsonl"))
    if not jsonl_files:
        sys.stderr.write(f"[INFO] No *.jsonl files found in {input_dir}\n")
        return

    for jp in jsonl_files:
        process_one_jsonl(jp, bed_index, chunk_size=chunk_size, compress=True)

    sys.stderr.write("[INFO] All files processed.\n")

if __name__ == "__main__":
    main()