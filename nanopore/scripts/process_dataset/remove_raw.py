import json
import sys
from pathlib import Path
from tqdm import tqdm

def process_one_file(in_path: Path):
    """Remove 'raw_signal' from each JSON line, streaming line-by-line."""
    out_path = in_path.with_suffix(in_path.suffix + ".clean.jsonl")

    # 先数一遍行数（用于进度条总数）
    try:
        with in_path.open("r", encoding="utf-8") as f:
            total = sum(1 for _ in f)
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to count lines in {in_path}: {e}\n")
        total = None

    with in_path.open("r", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
        for line in tqdm(fin, total=total, desc=f"Processing {in_path.name}"):
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError as e:
                sys.stderr.write(f"[WARN] JSON parse error in {in_path}: {e}\n")
                continue
            rec.pop("raw_signal", None)  # 删除 raw_signal
            fout.write(json.dumps(rec, ensure_ascii=False) + "\n")

def main():
    input_dir = Path("/mnt/d/nanopore_seq/nanopore/scripts/process_dataset")

    if not input_dir.exists() or not input_dir.is_dir():
        sys.stderr.write(f"[ERROR] Input dir does not exist or is not a directory: {input_dir}\n")
        sys.exit(1)

    files = list(input_dir.glob("*.jsonl"))
    if not files:
        sys.stderr.write(f"[INFO] No .jsonl files found in {input_dir}\n")
        return

    for fpath in files:
        process_one_file(fpath)

    sys.stderr.write("[INFO] All files processed.\n")

if __name__ == "__main__":
    main()