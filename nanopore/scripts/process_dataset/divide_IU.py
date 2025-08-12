import sys
from pathlib import Path
import numpy as np
from tqdm import tqdm

DEFAULT_CHUNK = 5000
MAX_QTS_LEN = 200

def save_chunk(out_dir: Path, tag: str, part_idx: int, recs, compress=True) -> Path:
    """
    写一个 NPZ 分片：
      <out_dir>/<tag>.partXXXX.npz   # tag 为 'I' 或 'U'
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{tag}.part{part_idx:04d}.npz"
    arr = np.array(recs, dtype=object)
    if compress:
        np.savez_compressed(out_path, records=arr)
    else:
        np.savez(out_path, records=arr)
    return out_path

def split_npz_by_type(src_dir: Path, out_dir: Path, chunk_size: int = DEFAULT_CHUNK):
    npz_files = sorted(src_dir.rglob("*.npz"))
    if not npz_files:
        print(f"[INFO] No .npz found under {src_dir}")
        return

    buf_I, buf_U = [], []
    part_I = part_U = 0
    written_I = written_U = 0
    seen = 0
    skipped = 0

    for npz_path in tqdm(npz_files, desc="Scanning NPZ"):
        try:
            data = np.load(npz_path, allow_pickle=True)
        except Exception as e:
            print(f"[WARN] Failed to load {npz_path}: {e}")
            continue

        if "records" not in data:
            print(f"[WARN] 'records' not found in {npz_path}, skipped.")
            continue

        recs = data["records"]  # object array
        for rec in recs:
            seen += 1
            try:
                typ = rec.get("damage", {}).get("type", None)
            except AttributeError:
                # rec 不是 dict
                skipped += 1
                continue
            if "query_to_signal" in rec and isinstance(rec["query_to_signal"], list):
                rec["query_to_signal"] = rec["query_to_signal"][:MAX_QTS_LEN]
            if typ == "I":
                buf_I.append(rec)
                if len(buf_I) >= chunk_size:
                    save_chunk(out_dir, "I", part_I, buf_I, compress=True)
                    buf_I.clear()
                    part_I += 1
                    written_I += 1
            elif typ == "U":
                buf_U.append(rec)
                if len(buf_U) >= chunk_size:
                    save_chunk(out_dir, "U", part_U, buf_U, compress=True)
                    buf_U.clear()
                    part_U += 1
                    written_U += 1
            else:
                # 没有 damage 或未知类型
                skipped += 1

        # 及时释放
        try:
            data.close()
        except Exception:
            pass

    # 收尾：写出剩余不足 chunk 的
    if buf_I:
        save_chunk(out_dir, "I", part_I, buf_I, compress=True)
        written_I += 1
    if buf_U:
        save_chunk(out_dir, "U", part_U, buf_U, compress=True)
        written_U += 1

    print(f"[INFO] Done. seen={seen}, skipped={skipped}, wrote I={written_I} chunks, U={written_U} chunks")
    print(f"[INFO] Output dir: {out_dir}")

def main():
    if len(sys.argv) < 3:
        print(
            f"Usage: {sys.argv[0]} <src_npz_dir> <out_dir> [chunk_size]\n"
            f"Desc : Recursively read NPZs under <src_npz_dir>, split records by damage.type (I/U),\n"
            f"       and write NPZ chunks (default 500 per chunk) into <out_dir>."
        )
        sys.exit(1)

    src_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    chunk_size = int(sys.argv[3]) if len(sys.argv) >= 4 else DEFAULT_CHUNK

    if not src_dir.is_dir():
        print(f"[ERROR] Not a directory: {src_dir}")
        sys.exit(1)

    split_npz_by_type(src_dir, out_dir, chunk_size)

if __name__ == "__main__":
    main()