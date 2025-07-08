import pysam
from collections import defaultdict
import pod5

###filte primary read
def read_is_primary(aln):
    return not aln.is_secondary and not aln.is_supplementary

def check_bed_coverage(bam_path, bed_path, pod5_path, output_file):
    bam_fh = pysam.AlignmentFile(bam_path, "rb")
    
    if not bam_fh.has_index():
        #print(f"Indexing BAM file: {bam_path}")
        pysam.index(bam_path)
        bam_fh = pysam.AlignmentFile(bam_path, "rb")

    print(f"Loading read IDs from {pod5_path}")
    pod5_read_ids = set()
    with pod5.Reader(pod5_path) as reader:
        for read in reader.reads():
            pod5_read_ids.add(read.read_id)

    uncovered_regions = []
    passed_regions = []

    with open(bed_path) as bed_file:
        for line in bed_file:
            if line.startswith("#") or line.strip() == "":
                continue

            fields = line.strip().split()
            ctg = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            all_reads = list(bam_fh.fetch(ctg, start, end))
            if not all_reads:
                uncovered_regions.append(f"{ctg}:{start}-{end} (no reads)")
                continue

            # 筛选 primary read
            primary_reads = [r for r in all_reads if read_is_primary(r)]
            if not primary_reads:
                uncovered_regions.append(f"{ctg}:{start}-{end} (no primary)")
                continue

            # 检查这些 primary read 是否出现在 pod5 中
            found = False
            for read in primary_reads:
                if read.query_name in pod5_read_ids:
                    found = True
                    break

            if found:
                passed_regions.append((ctg, start, end))
            else:
                uncovered_regions.append(f"{ctg}:{start}-{end} (primary read(s) not in pod5)")

    print("\n Regions covered and with valid pod5 reads: saved to passed_regions.bed")
    with open(output_file, "w") as f_out:
        for ctg, start, end in passed_regions:
            f_out.write(f"{ctg}\t{start}\t{end}\n")

    if uncovered_regions:
        print("\n Regions without usable reads:")
        #for region in uncovered_regions:
            #print(f"  {region}")
    else:
        print("\n All regions are covered by valid pod5 + BAM primary reads.")

def main(args):
    bam_path = args.bam_path
    bed_path = args.bed_path
    pod5_path = args.pod5_path
    output_file = args.output_file
    
    check_bed_coverage(bam_path, bed_path, pod5_path, output_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter BED regions based on BAM and POD5 read coverage")

    parser.add_argument("--bam_path", type=str, required=True, help="Path to input BAM file")
    parser.add_argument("--bed_path", type=str, required=True, help="Path to input BED file")
    parser.add_argument("--pod5_path", type=str, required=True, help="Path to input POD5 file")
    parser.add_argument("--output_file", type=str, required=True, help="Path to output filtered BED file")


