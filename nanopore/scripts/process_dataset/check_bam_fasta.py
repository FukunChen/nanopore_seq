import pysam
from Bio import SeqIO
import sys

def count_mapped_read(bampath):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    count = sum(1 for read in bamfile.fetch(until_eof=True))
    bamfile.close()
    return count

def count_U_I_reads(fasta_path):
    u_count = 0
    i_count = 0
    total = 0

    for record in SeqIO.parse(fasta_path,"fasta"):
        seq = str(record.seq).upper()
        total += 1
        if "U" in seq:
            u_count += 1
        if "I" in seq:
            i_count += 1

    return total, u_count, i_count

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(1)
    
    bam_path = sys.argv[1]
    fasta_path = sys.argv[2]

    print(f"Checking Bam file...")
    mapped_count = count_mapped_read(bam_path)
    print(f"Mapped reads in BAM file are {mapped_count} / mapped prima")

    print(f"\nChecking FASTA file...")
    total, u_count, i_count = count_U_I_reads(fasta_path)
    print(f"Total reads in fasta {total}")
    print(f"Total read of U {u_count}")
    print(f"Total reads of I {i_count}")
