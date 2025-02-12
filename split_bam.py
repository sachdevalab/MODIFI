### given a bam, extract the read alignment for each contig, and store it in a separate bam file


import pysam
import os
from multiprocessing import Pool
import sys
import argparse
from Bio import SeqIO

pbindex_bin = "/home/shuaiw/smrtlink/pbindex"

def split_bam(bam, split_bam_dir, whole_ref, threads=10, min_len=50000):
    # Ensure the output directory exists
    os.makedirs(split_bam_dir, exist_ok=True)
    bam_dir = split_bam_dir + "/bams/"
    os.makedirs(bam_dir, exist_ok=True)
    contig_dir = split_bam_dir + "/contigs/"
    os.makedirs(contig_dir, exist_ok=True)
    i = 1
    args = []
    # Loop through all contigs in the ref file us ing biopython
    for record in SeqIO.parse(whole_ref, "fasta"):
        contig = record.id
        contig_len = len(record.seq)
        if contig_len < min_len:
            continue
        contig_bam = bam_dir + contig + ".bam"
        ref = contig_dir + contig + ".fa"
        args.append((contig, ref, contig_bam, bam, whole_ref))
        i += 1
    # samfile.close()

    # Use multiprocessing to handle each contig in parallel
    with Pool(processes=threads) as pool:
        pool.starmap(handle_each_contig, args)
    print ("there are ", i, " contigs")

def handle_each_contig(contig,ref,contig_bam,bam,whole_ref):
    print (f"Processing {contig} {ref}")
    os.system(f"samtools faidx {whole_ref} {contig} > {ref}")
    os.system(f"samtools faidx {ref}")


    #  Create a new header that only includes the specific contig
    samfile = pysam.AlignmentFile(bam, "rb")
    new_header = samfile.header.to_dict().copy()
    
    new_header['SQ'] = [sq for sq in new_header['SQ'] if sq['SN'] == contig]

    ### covert the too large read group id to fit int32, and convert it to hexadecimal
    ### not working, need to fix it
    # for i, rg in enumerate(new_header['RG']):
    #     rg['ID'] = hex(i)

    contig_samfile = pysam.AlignmentFile(contig_bam, "wb", header=new_header)
    for read in samfile.fetch(contig):
        read.reference_id = 0
        contig_samfile.write(read)
    contig_samfile.close()
    samfile.close()
    ## index the bam file
    os.system(f"samtools index {contig_bam}")
    os.system(f"{pbindex_bin} {contig_bam}")

def main():
    parser = argparse.ArgumentParser(description="Split BAM file by reference.")
    parser.add_argument("--bam", help="Input BAM file")
    parser.add_argument("--whole_ref", help="Reference name")
    parser.add_argument("--work_dir", help="Directory to save split BAM files")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use (default: 1)")
    parser.add_argument("--min_len", type=int, default=50000, help="Minimum length of reads to include (default: 1000)")

    args = parser.parse_args()

    split_bam(args.bam, args.work_dir, args.whole_ref, args.threads, args.min_len)


if __name__ == "__main__":
    # bam = "/home/shuaiw/methylation/data/borg/all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
    # whole_ref="/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    # split_bam_dir = "/home/shuaiw/methylation/data/borg/split_bam_dir3/"
    # split_bam(bam, split_bam_dir)
    # bam = sys.argv[1]
    # ref = sys.argv[2]
    # split_bam_dir = sys.argv[3]
    # threads = int(sys.argv[4])
    # min_len = int(sys.argv[5])
    # split_bam(bam, split_bam_dir, ref, threads, min_len)
    main()

