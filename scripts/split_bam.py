### given a bam, extract the read alignment for each contig, and store it in a separate bam file


import pysam
import os
from multiprocessing import Pool
import sys
import argparse
from Bio import SeqIO
import numpy as np

MAP_Q = 20
MAX_DEPTH = 500

# pbindex_bin = "/home/shuaiw/smrtlink/pbindex"

def split_bam(bam, split_bam_dir, whole_ref, pbindex_bin, threads=10, min_len=50000, max_NM=3, min_dp=0, min_iden=0.97):
    # Ensure the output directory exists
    os.makedirs(split_bam_dir, exist_ok=True)
    bam_dir = os.path.join(split_bam_dir, "bams")
    os.makedirs(bam_dir, exist_ok=True)
    contig_dir = os.path.join(split_bam_dir, "contigs")
    os.makedirs(contig_dir, exist_ok=True)
    i = 1
    args = []
    # Loop through all contigs in the ref file us ing biopython
    for record in SeqIO.parse(whole_ref, "fasta"):
        contig = record.id
        contig_len = len(record.seq)
        if contig_len < min_len:
            continue
        contig_bam = os.path.join(bam_dir, contig + ".bam")
        print (f"Processing {contig} {contig_bam}")
        ref = os.path.join(contig_dir, contig + ".fa")
        args.append((contig, contig_len, ref, contig_bam, bam, whole_ref, pbindex_bin, max_NM, MAP_Q, MAX_DEPTH, min_dp, min_iden))
        i += 1
    # samfile.close()
    ctg_depth_dict = {}
    # Use multiprocessing to handle each contig in parallel
    with Pool(processes=threads) as pool:
        results = pool.starmap(handle_each_contig, args)
    for contig, mean_depth, contig_len in results:
        ctg_depth_dict[contig] = [mean_depth, contig_len]

    return ctg_depth_dict

def count_cigar(read, long_cutoff = 50):
    gap_num = 0
    map_len = 0
    long_gap = 0
    

    for ci in read.cigar:
        if ci[0] == 0:
            map_len += ci[1]
        elif ci[0] == 1 or ci[0] == 2:
            gap_num += ci[1]
            map_len += ci[1]
            if ci[1] > long_cutoff:
                long_gap += ci[1]
    return gap_num, map_len, long_gap


def parse_cigar_operations(read, cigar_tuples, contig_len, skip_indel_len=50):
    """
    Parse CIGAR operations to count matches, mismatches, insertions, and deletions
    
    CIGAR operations:
    0: M (match/mismatch)
    1: I (insertion)
    2: D (deletion)
    4: S (soft clipping)
    5: H (hard clipping)
    7: = (sequence match)
    8: X (sequence mismatch)
    """
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    soft_clip = 0
    hard_clip = 0

    i = 0
    cigar_num =len(cigar_tuples)
    
    for op, length in cigar_tuples:
        if op == 0:  # M (match/mismatch) - need NM tag to distinguish
            matches += length  # We'll adjust this with NM tag
        elif op == 1:  # I (insertion)
            if length > skip_indel_len:  ## skip long indels which might derive from HGT
                continue
            insertions += length
        elif op == 2:  # D (deletion)
            if length > skip_indel_len:  ## skip long indels which might derive from HGT
                continue
            deletions += length
        elif op == 7:  # = (sequence match)
            matches += length
        elif op == 8:  # X (sequence mismatch)
            mismatches += length
        elif op == 4:  # S (soft clipping)
            if i == 0 and read.reference_start == 0:   ### ignore clip at contig start, might led by circular
                continue
            elif i == cigar_num - 1 and read.reference_end == contig_len:  ### ignore clip at contig end, might led by circular
                continue
            soft_clip += length
        elif op == 5:  # H (hard clipping)
            if i == 0 and read.reference_start == 0:   ### ignore clip at contig start, might led by circular
                continue
            elif i == cigar_num - 1 and read.reference_end == contig_len:  ### ignore clip at contig end, might led by circular
                continue
            hard_clip += length
        i += 1
    clip = soft_clip + hard_clip
    # print (soft_clip, hard_clip, clip, matches)
    return matches, mismatches, insertions, deletions, clip

def calculate_identities(read, contig_len):
    """
    Calculate three types of identity based on CIGAR string and NM tag
    """
    cigar_tuples = read.cigartuples
    
    if not cigar_tuples:
        return 0
    
    # Parse CIGAR operations
    matches, mismatches, insertions, deletions, clip = parse_cigar_operations(read, cigar_tuples, contig_len)

    
    # 3. Alignment identity (fraction of aligned columns that are matches)
    total_aligned_columns = matches + mismatches + insertions + deletions + clip
    alignment_identity = matches / total_aligned_columns if total_aligned_columns > 0 else 0
    clip_ratio = clip / total_aligned_columns if total_aligned_columns > 0 else 0

    return alignment_identity, clip_ratio


def test_read(read, contig_len, max_NM, q, min_iden):
    if read.is_unmapped:
        return False
    if read.mapping_quality < q:
        return False
    if not read.has_tag('NM'):
        print (f"Read {read.query_name} has no NM tag")
        return False
    # if read.get_tag("NM") > max_NM:
    #     return False
    alignment_identity, clip_ratio = calculate_identities(read, contig_len)
    if alignment_identity < min_iden:
        return False
    # if clip_ratio > 0.1:
    #     return False
    return True

def handle_each_contig(contig,contig_len,ref,contig_bam,bam,whole_ref, pbindex_bin, max_NM, q=20, max_depth=500, min_dp=0, min_iden=0.97):
    print (f"Processing {contig} {ref}")

    #  Create a new header that only includes the specific contig
    samfile = pysam.AlignmentFile(bam, "rb")
    new_header = samfile.header.to_dict().copy()
    
    new_header['SQ'] = [sq for sq in new_header['SQ'] if sq['SN'] == contig]

    ## first count the number of reads that are valid
    raw_read_number = 0
    total_bases = 0
    for read in samfile.fetch(contig):
        raw_read_number += 1
        passed = test_read(read, contig_len, max_NM, q, min_iden)
        if not passed:
            continue
        
        total_bases += read.query_alignment_length

    mean_depth = total_bases / contig_len if contig_len > 0 else 0
    mean_depth = round(mean_depth, 2)
    print(f"Mean depth for {contig}: {mean_depth:.2f}")
    if mean_depth < min_dp:
        print(f"Mean depth {mean_depth} is less than minimum depth {min_dp}. Skipping {contig}.")
        return contig, mean_depth, contig_len
    
    if mean_depth == 0:
        downsample_rate = 1
    else:
        downsample_rate = min(float(max_depth) / mean_depth, 1.0)
    print (f"Sample rate is {downsample_rate*100}%.")


    contig_samfile = pysam.AlignmentFile(contig_bam, "wb", header=new_header)
    valid_num = 0
    for read in samfile.fetch(contig):
        read.reference_id = 0

        ## downsample the reads
        if np.random.rand() > downsample_rate:
            continue
        
        passed = test_read(read, contig_len, max_NM, q, min_iden)
        if not passed:
            continue

        contig_samfile.write(read)
        valid_num += 1
    print (f"Finished {contig}, {valid_num} reads are written, raw read number is {raw_read_number}.")

    contig_samfile.close()
    samfile.close()
    ## index the bam file
    os.system(f"samtools index {contig_bam}")
    os.system(f"{pbindex_bin} {contig_bam}")
    os.system(f"samtools faidx {whole_ref} {contig} > {ref}")
    os.system(f"samtools faidx {ref}")
    return contig, mean_depth, contig_len

def main():
    parser = argparse.ArgumentParser(description="Split BAM file by reference.")
    parser.add_argument("--bam", help="Input BAM file")
    parser.add_argument("--whole_ref", help="Reference name")
    parser.add_argument("--work_dir", help="Directory to save split BAM files")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use (default: 1)")
    parser.add_argument("--min_len", type=int, default=50000, help="Minimum length of reads to include (default: 1000)")
    parser.add_argument("--max_NM", type=int, help="<int> Max mismatch number in CCS reads.", default=3, metavar="\b")
    parser.add_argument("--min_iden", type=float, default=0.97,
                        help="Minimum identity allowed for read alignment.")

    args = parser.parse_args()

    split_bam(args.bam, args.work_dir, args.whole_ref, args.threads, args.min_len, args.max_NM, args.min_iden)


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

