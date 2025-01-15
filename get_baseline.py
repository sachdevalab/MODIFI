import pysam
import numpy as np
from collections import defaultdict
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import random
import pickle
import os

def read_subread_bam(bam_file):
    """
    see this for fi and ri tag : 
    https://pacbiofileformats.readthedocs.io/en/12.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-kinetic-information
    """
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    contig_forward_dict = defaultdict(list)
    contig_reverse_dict = defaultdict(list)
    print ("start reading")
    contig = "NC_000001.11"
    i = 0
    for read in samfile.fetch(contig, 1, 100000):

        IPD_info = read.get_tag("ip")

        # Get aligned positions on the reference for each read base
        aligned_pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is not None:
                if IPD_info[query_pos] != 0:
                    
                    if read.is_reverse == False:
                        contig_reverse_dict[ref_pos].append(IPD_info[query_pos])
                    else:
                        ## get the reverse index for query position
                        query_pos = len(IPD_info) - query_pos - 1
                        contig_forward_dict[ref_pos].append(IPD_info[query_pos])
        i += 1
        if i % 100 == 0:
            print (i)

    seq = seq_dict[contig]
    observed_IPD_list = get_IPD_list(len(seq), contig_forward_dict)
    observed_IPD_reverse_list = get_IPD_list(len(seq), contig_reverse_dict)
    kmer_baseline_dict = prepare_data(seq, observed_IPD_list)
    for kmer in kmer_baseline_dict:
        print (kmer, len(kmer_baseline_dict[kmer]), np.mean(kmer_baseline_dict[kmer]))
        # break

    # return observed_IPD_list, observed_IPD_reverse_list

def read_hifi_bam(bam_file):
    """
    see this for fi and ri tag : 
    https://pacbiofileformats.readthedocs.io/en/12.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-kinetic-information
    """
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    contig_forward_dict = defaultdict(list)
    contig_reverse_dict = defaultdict(list)
    for read in samfile.fetch("SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C", 1, 1000):
        ## if NM tag larger than 3, skip
        if read.get_tag("NM") > 3:
            continue

        forward_IPD_info = read.get_tag("fi")[::-1]   ## weired, why need to reverse
        reverse_IPD_info = read.get_tag("ri")

        # Get aligned positions on the reference for each read base
        aligned_pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is not None:
                if forward_IPD_info[query_pos] != 0:
                    contig_forward_dict[ref_pos].append(forward_IPD_info[query_pos])
                if reverse_IPD_info[query_pos] != 0:
                    contig_reverse_dict[ref_pos].append(reverse_IPD_info[query_pos])

    observed_IPD_list = get_IPD_list(1000, contig_forward_dict)
    observed_IPD_reverse_list = get_IPD_list(1000, contig_reverse_dict)
    return observed_IPD_list, observed_IPD_reverse_list


def get_IPD_list(contig_length, contig_dict):
    average_IPD = []
    for pos in range(contig_length):
        if pos not in contig_dict:
            average_IPD.append(0)
        else:
            average_IPD.append(round(np.mean(contig_dict[pos]), 1))
    return average_IPD
             
def extract_context(fasta):
    ## load the fasta using biopython
    print ("loading fasta")
    seq_dict = {}
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        # print(record.id)
        if record.id != "NC_000001.11":
            continue
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        raw_seq = seq.upper()
        seq_dict[record.id] = raw_seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict

def prepare_data(seq, control_list, up=8, down=4):
    # print(seq[:100])
    # seq, control_list = seq[:1000], control_list[:1000]
    save_kmer = {}  ## {pos:kmer}
    kmer_baseline_dict = defaultdict(list)
    # y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        if seq[i] == 'N':
            continue
        kmer = seq[i-up:i+down]
        ipd = control_list[i]
        if ipd == 0:
            continue
        kmer_baseline_dict[kmer].append(ipd)
        save_kmer[i] = kmer    
    return kmer_baseline_dict, save_kmer



def extract_ipd_ratio(file_path):
    ## open it using pandas
    contig = "NC_000001.11"
    contig_forward_dict = defaultdict(list)
    contig_reverse_dict = defaultdict(list)
    ipd_sum_for_control = {}
    ipd_sum_for_reverse = {}

    ## extract the IPD ratio
    i = 0
    f = open(file_path)
    f.readline()
    for line in f:
        field = line.strip().split(',')
        # print (field[0], field[1], field[2], field[3], field[4], field[5])
        if field[0] != '"NC_000001.11"':
            break
        # print (field[0], field[1], field[2], field[3], field[4], field[5])
        pos = int(field[1]) #- 1
        ipd = float(field[5])
        ipd_sum_control = field[7]
        strand = int(field[2])
        if strand == 0:
            contig_forward_dict[pos] = [ipd]
            ipd_sum_for_control[pos] = ipd_sum_control
            
        # elif strand == 1:
        #     contig_reverse_dict[pos] = [ipd]
        i+= 1
        # if i > 30000:
        #     break
    f.close()
    print (len(contig_forward_dict), len(contig_reverse_dict))

    seq = seq_dict[contig]
    observed_IPD_list = get_IPD_list(len(seq), contig_forward_dict)
    # observed_IPD_reverse_list = get_IPD_list(len(seq), contig_reverse_dict)
    kmer_baseline_dict, save_kmer = prepare_data(seq, observed_IPD_list)

    data = []
    for pos in save_kmer:
        kmer = save_kmer[pos]
        if kmer in kmer_baseline_dict and pos in contig_forward_dict:
            print (pos, kmer, len(kmer_baseline_dict[kmer]), np.mean(kmer_baseline_dict[kmer]), contig_forward_dict[pos][0], ipd_sum_for_control[pos])
            data.append([pos, kmer, len(kmer_baseline_dict[kmer]), np.mean(kmer_baseline_dict[kmer]), contig_forward_dict[pos][0], ipd_sum_for_control[pos]])
    df = pd.DataFrame(data, columns=['pos', 'kmer', 'count', 'estimated', 'real', 'ipd_sum'])
    df.to_csv(test_csv, index=False)

# ref = "/home/shuaiw/methylation/data/borg/hg38/GCF_000001405.40_GRCh38.p14_genomic.fasta"
ref = "/home/shuaiw/methylation/data/borg/hg38/NC_000001.11.fasta"
subread_bam = "/home/shuaiw/methylation/data/borg/human/human_000733.subreads.align.bam"
csv = "/home/shuaiw/methylation/data/borg/human/human_000733_subreads.csv"
test_csv = "/home/shuaiw/methylation/data/borg/human/test_result.csv"

seq_dict = extract_context(ref)
print ("ref loaded")
# read_subread_bam(subread_bam)
extract_ipd_ratio(csv)
print ("csv is loaded")
# print ("bam is loaded", np.median(observed_IPD_list), np.median(observed_IPD_reverse_list))





