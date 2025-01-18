import pysam
import numpy as np
from collections import defaultdict
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import random
import pickle
import os
from scipy.stats import pearsonr


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
        print(record.id)
        # if record.id != "NC_000001.11":
        #     continue
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)#[:50000000]
        ## convert to capital
        raw_seq = seq.upper()
        seq_dict[record.id] = raw_seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict

def prepare_data(seq, control_list, kmer_baseline_dict, kmer_num_dict, up=7, down=3):
    save_kmer = {}  ## {pos:kmer}
    
    # y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        kmer = seq[i-up:i+down]
        if 'N' in kmer:
            continue  # Skip kmers containing 'N'
        ipd = control_list[i]
        if ipd == 0:
            continue
        # kmer_baseline_dict[kmer].append(ipd)
        kmer_baseline_dict[kmer] += ipd
        kmer_num_dict[kmer] += 1
        save_kmer[i] = kmer    
    return kmer_baseline_dict, kmer_num_dict, save_kmer

class Contig:
    def __init__(self, name):
        self.name = name
        self.forward_dict = {}
        self.reverse_dict = {}
        self.forward_ipd_sum = {}
        self.reverse_ipd_sum = {}

def get_reverse_cmplement(kmer):
    kmer_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    complement_kmer = ''
    for i in range(len(kmer)):
        complement_kmer += kmer_dict[kmer[i]]
    # reverse_kmer = complement_kmer[::-1]
    # return reverse_kmer
    return complement_kmer

def extract_ipd_ratio_all(file_path):
    contig_forward_dict_dict = defaultdict(dict)
    contig_reverse_dict_dict = defaultdict(dict)
    ipd_sum_for_control_dict = defaultdict(dict)
    ipd_sum_for_control_reverse_dict = defaultdict(dict)

    # df = pd.read_csv(file_path, nrows=1000)
    df = pd.read_csv(file_path)
    print ("length of df", len(df))
    for index, row in df.iterrows():
        contig = row['refName']
        contig = contig.replace('"', '')
        # if len(contig_forward_dict_dict[contig]) > 50000000:
        #     continue
        # if len(contig_forward_dict_dict[contig]) == 1:
        #     print (contig)
        # if len(contig_forward_dict_dict) > 4:
        #     break 
        pos = int(row['tpl']) 
        ipd = float(row['tMean'])
        ipd_sum_control = row['modelPrediction']

        if int(row['strand']) == 1:
            contig_forward_dict_dict[contig][pos] = [ipd]
            ipd_sum_for_control_dict[contig][pos] = ipd_sum_control
        else:
            contig_reverse_dict_dict[contig][pos] = [ipd]
            ipd_sum_for_control_reverse_dict[contig][pos] = ipd_sum_control

    print ("ipd is loaded", len(contig_forward_dict_dict))
    return contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict

class IPD_sum:
    def __init__(self,forward_ipd, reverse_ipd, forward_ipd_sum, reverse_ipd_sum):
        self.forward_ipd = forward_ipd
        self.reverse_ipd = reverse_ipd
        self.forward_ipd_sum = forward_ipd_sum
        self.reverse_ipd_sum = reverse_ipd_sum

def process_chunk(chunk, contig_forward_dict_dict, contig_reverse_dict_dict, \
                  ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict):
    print("Processing chunk of size", len(chunk))
    for index, row in chunk.iterrows():
        contig = row['refName']
        contig = contig.replace('"', '')
        # Add your processing logic here
        pos = int(row['tpl']) 
        ipd = float(row['tMean'])
        ipd_sum_control = row['modelPrediction']

        if int(row['strand']) == 1:
            contig_forward_dict_dict[contig][pos] = [ipd]
            ipd_sum_for_control_dict[contig][pos] = ipd_sum_control
        else:
            contig_reverse_dict_dict[contig][pos] = [ipd]
            ipd_sum_for_control_reverse_dict[contig][pos] = ipd_sum_control
    return contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict

def extract_ipd_ratio_all_chunk(file_path, chunksize=10000):
    contig_forward_dict_dict = defaultdict(dict)
    contig_reverse_dict_dict = defaultdict(dict)
    ipd_sum_for_control_dict = defaultdict(dict)
    ipd_sum_for_control_reverse_dict = defaultdict(dict)

    chunk_iter = pd.read_csv(file_path, chunksize=chunksize,nrows=100000)
    for chunk in chunk_iter:
        contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, \
            ipd_sum_for_control_reverse_dict = process_chunk(chunk, contig_forward_dict_dict,\
                                                              contig_reverse_dict_dict, ipd_sum_for_control_dict, \
                                                                ipd_sum_for_control_reverse_dict)
    print ("ipd is loaded", len(contig_forward_dict_dict))
    return contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict

def count_kmer(contig_forward_dict_dict, seq, strand = 1):

    # kmer_baseline_dict = defaultdict(list)
    kmer_baseline_dict = defaultdict(int)   ## sum
    kmer_num_dict = defaultdict(int)  # count
    save_kmer_dict = {}
    for contig in contig_forward_dict_dict:
        contig_forward_dict = contig_forward_dict_dict[contig]

        seq = seq_dict[contig]
        if strand == 0:
            seq = get_reverse_cmplement(seq)
        observed_IPD_list = get_IPD_list(len(seq), contig_forward_dict)
        kmer_baseline_dict, kmer_num_dict, save_kmer = prepare_data(seq, observed_IPD_list, kmer_baseline_dict, kmer_num_dict, 8, 4)
        save_kmer_dict[contig] = save_kmer
    print ("kmer is counted")
    
    mean_dict, median_dict = {}, {}
    for kmer in kmer_baseline_dict:
        mean_dict[kmer] = round(kmer_baseline_dict[kmer]/kmer_num_dict[kmer], 3) ## calculate the mean of each kmer
        median_dict[kmer] = 0
    print ("mean and median is computed", len(mean_dict), 'kmers')
    return save_kmer_dict, mean_dict, median_dict, kmer_baseline_dict, kmer_num_dict

def align_kmer(contig_forward_dict_dict, ipd_sum_for_control_dict, \
               kmer_baseline_dict, kmer_num_dict, save_kmer_dict, mean_dict, median_dict, test_csv, strand = 1):
    data = []
    for contig in contig_forward_dict_dict:
        if contig not in ipd_sum_for_control_dict:
            continue
        contig_forward_dict = contig_forward_dict_dict[contig]
        ipd_sum_for_control = ipd_sum_for_control_dict[contig]
        save_kmer = save_kmer_dict[contig] 
        print (contig, len(contig_forward_dict), len(ipd_sum_for_control), len(save_kmer))
        for pos in save_kmer:
            kmer = save_kmer[pos]
            if kmer_num_dict[kmer] < 10:
                continue
            if kmer in kmer_baseline_dict and pos in contig_forward_dict:
                # print (pos, kmer, len(kmer_baseline_dict[kmer]), np.mean(kmer_baseline_dict[kmer]), contig_forward_dict[pos][0], ipd_sum_for_control[pos])
                data.append([contig, pos, strand, kmer, kmer_num_dict[kmer], \
                mean_dict[kmer], median_dict[kmer], \
                contig_forward_dict[pos][0], ipd_sum_for_control[pos]])

    df = pd.DataFrame(data, columns=['contig', 'pos', 'strand', 'kmer', 'count', 'estimated', 'median', 'real', 'ipd_sum'])
    df.to_csv(test_csv, index=False)


    

if __name__ == "__main__":   
    # ref = "/home/shuaiw/methylation/data/borg/hg38/NC_000001.11.fasta"
    # subread_bam = "/home/shuaiw/methylation/data/borg/human/human_000733.subreads.align.bam"

    ref = "/home/shuaiw/Methy/borg/contigs/break_contigs.fasta"
    csv = "/home/shuaiw/Methy/borg/break_contigs/break_contigs.csv"
    test_csv = "/home/shuaiw/methylation/data/borg/human/test_result6.csv"

    save_kmer_file = "/home/shuaiw/methylation/data/borg/human/kmer_baseline_dict.pkl"

    seq_dict = extract_context(ref)
    print ("ref loaded")

    contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict = extract_ipd_ratio_all_chunk(csv)
    print ("ipd is loaded")
    # contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict = extract_ipd_ratio_all(csv)
    # print ("ipd is loaded")

    save_kmer_dict, mean_dict, median_dict, kmer_baseline_dict, kmer_num_dict = count_kmer(contig_forward_dict_dict, seq_dict, 1)
    print ("kmer is counted")
    align_kmer(contig_forward_dict_dict, ipd_sum_for_control_dict, kmer_baseline_dict, kmer_num_dict, save_kmer_dict, mean_dict, median_dict, test_csv, 1)
    print ("kmer is aligned")

    # save_kmer_dict, mean_dict, median_dict, kmer_baseline_dict = count_kmer(contig_reverse_dict_dict, seq_dict, 0)
    # print ("kmer is counted")
    # align_kmer(contig_reverse_dict_dict, ipd_sum_for_control_reverse_dict, kmer_baseline_dict, save_kmer_dict, mean_dict, median_dict, test_csv, 0)
    # print ("kmer is aligned")
