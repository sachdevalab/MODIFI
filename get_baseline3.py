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
        seq = str(record.seq)[:50000000]
        ## convert to capital
        raw_seq = seq.upper()
        seq_dict[record.id] = raw_seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict

def prepare_data(seq, control_list, kmer_baseline_dict, up=7, down=3):
    # print(seq[:100])
    # seq, control_list = seq[:1000], control_list[:1000]
    save_kmer = {}  ## {pos:kmer}
    
    # y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        # if seq[i] == 'N':
        #     continue
        kmer = seq[i-up:i+down]
        if 'N' in kmer:
            continue  # Skip kmers containing 'N'
        ipd = control_list[i]
        if ipd == 0:
            continue
        kmer_baseline_dict[kmer].append(ipd)
        save_kmer[i] = kmer    
    return kmer_baseline_dict, save_kmer

class Contig:
    def __init__(self, name):
        self.name = name
        self.forward_dict = {}
        self.reverse_dict = {}
        self.forward_ipd_sum = {}
        self.reverse_ipd_sum = {}

def get_reverse_cmplement(kmer):
    kmer_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    complement_kmer = ''
    for i in range(len(kmer)):
        complement_kmer += kmer_dict[kmer[i]]
    reverse_kmer = complement_kmer[::-1]
    return reverse_kmer

def extract_ipd_ratio_all(file_path):
    ## open it using pandas
    # contig = "NC_000001.11"
    contig_forward_dict_dict = defaultdict(dict)
    contig_reverse_dict_dict = defaultdict(dict)
    ipd_sum_for_control_dict = defaultdict(dict)
    ipd_sum_for_control_reverse_dict = defaultdict(dict)

    df = pd.read_csv(file_path, nrows=10000)
    print ("length of df", len(df))
    for index, row in df.iterrows():
        contig = row['refName']
        contig = contig.replace('"', '')
        if len(contig_forward_dict_dict[contig]) > 50000000:
            continue
        if len(contig_forward_dict_dict[contig]) == 1:
            print (contig)
        if len(contig_forward_dict_dict) > 4:
            break 
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

def count_kmer(contig_forward_dict_dict, seq):

    kmer_baseline_dict = defaultdict(list)
    save_kmer_dict = {}
    for contig in contig_forward_dict_dict:
        contig_forward_dict = contig_forward_dict_dict[contig]

        seq = seq_dict[contig]
        observed_IPD_list = get_IPD_list(len(seq), contig_forward_dict)
        kmer_baseline_dict, save_kmer = prepare_data(seq, observed_IPD_list, kmer_baseline_dict, 8, 4)
        save_kmer_dict[contig] = save_kmer
    print ("kmer is counted")
    # import pickle
    # with open(save_kmer_file, 'wb') as f:
    #     pickle.dump(kmer_baseline_dict, f)
    
    mean_dict, median_dict = {}, {}
    for kmer in kmer_baseline_dict:
        mean_dict[kmer] = np.mean(kmer_baseline_dict[kmer])
        median_dict[kmer] = np.median(kmer_baseline_dict[kmer])
    print ("mean and median is computed")
    return save_kmer_dict, mean_dict, median_dict, kmer_baseline_dict

def align_kmer(contig_forward_dict_dict, ipd_sum_for_control_dict, \
               kmer_baseline_dict, save_kmer_dict, mean_dict, median_dict, test_csv):
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
            # if len(kmer_baseline_dict[kmer]) < 10:
            #     continue
            if kmer in kmer_baseline_dict and pos in contig_forward_dict:
                # print (pos, kmer, len(kmer_baseline_dict[kmer]), np.mean(kmer_baseline_dict[kmer]), contig_forward_dict[pos][0], ipd_sum_for_control[pos])
                data.append([contig, pos, kmer, len(kmer_baseline_dict[kmer]), \
                mean_dict[kmer], median_dict[kmer], \
                contig_forward_dict[pos][0], ipd_sum_for_control[pos]])

    df = pd.DataFrame(data, columns=['contig', 'pos', 'kmer', 'count', 'estimated', 'median', 'real', 'ipd_sum'])
    df.to_csv(test_csv, index=False)


    

if __name__ == "__main__":

    ref = "/home/shuaiw/methylation/data/borg/hg38/GCF_000001405.40_GRCh38.p14_genomic.fasta"
    # ref = "/home/shuaiw/methylation/data/borg/hg38/NC_000001.11.fasta"
    subread_bam = "/home/shuaiw/methylation/data/borg/human/human_000733.subreads.align.bam"
    csv = "/home/shuaiw/methylation/data/borg/human/human_000733_subreads.csv"
    test_csv = "/home/shuaiw/methylation/data/borg/human/test_result5.csv"
    save_kmer_file = "/home/shuaiw/methylation/data/borg/human/kmer_baseline_dict.pkl"

    seq_dict = extract_context(ref)
    print ("ref loaded")
    # read_subread_bam(subread_bam)
    # extract_ipd_ratio_all(csv)
    # print ("csv is loaded")
    # print ("bam is loaded", np.median(observed_IPD_list), np.median(observed_IPD_reverse_list))

    contig_forward_dict_dict, contig_reverse_dict_dict, ipd_sum_for_control_dict, ipd_sum_for_control_reverse_dict = extract_ipd_ratio_all(csv)
    print ("ipd is loaded")
    save_kmer_dict, mean_dict, median_dict, kmer_baseline_dict = count_kmer(contig_forward_dict_dict, seq_dict)
    print ("kmer is counted")
    align_kmer(contig_forward_dict_dict, ipd_sum_for_control_dict, kmer_baseline_dict, save_kmer_dict, mean_dict, median_dict, test_csv)
    print ("kmer is aligned")
