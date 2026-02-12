from sklearn.datasets import make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

def extract_ipd_ratio(file_path):
    ## open it using pandas
    IPD_ratio_list, reverse_IPD_ratio_list = [], []
    IPD_list, reverse_IPD_list = [], []
    control_list, reverse_control_list = [], []
    df = pd.read_csv(file_path)
    # print ("csv loaded")
    ### only retain the strand == 0
    df = df[df['strand'] == 1]
    ### convert the modelPrediction column to list
    # control_list = df['modelPrediction'].tolist()
    control_list = df['tMean'].tolist()
    return control_list

def extract_context(fasta):
    ## load the fasta using biopython
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        print(record.id)
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        seq = seq.upper()
        # seq = seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
        return seq

def prepare_data(seq, control_list, up=8, down=4):
    # print(seq[:100])
    X = []
    y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        X.append(seq[i-up:i+down])
    return X, y

def count_kmer(X, y):
    # X, y = X[:1000], y[:1000]
    kmer_dict = {}
    for i in range(len(X)):
        kmer = X[i]
        if kmer not in kmer_dict:
            kmer_dict[kmer] = []
        kmer_dict[kmer].append(y[i])
    ## sort the kmer_dict by the length of the value
    new_kmer_dict = {}
    for kmer in kmer_dict:
        if len(kmer_dict[kmer]) > 100:
            new_kmer_dict[kmer] = kmer_dict[kmer]
    kmer_dict = new_kmer_dict
    kmer_dict = dict(sorted(kmer_dict.items(), key=lambda item: max(item[1]), reverse=True))
    #kmer_dict = dict(sorted(kmer_dict.items(), key=lambda item: len(item[1]), reverse=True))
    plot_kmer(kmer_dict)

    # i = 0
    # data = []
    # for kmer in kmer_dict:
    #     # print (kmer, len(kmer_dict[kmer]), kmer_dict[kmer], np.mean(kmer_dict[kmer]), np.std(kmer_dict[kmer]))
    #     data.append([kmer, len(kmer_dict[kmer]), np.mean(kmer_dict[kmer]), np.std(kmer_dict[kmer]), min(kmer_dict[kmer]), max(kmer_dict[kmer])])
    #     i += 1
    #     if i > 100:
    #         break
    # df = pd.DataFrame(data, columns=['kmer', 'count', 'mean', 'std', 'min', 'max'])
    # df.to_csv("/home/shuaiw/Methy/borg/customized/kmer_ipd_distribution_k4.csv", index=False)
    # # return kmer_dict

def plot_kmer(kmer_dict):
    data = []
    i = 0
    ## sort the kmer_dict by the count of the value
    kmer_dict = dict(sorted(kmer_dict.items(), key=lambda item: len(item[1]), reverse=True))
    for kmer in kmer_dict:
        for y in kmer_dict[kmer]:
            data.append([kmer, y[0], y[1]])
        i += 1
        if i > 50000:
            break
    print ("IPD num", len(data))
    df = pd.DataFrame(data, columns=['kmer', 'ipd', 'contig'])
    df.to_csv("/home/shuaiw/borg/paper/workflow_fig/kmer_ipd_distribution.csv", index=False)
    ### plot a histogram distribution of the ipd for each kmer, plot each kmer in a subfigure 

    ## define figure size
    # plt.figure(figsize=(10, 50))
    # # Create a FacetGrid for subplots
    # g = sns.FacetGrid(df, col="kmer", col_wrap=1, sharex=True, sharey=False, aspect=3 )
    # g.map(sns.histplot, "ipd", bins=30)
    
    # # Adjust the layout
    # g.tight_layout()

    
    # # Save the plot
    # plt.savefig("/home/shuaiw/mGlu/tmp/figures/framework/kmer_ipd_histograms.pdf")
    
def run_single_contig(csv, fasta):
    seq = extract_context(fasta)
    control_list = extract_ipd_ratio(csv)
    # print ("length of control list", len(control_list), len(seq))
    if len(control_list) != len(seq):
        print ("different length of control list", len(control_list), len(seq), csv)
        return [], []
    X, y = prepare_data(seq, control_list, 7, 3)
    return X, y

def run_multiple_contigs(dir):
    all_X, all_y = [], []
    i = 0
    fasta_dir = os.path.join(dir, "contigs/")
    ipd_dir = os.path.join(dir, "ipd/")
    for file in os.listdir(ipd_dir):
        if file.endswith(".csv"):
            csv = os.path.join(ipd_dir, file)
            contig_name = file.replace(".ipd1.csv", "")
            fasta = os.path.join(fasta_dir, file.replace(".ipd1.csv", ".fa"))
            if not os.path.exists(fasta):
                continue
            X, y = run_single_contig(csv, fasta)
            if len(X) == 0:
                continue
            i += 1
            print ("length of X, y", len(X), len(y), i)
            all_X += X
            for value in y:
                all_y.append((value, contig_name))
            # all_y += y
            if i > 50:
                break
    return all_X, all_y

if __name__ == "__main__":
    seed = 8
    # bam_file = "borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam"
    # bam_file = "all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"

    # csv = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.csv"
    # fasta = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.fasta"
    # dir = "borg/split_bam_dir/"
    dir = "/home/shuaiw/borg/paper/run2/96plex/96plex_methylation4/"
    
    # X, y = run_single_contig(csv, fasta)
    X, y = run_multiple_contigs(dir)
    count_kmer(X, y)
    

    # model(X, y, seed)



