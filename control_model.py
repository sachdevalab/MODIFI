from sklearn.datasets import make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import os


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
    cal_correlation(df)
    return control_list

def cal_correlation(df):
    from scipy.stats import pearsonr
    if len(df['modelPrediction']) == len(df['tMean']):
        corr = pearsonr(df['modelPrediction'], df['tMean'])
        print ("correlation", corr)

def extract_context(fasta):
    ## load the fasta using biopython
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        print(record.id)
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        raw_seq = seq.upper()
        seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
        return seq, raw_seq

def prepare_data(seq, control_list, raw_seq, up=8, down=4):
    # print(seq[:100])
    # seq, control_list = seq[:1000], control_list[:1000]
    X, S = [], []
    y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        kmer = []
        for z in range(i-up, i+down):
            # print (z, i-up, i+down)
            kmer.append(int(seq[z]))
        S.append(raw_seq[i-up:i+down])
        # print (kmer)
        X.append(kmer)
        # break
        # X.append(seq[i-up:i+down])
    return X, y, S

def model(X, y, seed):
    # Assuming X and y are defined and properly shaped
    # X = np.array(X).reshape(-1, 1)
    X = np.array(X)
    y = np.array(y)
    print (len(X), len(y))
    for i in range(100):
        print (X[i], y[i])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=seed, test_size=0.1)

    # reg = HistGradientBoostingRegressor(
    #     learning_rate=0.15,
    #     max_iter=60000,
    #     max_depth=12,
    #     random_state=seed,
    # )

    reg = RandomForestRegressor(
            n_estimators=100,
            max_depth=12,
            random_state=seed,
    )

    # reg = SVR(    ### very slow
    #         kernel='rbf',
    #         C=1.0,
    #         epsilon=0.1
    #     )

    reg.fit(X_train, y_train)
    # reg.predict(X_test)
    result = reg.score(X_test, y_test)
    print (result)
    print (reg.predict(X_test[1:5]))
    print (y_test[1:5])

def run_single_contig(csv, fasta, up=8, down=4):
    seq, raw_seq = extract_context(fasta)
    control_list = extract_ipd_ratio(csv)
    # print ("length of control list", len(control_list), len(seq))
    if len(control_list) != len(seq):
        print ("different length of control list", len(control_list), len(seq), csv)
        return [], [], []
    X, y, S = prepare_data(seq, control_list,raw_seq, up, down)
    return X, y, S

def run_multiple_contigs(dir, up=8, down=4):
    all_X, all_y, all_S = [], [], []
    i = 0
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            csv = dir + file
            fasta = dir + file.replace(".csv", ".fasta")
            if not os.path.exists(fasta):
                continue
            X, y, S = run_single_contig(csv, fasta, up, down)
            if len(X) == 0:
                continue
            i += 1
            print ("length of X, y", len(X), len(y), i)
            all_X += X
            all_y += y
            all_S += S
            if i > 1:
                break
    return all_X, all_y, all_S

if __name__ == "__main__":
    seed = 8
    # bam_file = "borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam"
    # bam_file = "all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"

    csv = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.csv"
    fasta = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.fasta"
    dir = "borg/split_bam_dir/"
    
    # X, y = run_single_contig(csv, fasta)
    X, y = run_multiple_contigs(dir)

    model(X, y, seed)



