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
    score_list = []
    df = pd.read_csv(file_path)
    # print ("csv loaded")
    ### only retain the strand == 0
    df = df[df['strand'] == 1]
    ### convert the modelPrediction column to list
    # control_list = df['modelPrediction'].tolist()
    control_list = df['tMean'].tolist()
    score_list = df['score'].tolist()
    return control_list, score_list

def extract_context(fasta):
    ## load the fasta using biopython
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        print(record.id)
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        seq = seq.upper()
        seq = seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
        return seq

def prepare_data(seq, control_list, score_list, up=8, down=4):
    # print(seq[:100])
    X = []
    y = score_list[up:len(seq) - down]
    ### convert the element in control_list to string
    control_list = [str(i) for i in control_list]
    for i in range(up, len(seq) - down):
        X.append(seq[i-up:i+down] + ' ' + ' '.join(control_list[i-up:i+down]))
    new_y = []
    for i in range(len(y)):
        if y[i] > 10:
            new_y.append(1)
        else:
            new_y.append(0)
    y = new_y
    return X, y

def get_csv(X, y, csv):
    with open(csv, 'w') as f:
        for i in range(len(X)):
            f.write(str(y[i]) + "\t" + str(X[i]) + "\n")

def model(X, y, seed):
    # Assuming X and y are defined and properly shaped
    X = np.array(X).reshape(-1, 1)
    y = np.array(y)
    # 
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

def run_single_contig(csv, fasta):
    seq = extract_context(fasta)
    control_list,score_list = extract_ipd_ratio(csv)
    # print ("length of control list", len(control_list), len(seq))
    if len(control_list) != len(seq):
        print ("different length of control list", len(control_list), len(seq), csv)
        return [], []
    X, y = prepare_data(seq, control_list,score_list)
    return X, y

def run_multiple_contigs(dir):
    all_X, all_y = [], []
    i = 0
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            csv = dir + file
            fasta = dir + file.replace(".csv", ".fasta")
            if not os.path.exists(fasta):
                continue
            X, y = run_single_contig(csv, fasta)
            if len(X) == 0:
                continue
            i += 1
            print ("length of X, y", len(X), len(y), i)
            all_X += X
            all_y += y
            if i > 10:
                break
    return all_X, all_y

if __name__ == "__main__":
    seed = 8
    # bam_file = "borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam"
    # bam_file = "all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"

    csv = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.csv"
    fasta = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.fasta"
    dir = "borg/split_bam_dir/"
    
    # X, y = run_single_contig(csv, fasta)
    X, y = run_multiple_contigs(dir)
    csv = 'borg/train_data.csv'
    get_csv(X, y, csv)

    # model(X, y, seed)



