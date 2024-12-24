from sklearn.datasets import make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np

def extract_ipd_ratio(file_path):
    ## open it using pandas
    IPD_ratio_list, reverse_IPD_ratio_list = [], []
    IPD_list, reverse_IPD_list = [], []
    control_list, reverse_control_list = [], []
    df = pd.read_csv(file_path)
    print ("loaded")
    ### only retain the strand == 0
    df = df[df['strand'] == 1]
    ### convert the modelPrediction column to list
    # control_list = df['modelPrediction'].tolist()
    control_list = df['ipdRatio'].tolist()
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
        seq = seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
        return seq

def prepare_data(seq, control_list):
    print(seq[:100])
    X = []
    y = control_list[8:]
    for i in range(8, len(seq)):
        X.append(seq[i-8:i+4])
    return X, y

gff = "split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.reprocess.gff"
csv = "split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.csv"
fasta = "split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.fasta"
# bam_file = "all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
bam_file = "customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam"

# X, y = make_regression(random_state=0)
# print (X[:2], "y", y[:2])

seed = 8

seq = extract_context(fasta)
control_list = extract_ipd_ratio(csv)
print ("length of control list", len(control_list), len(seq))
X, y = prepare_data(seq, control_list)
print (len(X), len(y))
# Assuming X and y are defined and properly shaped
X = np.array(X).reshape(-1, 1)
y = np.array(y)
# 
X_train, X_test, y_train, y_test = train_test_split(
    X, y, random_state=seed, test_size=0.1)
# reg = HistGradientBoostingRegressor(random_state=0, n_jobs=-1)

reg = HistGradientBoostingRegressor(
    learning_rate=0.15,
    max_iter=60000,
    max_depth=12,
    random_state=seed,
)
reg.fit(X_train, y_train)
# reg.predict(X_test)
result = reg.score(X_test, y_test)
print (result)
print (reg.predict(X_test[1:5]))
print (y_test[1:5])


