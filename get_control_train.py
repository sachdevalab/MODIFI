import control_model

from sklearn.datasets import make_regression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import os




def run_multiple_contigs(dir, up=8, down=4):
    all_X, all_y = [], []
    i = 0
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            csv = dir + file
            fasta = dir + file.replace(".csv", ".fasta")
            if not os.path.exists(fasta):
                continue
            X, y = control_model.run_single_contig(csv, fasta, up, down)
            if len(X) == 0:
                continue
            i += 1
            print ("length of X, y", len(X), len(y), i)
            all_X += X
            all_y += y
            if i > 10:
                break
    return all_X, all_y

def out_data(X, y, file):
    df = pd.DataFrame(X, columns=['x'+str(i) for i in range(len(X[0]))])
    df['y'] = y
    print (len(df))
    ## randomly subsample
    df = df.sample(frac=0.1).reset_index(drop=True)
    print ("after subsample", len(df))
    df.to_csv(file, index=False)

if __name__ == "__main__":
    dir = "borg/split_bam_dir/"
    control_train_file = dir + '/control_train_data.csv'
    X, y = run_multiple_contigs(dir)
    out_data(X, y, control_train_file)

