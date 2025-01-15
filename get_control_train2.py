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
    ## build a dataframe
    df1 = pd.DataFrame()
    i = 0
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            csv = dir + file
            fasta = dir + file.replace(".csv", ".fasta")
            if not os.path.exists(fasta):
                continue
            X, y, S = control_model.run_single_contig(csv, fasta, up, down)
            if len(X) == 0:
                continue
            ### add X, y, S to the dataframe
            df1 = pd.DataFrame(X, columns=['x'+str(i) for i in range(len(X[0]))])
            df1['sequence'] = S
            df1['label'] = y
            ## randomly subsample df1
            df1 = df1.sample(frac=0.01).reset_index(drop=True)
            ## append to the main dataframe
            if i == 0:
                df = df1
            else:
                df = pd.concat([df, df1], ignore_index=True)

            i += 1
            print ("length of X, y", len(X), len(y), i)

            if i > 1:
                break
    print ("after subsample", len(df))
    df.to_csv(control_train_file, index=False)
    

if __name__ == "__main__":
    dir = "borg/split_bam_dir/"
    up, down = 250, 250
    control_train_file = dir + f'/control_train_data_{up}_{down}.csv'
    run_multiple_contigs(dir, up, down)


