"""
Given the checkM's result, plot dot plot where x is contamination and y is completeness
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple
from Bio import SeqIO




def read_report(report):
    df = pd.read_csv(report, sep="\t")
    ## plot contamination and completeness
    df = df[["Contamination", "Completeness"]]
    df.columns = ["contamination", "completeness"]
    df = df.dropna()
    ## plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df, x="contamination", y="completeness")
    plt.xlabel("Contamination")
    plt.ylabel("Completeness")
    plt.title("Contamination vs Completeness")
    ## save figure
    plt.savefig(report.replace(".tsv", ".png"))
    plt.close()

    print (len(df))
    high_quality_bins = len(df[(df["completeness"] > 99) & (df["contamination"] < 1)])
    print ("Number of high quality bins: ", high_quality_bins,high_quality_bins/len(df)) 
    ## count the number of bins with completeness < 30
    df["completeness"] = df["completeness"].astype(float)
    df["contamination"] = df["contamination"].astype(float)
    df1 = df[df["completeness"] < 30]
    print("Number of bins with completeness < 30: ", len(df1))
    ### count the number of bins with contamination > 20
    df2 = df[df["contamination"] > 20]
    print("Number of bins with contamination > 20: ", len(df2))

# report = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break/checkM2/quality_report.tsv"
report = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break/checkM_metabat/quality_report.tsv"
read_report(report)