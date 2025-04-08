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

report = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break/checkM2/quality_report.tsv"
# report = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break/checkM_metabat/quality_report.tsv"
read_report(report)