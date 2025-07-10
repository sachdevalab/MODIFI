


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity

from Bio.Seq import Seq
import sys


def get_circular_ctgs(fai, min_len = 1000000):
    """
    Get the best contig based on length from a fasta file.
    """
    circular_ctgs = []
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            if ctg[-1] == "C" and length >= min_len:
                circular_ctgs.append(ctg)
    print (f"Total {len(circular_ctgs)} circular contigs with length >= {min_len} found.")
    return circular_ctgs

def get_high_quality_ctgs(checkM_report):
    """
    Get the high quality contigs from a CheckM report.
    """
    high_quality_ctgs = []
    df = pd.read_csv(checkM_report, sep="\t", header=0)
    df = df[df['Completeness'] >= 50]
    df = df[df['Contamination'] <= 5]
    for ctg in df['Name']:
        high_quality_ctgs.append(ctg)

    print (f"Total {len(high_quality_ctgs)} high quality contigs found.")
    return high_quality_ctgs

def get_best_ctg(fai, checkM_report, high_quality_file, min_len = 1000000):
    high_quality_ctgs = get_high_quality_ctgs(checkM_report)
    circular_ctgs = get_circular_ctgs(fai, min_len)
    selected_ctgs = list(set(high_quality_ctgs) | set(circular_ctgs))
    print (f"Total {len(selected_ctgs)} contigs selected.")
    # return selected_ctgs
    f = open(high_quality_file, "w")
    for ctgs in selected_ctgs:
        f.write(f"{ctgs}\t{ctgs}\n")
    f.close()


if __name__ == "__main__":
    # Example usage
    # workdir = "/home/shuaiw/borg/assembly/96plex/96plex_p5_2"
    # prefix = "96plex_p5"

    workdir = sys.argv[1]
    prefix = sys.argv[2]

    fai = f"{workdir}/{prefix}.final.fa.fai"
    checkM_report = f"{workdir}/checkM2/quality_report.tsv"
    high_quality_file = f"{workdir}/checkM2/high_quality_ctgs.tab"
    selected_ctgs = get_best_ctg(fai, checkM_report, high_quality_file, min_len=1000000)
    # print(selected_ctgs)




