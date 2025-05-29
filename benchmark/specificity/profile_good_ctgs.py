import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE


def get_best_ctg(min_len = 1000000):
    """
    Get the best contig based on length from a fasta file.
    """
    best_ctgs = []
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            if ctg[-1] == "C" and length >= min_len:
                best_ctgs.append(ctg)
    print (f"Total {len(best_ctgs)} contigs with length >= {min_len} found.")
    return best_ctgs

def plot_heatmap(profile_file, min_frac=0.4):
    profiles = pd.read_csv(profile_file, index_col=0)
    ## only keep the columns that are in the best contigs
    best_ctgs = get_best_ctg()
    profiles = profiles.loc[:, profiles.columns.isin(best_ctgs)]
    print ("original shape", profiles.shape)                                    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    print ("filtered shape 1", profiles.shape)
    # Filter columns where any value is greater than 0.5
    profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]

    print ("filtered shape", profiles.shape)

    df = profiles.T

    sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(20, 15))
    plt.savefig("../../tmp/results/heatmap.png", dpi=300, bbox_inches='tight')
    plt.clf()


    # matrix = df.to_numpy()
    # mask = matrix == 0
    # matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())

    # X_embedded = TSNE(n_components=2).fit_transform(matrix)
    
    # ## define fig size
    # plt.figure(figsize=(10, 10))
    # ## plot the cluster result using seaborn
    # scatter_plot = sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], palette="viridis")
    # plt.savefig("../../tmp/results/tsne_plot.png", dpi=300, bbox_inches='tight')
    # plt.clf()


fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"
# get_best_ctg()
plot_heatmap(profile_file)