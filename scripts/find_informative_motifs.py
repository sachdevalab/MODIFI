
import pandas as pd
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster, linkage

def hierarchical_clustering(df, tree_fig, cluster_fig, cutoff=1.6):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T
    print (matrix)

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    # mask = matrix == 0
    # matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())

    my_linkage = linkage(matrix, method='average', metric='euclidean')
    cluster_labels = fcluster(my_linkage, t=cutoff, criterion='distance')

    ## calculate how many clusters
    n_clusters = len(set(cluster_labels))
    print (n_clusters, "clusters detected hierarchical_clustering.")
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(cluster_labels)):
            if cluster_labels[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
                print (df.columns[j], i)
    cluster_result = pd.DataFrame(data, columns = ['motifs', 'cluster'])
    cluster_result.to_csv(cluster_fig.replace(".pdf", ".h.csv"), index=False)
    # Plot dendrogram
    try:
        plt.figure(figsize=(20, 15))
        dendrogram(my_linkage, labels=df.columns, orientation='left', leaf_rotation=0)
        # plt.axhline(y=cutoff, color='r', linestyle='--')  # Show the cutoff threshold
        plt.axvline(x=cutoff, color='r', linestyle='--')
        plt.title(f"Dendrogram with Distance Threshold = {cutoff}")
        plt.xlabel("Sample Index")
        plt.ylabel("Cluster Distance")
        plt.savefig(tree_fig)
    except Exception as e:
        print(f"Failed to create dendrogram: {e}")
    print ("hierarchical clustering done.")

def find_informative_motifs(motif_profile):
    tree_fig = "tmp/tree.png"
    cluster_fig = "tmp/cluster.png"
    # Read the motif profile file
    df = pd.read_csv(motif_profile)
    ## transpose the dataframe
    ## set motifString as index
    df = df.set_index('motifString')
    df = df.T
    print (df)
    ## cluster the motifs which 
    # hierarchical_clustering(df, tree_fig, cluster_fig, cutoff=45)

def count_motif(motif_profile):
    # Read the motif profile file
    df = pd.read_csv(motif_profile)
    ## count the average of each motif in all contigs
    df = pd.read_csv(motif_profile)
    df = df.set_index('motifString')
    df = df.T

    ## calculate the average of each motif
    motif_avg = df.mean(axis=0)
    motif_avg = pd.DataFrame(motif_avg)
    motif_avg = motif_avg.reset_index()
    motif_avg.columns = ['motifString', 'average']
    motif_avg = motif_avg.sort_values(by='average', ascending=False)
    motif_avg.to_csv(motif_profile.replace(".csv", ".motif_avg.csv"), index=False)

    ## plot the average of each motif
    plt.figure(figsize=(20, 15))
    ## using histogram to plot the average of each motif using seaborn
    import seaborn as sns
    sns.set(style="whitegrid")
    ## plot distribution of the average of each motif
    sns.histplot(motif_avg['average'], bins=100, kde=True)
    plt.xlabel("Average of Each Motif")
    plt.ylabel("Frequency")
    plt.savefig(motif_profile.replace(".csv", ".motif_avg.png"))
    plt.close()





motif_profile = "/home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/motif_profile.csv"
# find_informative_motifs(motif_profile)
count_motif(motif_profile)

