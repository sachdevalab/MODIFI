import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics import jaccard_score
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

import sys
import argparse
from adjustText import adjust_text
import numpy as np


def UMAP(df):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T


    # print (matrix)
    ## normlaize the matrix using zscore
    from scipy.stats import zscore
    matrix = zscore(matrix)

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())



    try:
        import umap
        X_embedded = umap.UMAP().fit_transform(matrix, 
                                               n_neighbors=1,
                                               min_dist=0.1,
                                               metric='euclidean')
        # import umap.plot
        # umap.plot.points(X_embedded)
        # ## save the umap result
        # plt.savefig("/tmp/umap.pdf")
    except Exception as e:
        print(f"Failed to create UMAP: {e}")
        return
    
    # from sklearn.mixture import GaussianMixture
    # bic = []
    # max_cluster = 100
    # if len(X_embedded) < max_cluster:
    #     max_cluster = len(X_embedded)
    # for i in range(1, max_cluster+1):
    #     gmm = GaussianMixture(n_components=i)
    #     gmm.fit(X_embedded)
    #     bic.append(gmm.bic(X_embedded))
    # n_clusters = np.argmin(bic) + 1
    # print (n_clusters, "clusters detected in UMAP.")

    # # Fit the GaussianMixture model with the best number of clusters
    # gmm = GaussianMixture(n_components=n_clusters)
    # gmm.fit(X_embedded)
    # cluster_labels = gmm.predict(X_embedded)  # Get cluster labels

    # # Output the cluster result
    # data = []
    # for i in range(n_clusters):
    #     for j in range(len(cluster_labels)):
    #         if cluster_labels[j] == i:
    #             data.append([df.columns[j], i])
    # # Save the cluster result
    # cluster_result = pd.DataFrame(data, columns=['contigs', 'cluster'])
    # cluster_result.to_csv(result_file, index=False)


    clustering = DBSCAN(eps=0.4, min_samples=1).fit(X_embedded)
    # print (clustering.labels_)
    # calculate how many clusters
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected in UMAP.")

    ## output the cluster result, the elements with same cluster label are output together
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(clustering.labels_)):
            if clustering.labels_[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    ## save the cluster result
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(result_file, index=False)

def Hierachy(df, result_file):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    # from scipy.stats import zscore
    # matrix = zscore(matrix)

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())

    my_linkage = linkage(matrix, method='average', metric='euclidean')
    cluster_labels = fcluster(my_linkage, t=1.3, criterion='distance')

    # from sklearn.mixture import GaussianMixture

    ## use bic to determine the number of clusters
    # bic = []
    # max_cluster = 30
    # # if len(matrix) < max_cluster:
    # max_cluster = len(matrix)
    # for i in range(1, max_cluster+1):
    #     gmm = GaussianMixture(n_components=i)
    #     gmm.fit(matrix)
    #     bic.append(gmm.bic(matrix))
    # n_clusters = np.argmin(bic) + 1
    # print (n_clusters, "clusters detected in hierarchical_clustering.")
    # # n_clusters = 100
    # gmm = GaussianMixture(n_components=n_clusters)
    # gmm.fit(matrix)
    # cluster_labels = gmm.predict(matrix)


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
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(result_file, index=False)
    # Plot dendrogram
    # plt.figure(figsize=(20, 10))
    # dendrogram(my_linkage, labels=df.columns, orientation='left', leaf_rotation=0)#, leaf_font_size=8
    # ## show cut line
    # plt.axvline(x=1.5, color='r', linestyle='--')
    # plt.title(f"Dendrogram with Distance Threshold = {1.5}")
    # plt.xlabel("Sample Index")
    # plt.ylabel("Cluster Distance")
    # plt.savefig("tmp/tree.pdf")
    # print ("hierarchical clustering done.")

def JC(df, result_file):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T
    # cutoff = 0.45
    cutoff = 0.45

    ## binarize the matrix
    matrix = (matrix > 0.5).astype(int) 
    my_linkage = linkage(matrix, method='average', metric='jaccard')
    cluster_labels = fcluster(my_linkage, t=cutoff, criterion='distance')
    n_clusters = len(set(cluster_labels))
    print (n_clusters, "clusters detected in JC.")
    # print (cluster_labels)
    data = []
    for i in range(n_clusters+1):
        # print ("cluster", i)
        for j in range(len(cluster_labels)):
            if cluster_labels[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(result_file, index=False)
    # Plot dendrogram
    plt.figure(figsize=(20, 10))
    dendrogram(my_linkage, labels=df.columns, orientation='left', leaf_rotation=0)#, leaf_font_size=8
    ## show cut line
    plt.axvline(x=cutoff, color='r', linestyle='--')
    plt.title(f"Dendrogram with Distance Threshold = {cutoff}")
    plt.xlabel("Sample Index")
    ## save 
    plt.savefig("tmp/tree.j.pdf")
    print ("JC clustering done.")



# profile_file = "/home/shuaiw/methylation/data/borg/bench/zymo2/motif_profile2.csv"
# profile_file = "/home/shuaiw/methylation/data/borg/bench/zymo6_NM200/motif_profile.csv"
profile_file = "/home/shuaiw/methylation/data/borg/bench/zymo6_NM3/motif_profile.csv"
result_file = "tmp/zymo.u.csv"
# profile_file = "/home/shuaiw/methylation/data/borg/bench/all_break/motif_profile.csv"
# result_file = "tmp/all.u.csv"
# profile_file = "/home/shuaiw/methylation/data/borg/all_test_ccs2/motif_profile.csv"
# profile_file = "/home/shuaiw/methylation/data/borg/bench/all_subreads/motif_profile2.csv"
# result_file = "tmp/real.u.csv"

profiles = pd.read_csv(profile_file, index_col=0)
min_frac = 0.5
profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
# Filter columns where any value is greater than 0.5
profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]
print (profiles)
# Hierachy(profiles, result_file)
JC(profiles, result_file)
os.system("python assess_clustering.py")


