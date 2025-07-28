import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from typing import Dict, List
import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import sys
import argparse
from adjustText import adjust_text
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import umap



def load_cluster_assignments(cluster_csv: str) -> Dict[int, List[str]]:
    """Load contig-cluster assignments from a CSV file."""
    cluster_df = pd.read_csv(cluster_csv)
    cluster_map = defaultdict(list)
    for _, row in cluster_df.iterrows():
        contig = str(row["contigs"])
        cluster = int(row["cluster"])
        cluster_map[cluster].append(contig)
    return cluster_map


def load_cluster_assignments_df(cluster_df: pd.DataFrame) -> pd.DataFrame:
    """Load contig-cluster assignments from a DataFrame."""
    cluster_map = defaultdict(list)
    for _, row in cluster_df.iterrows():
        contig = str(row["contigs"])
        cluster = int(row["cluster"])
        cluster_map[cluster].append(contig)
    return cluster_map

def load_fasta_sequences(fasta_path: str) -> Dict[str, SeqIO.SeqRecord]:
    """Load all contig sequences into a dictionary."""
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))


def write_bins_to_fasta(cluster_map: Dict[int, List[str]],
                        sequences: Dict[str, SeqIO.SeqRecord],
                        output_prefix: str = "bin", max_len: int = 20000000) -> None:
    """Write each cluster/bin as a separate FASTA file."""
    for cluster_id, contigs in cluster_map.items():
        output_file = f"{output_prefix}_{cluster_id}.fasta"
        records = [sequences[ctg] for ctg in contigs if ctg in sequences]
        # Filter out bins that are too long
        total_length = sum(len(record.seq) for record in records)
        if total_length > max_len:
            print(f"Skipping bin {cluster_id} with total length {total_length} > {max_len}")
            continue
        SeqIO.write(records, output_file, "fasta")
        print(f"✔️  Wrote {len(records)} sequences to {output_file}")


def bin_contigs_to_fastas(cluster_csv: str,
                          fasta_path: str,
                          output_prefix: str = "bin") -> None:
    """Main function to process input and create bin FASTA files."""
    ## check if the input files exist
    if not os.path.exists(cluster_csv):
        print (f"Cluster CSV file not found: {cluster_csv}")
        return
    cluster_map = load_cluster_assignments(cluster_csv)
    sequences = load_fasta_sequences(fasta_path)
    write_bins_to_fasta(cluster_map, sequences, output_prefix)
    print("✅ All bins written to FASTA files.")

def bin_contigs_to_fastas_df(cluster_df: pd.DataFrame,
                          fasta_path: str,
                          output_prefix: str = "bin") -> None:
    cluster_map = load_cluster_assignments_df(cluster_df)
    sequences = load_fasta_sequences(fasta_path)
    write_bins_to_fasta(cluster_map, sequences, output_prefix)
    print("✅ All bins written to FASTA files.")


def TSE(df):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())
    # print (matrix)

    try:
        ## reduce dimention using t-SNE
        X_embedded = TSNE(n_components=2).fit_transform(matrix)
    except Exception as e:
        print(f"Failed to create t-SNE: {e}")
        return
    
    clustering = DBSCAN(eps=0.2, min_samples=1).fit(X_embedded)
    # print (clustering.labels_)
    # calculate how many clusters
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected in TSNE.")

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
    return cluster_result

def JC_hierarchical_clustering(df, cutoff=0.45):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    matrix = (matrix > 0.5).astype(int) 
    my_linkage = linkage(matrix, method='average', metric='jaccard')
    cluster_labels = fcluster(my_linkage, t=cutoff, criterion='distance')
    n_clusters = len(set(cluster_labels))
    print (n_clusters, "clusters detected in JC.")

    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(cluster_labels)):
            if cluster_labels[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    return cluster_result

def hierarchical_clustering(df, cutoff=1.6):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())

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
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    return cluster_result

def UMAP(df):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())
    # print (matrix)

    try:
        X_embedded = umap.UMAP().fit_transform(matrix, 
                                               n_neighbors=1,
                                               min_dist=0.1,
                                               metric='euclidean')
    except Exception as e:
        print(f"Failed to create UMAP: {e}")
        return
    
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
    return cluster_result

def PCA_plot(df):
    print ("start PCA...")
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    ## transpose the profiles
    # profiles = profiles.T
    ## add pseudo values to zero values
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())
    X_embedded = pca.fit_transform(matrix)

    ## cluster the pca result
    clustering = DBSCAN(eps=0.1, min_samples=1).fit(X_embedded)
    # print (clustering.labels_)
    ## save the cluster result in a dataframe, and plot it like in tse function
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected after PCA.")
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(clustering.labels_)):
            if clustering.labels_[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    ## save the cluster result
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    return cluster_result


if __name__ == "__main__":
    # Example usage:
    # cluster_csv = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_NM3/motif_cluster.h.csv"  # Path to the cluster CSV file
    # fasta_path = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa"  # Path to the input FASTA file
    # output_prefix = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_NM3/bins/bin"  # Prefix for output files

    # cluster_csv = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/motif_cluster.csv"  # Path to the cluster CSV file
    # fasta_path = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_200ppm_r1_LR_scaffold.fa"  # Path to the input FASTA file
    # output_prefix = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/bins/bin"  # Prefix for output files

    # bin_contigs_to_fastas(cluster_csv, fasta_path, output_prefix)

    whole_ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    bin_dir = "/home/shuaiw/borg/bench/soil/run2/test_bin/"
    profile_file =  "/home/shuaiw/borg/bench/soil/run2/motif_profile.csv"
    min_frac = 0.4

    ## build the dir if it is not exist
    if not os.path.exists(bin_dir):
        os.makedirs(bin_dir)


    profiles = pd.read_csv(profile_file, index_col=0)
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    output_prefix = os.path.join(bin_dir, "bin")
    
    # Filter columns where any value is greater than 0.5
    profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]
    print ("filtered shape", profiles.shape)
    print (profiles)
    # cluster_df = TSE(profiles)
    cluster_df = JC_hierarchical_clustering(profiles)
    bin_contigs_to_fastas_df(cluster_df, whole_ref, output_prefix)

