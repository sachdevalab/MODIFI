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

def out_best_ctgs(ref, best_ref, best_ctgs):
    ## use biopython to extract the best contigs from the reference fasta file
    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    with open(ref, "r") as f_in, open(best_ref, "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            if record.id in best_ctgs:
                SeqIO.write(record, f_out, "fasta")
    print(f"Extracted {len(best_ctgs)} best contigs to {best_ref}")

def count_motifs(depth_file, best_ctgs, work_dir):
    ## read depth file
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    for index, row in depth_df.iterrows():
        if row['depth'] >= 10:
            good_depth[row['contig']] = row['depth']
    best_depth_ctg = []
    for ctg in best_ctgs:
        if ctg in good_depth:
            best_depth_ctg.append(ctg)
    print(f"Total {len(best_depth_ctg)} contigs with depth >= 10 found.")
    has_motif_ctg_num = 0
    motif_num_list = []
    for ctg in best_depth_ctg:
        motif_file = os.path.join(work_dir,"motifs", f"{ctg}.motifs.csv")
        # print (motif_file)
        if os.path.exists(motif_file):
            df_motif = pd.read_csv(motif_file)
            ## only keep themotifs with fraction >= 0.4, and nDetected >=100
            df_motif = df_motif[(df_motif['fraction'] >= 0.4) & (df_motif['nDetected'] >= 100)]
            if not df_motif.empty:
                has_motif_ctg_num += 1
                motif_num_list.append(len(df_motif))
                pass
        else:
            print(f"No motifs found for {ctg}.")
    print (f"Total {has_motif_ctg_num} contigs with motifs found in the best contigs with depth >= 10.")
    print (has_motif_ctg_num/ len(best_depth_ctg) * 100, "% of the best contigs with depth >= 10 have motifs.")
    ## print mean and median of motif numbers
    print(f"Mean motif number: {np.mean(motif_num_list)}")
    print(f"Median motif number: {np.median(motif_num_list)}")
    print(f"Max motif number: {np.max(motif_num_list)}")
    print(f"Min motif number: {np.min(motif_num_list)}")
    ## plot the distribution of motif numbers
    plt.figure(figsize=(10, 6))
    plt.hist(motif_num_list, bins=50, color='blue', alpha=0.7)
    plt.title('Distribution of Motif Numbers in Best Contigs')
    plt.xlabel('Number of Motifs')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join("../../tmp/results", "motif_num_distribution.png"), dpi=300, bbox_inches='tight')


fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
best_ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.circular.contigs.fa"
profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"

work_dir = "/home/shuaiw/borg/bench/soil/run1/"
depth_file = os.path.join(work_dir, "mean_depth.csv")
best_ctgs = get_best_ctg()
count_motifs(depth_file, best_ctgs, work_dir)
# out_best_ctgs(ref, best_ref, best_ctgs)
# plot_heatmap(profile_file)