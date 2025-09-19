import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity
import re

from Bio.Seq import Seq

def read_motif_freq(motif_freq_file, prefix, all_motif_freq):
    motif_freq_df = pd.read_csv(motif_freq_file)
    ## make the column percentage_of_profile * 0.01
    motif_freq_df["percentage_of_profile"] = motif_freq_df["percentage_of_profile"] * 0.01
    ## only keep the rows with percentage_of_profile > 0.01
    motif_freq_df = motif_freq_df[motif_freq_df["percentage_of_profile"] > 0.04]
    ## remove the rows with motifString's and start end letter is not in GATC
    motif_freq_df = motif_freq_df[motif_freq_df["motifString"].str.endswith(('G', 'A', 'T', 'C'))]
    motif_freq_df = motif_freq_df[motif_freq_df["motifString"].str.startswith(('G', 'A', 'T', 'C'))]
    ## percentage_of_profile > 0.1, set it as 0.1
    # motif_freq_df.loc[motif_freq_df["percentage_of_profile"] > 0.1, "percentage_of_profile"] = 0.1
    motif_freq_df["prefix"] = prefix
    for i, row in motif_freq_df.iterrows():
        all_motif_freq.append([row["percentage_of_profile"], prefix, row["motifString"]])
    return motif_freq_df

def collect_freq():
    all_motif_freq = []
    all_dir = "/home/shuaiw/borg/paper/run2/"
    i = 0
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if re.search("sludge", prefix):
            continue
        print (f"Processing {prefix}...")
        work_dir = f"{all_dir}/{prefix}/{prefix}_methylation3"
        fai = f"{all_dir}/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
        map_sum = f"{all_dir}/{prefix}/{prefix}.align.count.csv"
        all_host_file = f"{all_dir}/{prefix}/all_host_ctgs.tsv"
        depth_file = os.path.join(work_dir, "mean_depth.csv")
        host_sum_file = os.path.join(work_dir, "host_summary.csv")
        orphan_file = os.path.join(work_dir, "regulatory_motif_enrichment.csv")
        motif_freq_file = os.path.join(work_dir, "motif_length_stats.csv")

        
        if not os.path.exists(fai):
            print(f"Skipping {prefix} as fai file does not exist.")
            continue
        ## skip if all_host_file does not exist
        if not os.path.exists(all_host_file):
            print(f"Skipping {prefix} as all_host_file does not exist.")
            continue

        motif_freq_df = read_motif_freq(motif_freq_file, prefix, all_motif_freq)
        i += 1
        # if i > 10:
        #     break
    ## convert to dataframe
    all_motif_freq_df = pd.DataFrame(all_motif_freq, columns=["frequency", "prefix", "motif"])
    print (all_motif_freq_df)
    plot_motif_freq(all_motif_freq_df)

def plot_motif_freq(all_motif_freq_df):
    ## plot clustered heatmap
    motif_pivot = all_motif_freq_df.pivot(index='prefix', columns='motif', values='frequency').fillna(0)
    
    # Print data statistics to understand the range
    print(f"Data shape: {motif_pivot.shape}")
    print(f"Data range: {motif_pivot.min().min():.4f} to {motif_pivot.max().max():.4f}")
    print(f"Data mean: {motif_pivot.mean().mean():.4f}")
    print("Sample of data:")
    print(motif_pivot.head())
    
    # # Normalize the data for better visualization
    motif_pivot = motif_pivot.div(motif_pivot.sum(axis=0), axis=1)
    # Create a clustered heatmap without row standardization
    # plt.figure(figsize=(15, 30))
    g = sns.clustermap(motif_pivot, cmap='viridis', metric='cosine', method='average', 
                       standard_scale=1, figsize=(30, 10), 
                       xticklabels=True, yticklabels=True,  # Ensure labels are shown
                       row_cluster=True, col_cluster=True,  # Keep clustering
                       dendrogram_ratio=(0.01, 0.01))  # Very small dendrograms instead of 0
    
    # Hide the dendrograms by making them invisible
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    
    # Rotate x-axis labels for better readability
    g.ax_heatmap.tick_params(axis='x', labelsize=12)
    g.ax_heatmap.tick_params(axis='y', rotation=45, labelsize=12)
    
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)
    
    # Remove the colorbar/legend
    g.cax.set_visible(False)
    
    plt.title('Motif Frequency Heatmap')
    plt.savefig("../../tmp/results2/motif_frequency_heatmap.pdf", dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    # meta_file = "/home/shuaiw/Methy/assembly_pipe/prefix_table.tab"
    # sample_env_dict = read_metadata(meta_file)
    collect_freq()