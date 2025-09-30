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
    # motif_freq_df = motif_freq_df[motif_freq_df["percentage_of_profile"] > 0.04]
    ## remove the rows with motifString's and start end letter is not in GATC
    motif_freq_df = motif_freq_df[motif_freq_df["motifString"].str.endswith(('G', 'A', 'T', 'C'))]
    motif_freq_df = motif_freq_df[motif_freq_df["motifString"].str.startswith(('G', 'A', 'T', 'C'))]
    ## percentage_of_profile > 0.1, set it as 0.1
    # motif_freq_df.loc[motif_freq_df["percentage_of_profile"] > 0.1, "percentage_of_profile"] = 0.1
    motif_freq_df["prefix"] = prefix
    for i, row in motif_freq_df.iterrows():
        all_motif_freq.append([row["percentage_of_profile"], prefix, row["motifString"]])
    return motif_freq_df, all_motif_freq

def read_motif_freq_ctg(profile, all_motif_freq, mge_dict):
    df = pd.read_csv(profile)
    ## get the column names
    cols = df.columns.tolist()
    for i, row in df.iterrows():
        # print (row)
        for j in range(1, len(cols)):
            if cols[j] in mge_dict:
                ctg_type = "MGE"
            else:
                ctg_type = "Host"
            all_motif_freq.append([row[cols[j]], cols[j], row["motifString"], ctg_type])
    return all_motif_freq

def read_gtdb(gtdb):
    ## collect the species name of each contig, if no species name, use genus name, if no genus name, use family name
    gtdb_dict = {}
    if not os.path.exists(gtdb):
        print (f"[⚠️] GTDB file not found: {gtdb}")
        return gtdb_dict
    df = pd.read_csv(gtdb, sep="\t")
    for i, row in df.iterrows():
        classification = row["classification"]
        contig = row["user_genome"]
        fields = classification.split(";")

        # print (contig, classification)
        if len(fields) < 1:
            taxon = "Unclassified"
        else:
            for i in range(-1, -len(fields)-1, -1):
                if fields[i].strip():
                    if len(fields[i].strip()) > 3:
                        taxon = fields[i].strip()

                        break
        # print (contig, taxon)
        gtdb_dict[contig] = taxon   
    return gtdb_dict

## read length of each contig
def read_fai(fai):
    contig_len = {}
    if not os.path.exists(fai):
        print (f"[⚠️] FAI file not found: {fai}")
        return contig_len
    with open(fai, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            contig = fields[0]
            length = int(fields[1])
            contig_len[contig] = length
    return contig_len

def check_pure(gtdb_dict, contig_len_dict, mge_dict):
    """
    Check if the assembly is pure based on GTDB classification.
    """
    # 1) Length-weighted GTDB consistency (core test)

    # Make a table with: contig, length, GTDB_species (or genus if species is unavailable), and any confidence/support score.

    # Compute the fraction of total assembly length assigned to each species.

    # Quick call:

    # Pure: ≥99% of total assembly length → one species.

    # Likely pure: 95–99% one species and the leftovers are short contigs (<10–20 kb) that look like plasmid/phage.

    # Mixed: ≥2 species each >5% of total length (or any minority block >1 Mb from a different genus).
    species_len = {}
    total_len = 0
    for contig, taxon in gtdb_dict.items():
        if contig in mge_dict:  ## skip the contigs that are MGEs
            continue
        length = contig_len_dict.get(contig, 0)
        total_len += length
        if taxon not in species_len:
            species_len[taxon] = 0
        species_len[taxon] += length
    species_fraction = {taxon: length / total_len for taxon, length in species_len.items()}
    # print (species_fraction)      
    sorted_species = sorted(species_fraction.items(), key=lambda x: x[1], reverse=True)
    print (sorted_species)
    if len(sorted_species) == 0:
        return "unknown"
    top_species, top_fraction = sorted_species[0]
    if top_fraction >= 0.99:
        return "pure"
    elif top_fraction >= 0.95:
        # check if the other species are short contigs
        for taxon, fraction in sorted_species[1:]:
            if fraction * total_len > 20000:
                return "mixed"
        return "likely pure"
    else:
        for taxon, fraction in sorted_species[1:]:
            if fraction > 0.05:
                return "mixed"
            if fraction * total_len > 1000000:
                return "mixed"
        return "likely pure"

def check_pure2(checkm):
    ## read checkm file
    if not os.path.exists(checkm):
        print (f"[⚠️] CheckM file not found: {checkm}")
        return "unknown"
    df = pd.read_csv(checkm, sep="\t")
    if "Completeness" not in df.columns or "Contamination" not in df.columns:
        print (f"[⚠️] CheckM file does not have required columns.")
        return "unknown"
    completeness = df["Completeness"].values[0]
    contamination = df["Contamination"].values[0]
    if contamination <= 5:
        return "pure"
    else:
        return "mixed"

def read_MGE(mge_file):
    mge_dict = {}
    f = open(mge_file, "r")
    for line in f:
        if line.startswith("seq_name"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 2:
            continue
        contig = fields[0]
        mge_type = fields[1]
        mge_dict[contig] = mge_type
    return mge_dict

def read_host(host_sum_file, min_dp = 5):
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    all_df = pd.DataFrame()
    
   


    

def collect_freq():
    all_motif_freq = []
    all_link_df = pd.DataFrame()
    all_dir = "/home/shuaiw/borg/paper/isolation/bacteria/"
    i = 0
    pure_num = 0
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if re.search("sludge", prefix):
            continue
        # if prefix != "ERR10820930":
        #     continue
        # print (f"Processing {prefix}...")
        work_dir = f"{all_dir}/{prefix}/{prefix}_methylation"
        fai = f"{all_dir}/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
        map_sum = f"{all_dir}/{prefix}/{prefix}.align.count.csv"
        all_host_file = f"{all_dir}/{prefix}/all_host_ctgs.tsv"
        depth_file = os.path.join(work_dir, "mean_depth.csv")
        host_sum_file = os.path.join(work_dir, "host_summary.csv")
        orphan_file = os.path.join(work_dir, "regulatory_motif_enrichment.csv")
        motif_freq_file = os.path.join(work_dir, "motif_length_stats.csv")
        profile = os.path.join(work_dir, "motif_profile.csv")
        gtdb = os.path.join(work_dir, "../GTDB_2/gtdbtk.bac120.summary.tsv")
        checkm = os.path.join(work_dir, "../checkM2/quality_report.tsv")
        mge_file = f"{all_dir}/{prefix}/all_mge.tsv"
        gtdb_dict = read_gtdb(gtdb)
        contig_len_dict = read_fai(fai)
        mge_dict = read_MGE(mge_file)
        # pure_flag = check_pure(gtdb_dict, contig_len_dict, mge_dict)
        pure_flag = check_pure2(checkm)
        if pure_flag == "pure":
            pure_num += 1
        else:
            continue
        # print (f"{prefix} is {pure_flag}")
        if len(mge_dict) < 1:
            continue

        if not os.path.exists(motif_freq_file):
            print(f"Skipping {prefix} as no motifs.")
            continue
        ## if host_sum_file is not empty
        line_num = sum(1 for line in open(host_sum_file) if line.strip())
        if line_num > 1:
            host_sum = pd.read_csv(host_sum_file) 
            all_link_df = pd.concat([all_link_df, host_sum], axis=0)

        all_motif_freq = read_motif_freq_ctg(profile, all_motif_freq, mge_dict)
        i += 1
        # if i > 10:
        #     break
    ## convert to dataframe
    # all_motif_freq_df = pd.DataFrame(all_motif_freq, columns=["frequency", "prefix", "motif", "contig_type"])
    # print (all_motif_freq_df)
    # plot_motif_freq(all_motif_freq_df)
    # print (f"Total {i} samples, {pure_num} pure samples.")

    analyze_link(all_link_df)

def plot_motif_freq(all_motif_freq_df):
    ## plot clustered heatmap
    motif_pivot = all_motif_freq_df.pivot(index='prefix', columns='motif', values='frequency').fillna(0)
    
    # Create contig type annotation for rows
    # Get the contig type for each contig (prefix)
    contig_type_df = all_motif_freq_df[['prefix', 'contig_type']].drop_duplicates().set_index('prefix')
    
    # Create a color map for contig types
    contig_type_colors = {'Host': '#1f77b4', 'MGE': '#ff7f0e'}  # Blue for Host, Orange for MGE
    row_colors = contig_type_df['contig_type'].map(contig_type_colors)
    
    # Print data statistics to understand the range
    print(f"Data shape: {motif_pivot.shape}")
    print(f"Data range: {motif_pivot.min().min():.4f} to {motif_pivot.max().max():.4f}")
    print(f"Data mean: {motif_pivot.mean().mean():.4f}")
    print("Sample of data:")
    print(motif_pivot.head())
    print(f"Contig types: {contig_type_df['contig_type'].value_counts()}")
    
    # Clean the data to handle edge cases that cause NaN in clustering
    # Remove columns (motifs) that are all zeros
    motif_pivot = motif_pivot.loc[:, motif_pivot.sum(axis=0) > 0]
    
    # Remove rows (samples) that are all zeros
    motif_pivot = motif_pivot.loc[motif_pivot.sum(axis=1) > 0, :]
    
    # Align row colors with the cleaned data
    row_colors = row_colors.loc[motif_pivot.index]
    
    # Add small epsilon to avoid zero vectors (which cause issues with cosine distance)
    epsilon = 1e-10
    motif_pivot = motif_pivot + epsilon
    
    # Normalize by row (each sample sums to 1)
    # motif_pivot = motif_pivot.div(motif_pivot.sum(axis=1), axis=0)
    
    # Check for any remaining NaN values
    if motif_pivot.isna().any().any():
        print("Warning: Found NaN values after preprocessing, filling with zeros")
        motif_pivot = motif_pivot.fillna(0)
    
    print(f"After preprocessing - Data shape: {motif_pivot.shape}")
    print(f"After preprocessing - Data range: {motif_pivot.min().min():.6f} to {motif_pivot.max().max():.6f}")
    
    # Create a clustered heatmap without row standardization (since we already normalized)
    # Use euclidean distance instead of cosine to avoid numerical issues
    g = sns.clustermap(motif_pivot, cmap='viridis', metric='euclidean', method='average', 
                       standard_scale=None, figsize=(20, 20),  # No additional standardization
                       xticklabels=True, yticklabels=True,  # Ensure labels are shown
                       row_cluster=True, col_cluster=True,  # Keep clustering
                       dendrogram_ratio=(0.01, 0.01),  # Very small dendrograms
                       cbar_pos=(0.92, 0.3, 0.03, 0.4),  # Position colorbar outside: (x, y, width, height)
                       row_colors=row_colors,  # Add row annotation for contig types
                       colors_ratio=0.01)  # Minimize the width of the row color bar
    
    # # Hide the dendrograms by making them invisible
    # g.ax_row_dendrogram.set_visible(False)
    # g.ax_col_dendrogram.set_visible(False)
    
    # Rotate x-axis labels for better readability
    g.ax_heatmap.tick_params(axis='x', labelsize=12)
    g.ax_heatmap.tick_params(axis='y', rotation=45, labelsize=12)
    
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # Remove the colorbar/legend
    g.cax.set_visible(False)
    
    # Colorbar is now positioned outside via cbar_pos parameter
    # Add a label to the colorbar
    g.cax.set_ylabel('Motif Frequency', rotation=270, labelpad=20)
    
    # Create a custom legend for contig types
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=contig_type_colors['Host'], label='Host'),
                      Patch(facecolor=contig_type_colors['MGE'], label='MGE')]
    g.ax_heatmap.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))
    
    plt.title('Motif Frequency Heatmap')
    plt.savefig("../../tmp/results2/motif_freq_iso.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def analyze_link(all_link_df):
    ## keep these with depth > 10
    # all_link_df = all_link_df[all_link_df['MGE_cov'] > 10]
    all_link_df = all_link_df[all_link_df['MGE_len'] > 10000]
    ## plot the distribution of final_score in all_link_df
    plt.figure(figsize=(10, 6))
    sns.histplot(all_link_df['final_score'], bins=50, kde=True, color='skyblue')
    plt.title('Distribution of Final Scores', fontsize=16, fontweight='bold')
    plt.xlabel('Final Score', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.savefig("../../tmp/results2/final_score_distribution.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    ## output these with final_score == 0
    zero_score_df = all_link_df[all_link_df['final_score'] == 0]
    for  i, row in zero_score_df.iterrows():
        print (row["MGE"], row["host"], row["final_score"])
    
    ## three subplots, one is scatter plot of MGE_cov vs host_cov, one is MGE_gc vs host_gc, one is distribution of cos_sim
    print (len(all_link_df), "total links")
    fig, axs = plt.subplots(1, 3, figsize=(21, 7))
    # Scatter plot of MGE_cov vs host_cov
    axs[0].scatter(all_link_df['MGE_cov'], all_link_df['host_cov'], alpha=0.6, color='teal')
    axs[0].set_title('MGE Coverage vs Host Coverage', fontsize=16, fontweight='bold')
    axs[0].set_xlabel('MGE Coverage', fontsize=14)
    axs[0].set_ylabel('Host Coverage', fontsize=14)
    axs[0].grid(axis='both', alpha=0.3)
    # Set equal limits and add diagonal line for coverage plot
    cov_min = min(all_link_df['MGE_cov'].min(), all_link_df['host_cov'].min())
    cov_max = max(all_link_df['MGE_cov'].max(), all_link_df['host_cov'].max())
    axs[0].set_xlim(cov_min, cov_max)
    axs[0].set_ylim(cov_min, cov_max)
    axs[0].plot([cov_min, cov_max], [cov_min, cov_max], 'k--', alpha=0.7, linewidth=1.5, label='x = y')
    axs[0].legend()
    
    # Scatter plot of MGE_gc vs host_gc
    axs[1].scatter(all_link_df['MGE_gc'], all_link_df['host_gc'], alpha=0.6, color='coral')
    axs[1].set_title('MGE GC% vs Host GC%', fontsize=16, fontweight='bold')
    axs[1].set_xlabel('MGE GC%', fontsize=14)
    axs[1].set_ylabel('Host GC%', fontsize=14)
    axs[1].grid(axis='both', alpha=0.3)
    # Set equal limits and add diagonal line for GC plot
    gc_min = min(all_link_df['MGE_gc'].min(), all_link_df['host_gc'].min())
    gc_max = max(all_link_df['MGE_gc'].max(), all_link_df['host_gc'].max())
    axs[1].set_xlim(gc_min, gc_max)
    axs[1].set_ylim(gc_min, gc_max)
    axs[1].plot([gc_min, gc_max], [gc_min, gc_max], 'k--', alpha=0.7, linewidth=1.5, label='x = y')
    axs[1].legend()
    # Distribution of cos_sim
    sns.histplot(all_link_df['cos_sim'], bins=50, kde=True, color='orchid', ax=axs[2])
    axs[2].set_title('Distribution of Cosine Similarity', fontsize=16, fontweight='bold')
    axs[2].set_xlabel('Cosine Similarity', fontsize=14)
    axs[2].set_ylabel('Count', fontsize=14)
    axs[2].grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig("../../tmp/results2/link_features.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    for i, row in all_link_df.iterrows():
        if row["cos_sim"] < 0.5:
            print (row["MGE"], row["host"], row["MGE_gc"], row["host_gc"])


if __name__ == "__main__":
    # meta_file = "/home/shuaiw/Methy/assembly_pipe/prefix_table.tab"
    # sample_env_dict = read_metadata(meta_file)
    collect_freq()