import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity
import re
import sys
from Bio.Seq import Seq

from sample_object import My_sample, Isolation_sample

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

def quantify_sharing(sample_obj, mge_dict, depth_dict, length_dict, unique_motifs):
    """
    Add contig annotation columns (type, depth, length) to motif profile data.
    
    Args:
        profile (str): Path to motif profile CSV file
        mge_dict (dict): Dictionary mapping contig names to MGE types
        depth_dict (dict): Dictionary mapping contig names to depth values
        length_dict (dict): Dictionary mapping contig names to length values
    
    Returns:
        pd.DataFrame: Enhanced dataframe with additional annotation columns
    """
    print (sample_obj.profile)
    df = pd.read_csv(sample_obj.profile)
    
    # Get contig names (all columns except the first one which is motifString)
    contig_cols = df.columns.tolist()[1:]  # Skip the motifString column
    
    # Create annotation dictionaries for each contig
    contig_annotations = []
    
    for contig in contig_cols:
        # Determine contig type
        if contig in mge_dict:
            ctg_type = "MGE"
            mge_type = mge_dict[contig]
        else:
            ctg_type = "Host"
            mge_type = "N/A"
        
        # Get depth and length
        depth = depth_dict.get(contig, 0.0)
        length = length_dict.get(contig, 0)
        
        contig_annotations.append({
            'contig': contig,
            'contig_type': ctg_type,
            'mge_type': mge_type,
            'depth': depth,
            'length': length,
            'prefix': sample_obj.prefix,
            'phylum': sample_obj.phylum,
            'lineage': sample_obj.lineage
        })
    
    # Convert to DataFrame
    annotations_df = pd.DataFrame(contig_annotations)
    
    # Transform the motif profile from wide to long format
    df_melted = df.melt(id_vars=['motifString'], 
                        var_name='contig', 
                        value_name='motif_frequency')
    
    # Merge with annotations
    df_enhanced = df_melted.merge(annotations_df, on='contig', how='left')
    
    # Reorder columns for better readability
    # df_enhanced = df_enhanced[['motifString', 'contig', 'contig_type', 'mge_type', 
    #                            'depth', 'length', 'motif_frequency']]
    ## only consider the motif in unique_motifs
    df_enhanced = df_enhanced[df_enhanced['motifString'].isin(unique_motifs)]

    return cal_jaccard(df_enhanced, sample_obj.depth_cutoff, sample_obj.length_cutoff)
    # return df_enhanced

def cal_jaccard(df_enhanced, depth_cutoff = 10, length_cutoff = 5000, bin_freq = 0.3):
    ## filter the dataframe based on depth and length cutoff
    filtered_df = df_enhanced[
        (df_enhanced['depth'] >= depth_cutoff) &
        (df_enhanced['length'] >= length_cutoff)
    ]
    ## binary the motif_frequency based on bin_freq
    filtered_df['motif_frequency'] = filtered_df['motif_frequency'].apply(lambda x: 1 if x >= bin_freq else 0)
    # print (filtered_df)
    ## calculate jaccard similarity between each pair of MGE and host
    ## Only consider motifs with frequency = 1 (present motifs)
    present_motifs = filtered_df[filtered_df['motif_frequency'] == 1]
    
    mge_motifs = present_motifs[present_motifs['contig_type'] == 'MGE'][['motifString', 'contig', 'prefix', 'phylum', 'lineage']]
    host_motifs = present_motifs[present_motifs['contig_type'] == 'Host'][['motifString', 'contig', 'prefix', 'phylum', 'lineage']]

    # print(f"[✔] Found {len(mge_motifs)} MGE motifs and {len(host_motifs)} Host motifs with frequency = 1")
    
    jaccard_scores = []
    for mge_contig in mge_motifs['contig'].unique():
        mge_motifs_set = set(mge_motifs[mge_motifs['contig'] == mge_contig]['motifString'])
        for host_contig in host_motifs['contig'].unique():
            host_motifs_set = set(host_motifs[host_motifs['contig'] == host_contig]['motifString'])
            intersection = mge_motifs_set.intersection(host_motifs_set)
            union = mge_motifs_set.union(host_motifs_set)
            if len(union) == 0:
                jaccard = 0.0
            else:
                jaccard = len(intersection) / len(union)

            if mge_motifs['prefix'].values[0] == host_motifs['prefix'].values[0] :
                relation = 'same_sample'
            else:
                relation = esti_relation(
                    mge_motifs[mge_motifs['contig'] == mge_contig]['lineage'].values[0],
                    host_motifs[host_motifs['contig'] == host_contig]['lineage'].values[0]
                )
            ## count number of motifs specifically in mge
            mge_specific = mge_motifs_set - host_motifs_set
            mge_specific_num = len(mge_specific)

            jaccard_scores.append({
                'mge_contig': mge_contig,
                'host_contig': host_contig,
                'jaccard_similarity': jaccard,
                'share_num': len(intersection),
                'all_num': len(union),
                'relation': relation,
                'mge_specific_num': mge_specific_num,
            })
    # print (jaccard_scores)
    return pd.DataFrame(jaccard_scores), present_motifs

def cross_sample_jaccard(present_motifs_all, random_ctg_num=500):
    mge_motifs = present_motifs_all[present_motifs_all['contig_type'] == 'MGE'][['motifString', 'contig', 'prefix', 'phylum', 'lineage']]
    host_motifs = present_motifs_all[present_motifs_all['contig_type'] == 'Host'][['motifString', 'contig', 'prefix', 'phylum', 'lineage']]
    host_contig_list = host_motifs['contig'].unique()
    print(f"[✔] Found {len(mge_motifs)} MGE motifs and {len(host_motifs)} Host motifs with frequency = 1")
    
    jaccard_scores = []
    for mge_contig in mge_motifs['contig'].unique():
        mge_motifs_set = set(mge_motifs[mge_motifs['contig'] == mge_contig]['motifString'])
        ## randomly select 100 host contigs
        # for host_contig in np.random.choice(host_contig_list, size=min(random_ctg_num, len(host_contig_list)), replace=False):
        for host_contig in host_contig_list:
            host_motifs_set = set(host_motifs[host_motifs['contig'] == host_contig]['motifString'])
            intersection = mge_motifs_set.intersection(host_motifs_set)
            union = mge_motifs_set.union(host_motifs_set)
            if len(union) == 0:
                jaccard = 0.0
            else:
                jaccard = len(intersection) / len(union)

            # Get the prefix for the current MGE and host contigs
            mge_prefix = mge_motifs[mge_motifs['contig'] == mge_contig]['prefix'].iloc[0]
            host_prefix = host_motifs[host_motifs['contig'] == host_contig]['prefix'].iloc[0]
            
            if mge_prefix == host_prefix:
                relation = 'same_sample'
                # continue
            else:
                relation = esti_relation(
                    mge_motifs[mge_motifs['contig'] == mge_contig]['lineage'].values[0],
                    host_motifs[host_motifs['contig'] == host_contig]['lineage'].values[0]
                )
            jaccard_scores.append({
                'mge_contig': mge_contig,
                'host_contig': host_contig,
                'jaccard_similarity': jaccard,
                'share_num': len(intersection),
                'all_num': len(union),
                'relation': relation,
            })

    return pd.DataFrame(jaccard_scores)

def esti_relation(lineage1, lineage2):
    """
    Estimate the relationship between two lineages based on taxonomic levels.
    
    Args:
        lineage1 (str): Taxonomic lineage of the first contig
        lineage2 (str): Taxonomic lineage of the second contig
    
    Returns:
        str: Estimated relationship (e.g., 'same_species', 'same_genus', etc.)
    """
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ## either lineage is unknown, return different_lineage

    if re.search(r'Unclassified Bacteria', lineage1) or re.search(r'Unclassified Bacteria', lineage2):
        return 'different_lineage'

    taxa1 = lineage1.split(';')
    taxa2 = lineage2.split(';')
    # print ("###", lineage1, lineage2)
    
    for i in range(len(levels)-1, -1, -1):
        if taxa1[i] == taxa2[i] and taxa1[i] != 'unknown' and len(taxa1[i]) > 3:
            return f'same_{levels[i]}'
    
    return 'different_lineage'


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
    plt.savefig("../../tmp/results2/motif_freq_iso2.pdf", dpi=300, bbox_inches="tight")
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

def plot_jaccard(jaccard_all, fig_dir):
    ## box plot hue to phylum
    plt.figure(figsize=(5, 5))
    sns.boxplot(data=jaccard_all, x='phylum', y='jaccard_similarity', hue='phylum', palette='Set3', legend=False)
    # plt.title('Jaccard Similarity between MGE and Host by Phylum', fontsize=16, fontweight='bold')
    plt.xlabel('Phylum', fontsize=14)
    plt.ylabel('Jaccard Similarity', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_by_phylum.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def gradient_plot(jaccard_all, fig_dir):
    ## plot the proportion of MGE-host pairs with at least 0, 1, 2, ... shared motifs  using bar plot
    max_share = jaccard_all['share_num'].max()
    share_counts = []
    for i in range(0, int(max_share) + 1):
        count = len(jaccard_all[jaccard_all['share_num'] >= i])
        share_counts.append({'shared_motifs': i, 'count': count, 'proportion': count / len(jaccard_all)})
    share_counts_df = pd.DataFrame(share_counts)
    plt.figure(figsize=(5, 5))
    sns.barplot(data=share_counts_df, x='shared_motifs', y='proportion', color='steelblue')
    # plt.title('Number of MGE-Host Pairs by Shared Motifs', fontsize=16, fontweight='bold')
    plt.xlabel('Number of Shared Motifs', fontsize=14)
    plt.ylabel('Proportion of MGE-Host Pairs', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/mge_host_shared_motifs.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def cross_taxa_plot(jaccard_all, fig_dir):
    ## box plot hue to relation with ordered x-axis
    # Define the desired order from most closely related to least related
    relation_order = ['same_sample', 'same_species', 'same_genus', 'same_family', 'same_order', 
                     'same_class', 'same_phylum', 'same_domain']
    
    plt.figure(figsize=(5, 5))
    sns.boxplot(data=jaccard_all, x='relation', y='jaccard_similarity', 
                hue='relation', palette='Set2', legend=False, order=relation_order)
    # plt.title('Jaccard Similarity between MGE and Host by Taxonomic Relation', fontsize=16, fontweight='bold')
    plt.xlabel('Taxonomic Relation', fontsize=14)
    plt.ylabel('Jaccard Similarity', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_by_taxa_relation.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def count_jaccard(same_sample_df, jaccard_all_sample):
    ## print all rows with jaccard similarity < 0.5 for check
    low_jaccard = same_sample_df[same_sample_df['jaccard_similarity'] < 0.5]
    print("Low Jaccard Similarity Samples:")
    for index, row in low_jaccard.iterrows():
        print(f"Prefix: {row['prefix']}, MGE Contig: {row['mge_contig']}, Host Contig: {row['host_contig']}, Jaccard Similarity: {row['jaccard_similarity']:.4f}")
    
    ## count proportion of MGE-host pairs with jaccard similarity =1
    jaccard_1 = same_sample_df[same_sample_df['jaccard_similarity'] == 1]
    print(f"Proportion of MGE-host pairs with Jaccard similarity = 1: {len(jaccard_1) / len(same_sample_df):.4f}")
    ## compare cross-sample and same-sample jaccard similarity using boxplot, jaccard_all_sample have both same_sample and cross_sample
    combined_df = pd.concat([same_sample_df.assign(sample_type='same_sample'),
                             jaccard_all_sample[jaccard_all_sample['relation'] != 'same_sample'].assign(sample_type='cross_sample')])
    plot_compare(combined_df, fig_dir)
    ## print the rows with mge_specific_num > 0
    mge_specific = same_sample_df[same_sample_df['mge_specific_num'] > 0]
    print("MGE-specific motifs in same-sample pairs:")
    for index, row in mge_specific.iterrows():
        print(f"Prefix: {row['prefix']}, MGE Contig: {row['mge_contig']}, Host Contig: {row['host_contig']}, MGE-specific Motifs: {row['mge_specific_num']}")

def plot_compare(combined_df, fig_dir):
    plt.figure(figsize=(5, 5))
    sns.boxplot(data=combined_df, x='sample_type', y='jaccard_similarity', hue='sample_type', palette='Set1', legend=False)
    plt.xlabel('Sample Type', fontsize=14)
    plt.ylabel('Jaccard Similarity', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    ## add statistical annotation
    from scipy.stats import mannwhitneyu
    same_sample_values = combined_df[combined_df['sample_type'] == 'same_sample']['jaccard_similarity']
    cross_sample_values = combined_df[combined_df['sample_type'] == 'cross_sample']['jaccard_similarity']
    stat, p_value = mannwhitneyu(same_sample_values, cross_sample_values, alternative='two-sided')
    ## add p-value to the plot
    plt.text(0.5, max(combined_df['jaccard_similarity']) * 0.95, f'p-value = {p_value:.4e}', ha='center', fontsize=12)
    plt.savefig(f"{fig_dir}/jaccard_similarity_same_vs_cross_sample.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def plot_jaccard_distribution(same_sample_df, fig_dir):
    ## plot the distribution of jaccard similarity using histogram
    plt.figure(figsize=(5, 5))
    sns.histplot(same_sample_df['jaccard_similarity'], bins=20, kde=True, color='lightgreen')
    # plt.title('Distribution of Jaccard Similarity (Same Sample)', fontsize=16, fontweight='bold')
    plt.xlabel('Jaccard Similarity', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_distribution_same_sample.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def main(all_dir, fig_dir):
    i = 0
    pure_num, mge_num, has_motif_num = 0, 0, 0
    jaccard_all = pd.DataFrame()
    present_motifs_all = pd.DataFrame()
    for prefix in os.listdir(all_dir):

        sample_obj = Isolation_sample(prefix, all_dir)
        if not os.path.exists(sample_obj.profile):
            continue

        sample_obj.get_phylum()
        mge_dict = sample_obj.read_MGE()
        depth_dict, length_dict = sample_obj.read_depth()
        pure_flag = sample_obj.check_pure2()
        unique_motif_num, unique_motifs = sample_obj.get_unique_motifs()
        sample_obj.explore_specific_motifs()

        ### continue if line number < 2 sample_obj.profile
        line_count = sum(1 for line in open(sample_obj.profile)) - 1  # Subtract 1 for header
        if line_count < 2:
            print(f"Skipping {prefix} as insufficient motif data.")
            continue
        
        i += 1
        if pure_flag == "pure":
            pure_num += 1
        else:
            continue
        # print (f"{prefix} is {pure_flag}")
        if len(mge_dict) < 1:
            continue
        else:
            mge_num += 1

        if unique_motifs == 0:
            print(f"Skipping {prefix} as no motifs.")
            continue
        else:
            has_motif_num += 1
        

        
        # Enhanced motif sharing analysis with additional annotations
        jaccard_scores, present_motifs = quantify_sharing(sample_obj, mge_dict, depth_dict, length_dict, unique_motifs)
        jaccard_scores['prefix'] = prefix
        jaccard_scores['phylum'] = sample_obj.phylum
        
        # Use robust concatenation that handles empty DataFrames
        if jaccard_all.empty:
            jaccard_all = jaccard_scores.copy()
        else:
            jaccard_all = pd.concat([jaccard_all, jaccard_scores], axis=0, ignore_index=True)
            
        if present_motifs_all.empty:
            present_motifs_all = present_motifs.copy()
        else:
            present_motifs_all = pd.concat([present_motifs_all, present_motifs], axis=0, ignore_index=True)

    
    same_sample_df = jaccard_all[jaccard_all['relation'] == 'same_sample']
    print (same_sample_df)
    # analyze_link(all_link_df)
    plot_jaccard(same_sample_df, fig_dir)
    gradient_plot(same_sample_df, fig_dir)
    

    jaccard_all_sample = cross_sample_jaccard(present_motifs_all)
    cross_taxa_plot(jaccard_all_sample, fig_dir)
    
    count_jaccard(same_sample_df, jaccard_all_sample)
    plot_jaccard_distribution(same_sample_df, fig_dir)
    print (f"Total {i} samples, {pure_num} pure samples, {mge_num} samples with MGEs, {has_motif_num} samples with motifs.")

if __name__ == "__main__":
    # meta_file = "/home/shuaiw/Methy/assembly_pipe/prefix_table.tab"
    # sample_env_dict = read_metadata(meta_file)
    # all_dir = "/home/shuaiw/borg/paper/isolation/bacteria/"
    all_dir = "/home/shuaiw/borg/paper/isolation//batch2_results/"
    fig_dir = "../../tmp/figures/motif_sharing/"
    main(all_dir, fig_dir)