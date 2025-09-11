#!/usr/bin/env python3

"""
Enhanced motif dereplication script with similarity tracking and clustering
"""

import pandas as pd
import os
import sys
import argparse
from collections import defaultdict
from Bio.Seq import Seq
import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns

def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences of equal length"""
    if len(seq1) != len(seq2):
        return float('inf')
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def sequence_similarity(seq1, seq2):
    """Calculate sequence similarity (1 - normalized hamming distance)"""
    if len(seq1) != len(seq2):
        return 0.0
    ham_dist = hamming_distance(seq1, seq2)
    return 1.0 - (ham_dist / len(seq1))

def reverse_complement_similarity(seq1, seq2):
    """Calculate similarity including reverse complement"""
    # Direct similarity
    sim1 = sequence_similarity(seq1, seq2)
    
    # Reverse complement similarity
    rev_comp = str(Seq(seq2).reverse_complement())
    sim2 = sequence_similarity(seq1, rev_comp)
    
    return max(sim1, sim2)

def create_similarity_matrix(motifs):
    """Create similarity matrix for motifs"""
    n = len(motifs)
    similarity_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i, n):
            if i == j:
                similarity_matrix[i, j] = 1.0
            else:
                sim = reverse_complement_similarity(motifs[i], motifs[j])
                similarity_matrix[i, j] = sim
                similarity_matrix[j, i] = sim
    
    return similarity_matrix

def cluster_motifs(motifs, similarity_threshold=0.8):
    """Cluster motifs based on similarity threshold"""
    if len(motifs) <= 1:
        return [0] * len(motifs)
    
    # Create similarity matrix
    sim_matrix = create_similarity_matrix(motifs)
    
    # Convert to distance matrix
    distance_matrix = 1 - sim_matrix
    
    # Perform clustering with updated parameters
    try:
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=1-similarity_threshold,
            metric='precomputed',
            linkage='average'
        )
    except TypeError:
        # Fallback for older scikit-learn versions
        clustering = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=1-similarity_threshold,
            affinity='precomputed',
            linkage='average'
        )
    
    cluster_labels = clustering.fit_predict(distance_matrix)
    return cluster_labels

def read_motifs(motif_file, min_frac):
    """Read motifs from CSV file and return motif dictionary and samples"""
    df = pd.read_csv(motif_file)
    ## each row is a motif, the first column is the motif name, the rest are the samples

    ## remove the rows where all fractions are below min_frac
    numeric_cols = df.columns[1:]  # All columns except first
    df = df[(df[numeric_cols] > min_frac).any(axis=1)]
    
    ## remove the columns where all fractions are below min_frac
    # Vectorized column filtering
    cols_to_keep = [df.columns[0]]  # Always keep first column
    for col in numeric_cols:
        if (df[col] > min_frac).any():
            cols_to_keep.append(col)
    
    df = df[cols_to_keep]
    print ("filtered shape 1", df.shape)

    # Vectorized creation of motif dictionary
    motif_dict = {}
    samples = list(df.columns)[1:]
    
    # Convert to dictionary more efficiently
    for motif, *values in df.values:
        motif_dict[motif] = values
    
    return motif_dict, samples

def calculate_cosine_similarity_matrix(motif_dict):
    """Calculate cosine similarity matrix between motifs"""
    motifs = list(motif_dict.keys())
    # Create matrix more efficiently
    values_matrix = np.array([motif_dict[motif] for motif in motifs], dtype=np.float32)
    
    # Calculate cosine similarity (this is already optimized in sklearn)
    cos_sim_matrix = cosine_similarity(values_matrix)
    
    return cos_sim_matrix, motifs

def plot_cosine_similarity_heatmap(cos_sim_matrix, motifs, output_dir):
    """Plot cosine similarity heatmap with clustering"""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    
    plt.figure(figsize=(14, 12))
    
    # Convert similarity to distance matrix for clustering
    distance_matrix = 1 - cos_sim_matrix
    
    # Perform hierarchical clustering if we have enough motifs
    if len(motifs) > 2:
        try:
            # Convert to condensed distance matrix for linkage
            condensed_dist = squareform(distance_matrix)
            
            # Perform clustering
            row_linkage = linkage(condensed_dist, method='average')
            col_linkage = row_linkage  # Same for both axes since matrix is symmetric
            
            # Create clustered heatmap
            g = sns.clustermap(cos_sim_matrix, 
                             row_linkage=row_linkage,
                             col_linkage=col_linkage,
                             xticklabels=motifs,
                             yticklabels=motifs,
                             cmap='viridis',
                             annot=False,
                             square=True,
                             cbar_kws={'label': 'Cosine Similarity'},
                             figsize=(14, 12))
            
            # Rotate labels for better readability
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), rotation=45, ha='right')
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), rotation=0)
            
            plt.suptitle('Clustered Motif Cosine Similarity Matrix', y=0.95)
            
            heatmap_path = os.path.join(output_dir, 'cosine_similarity_heatmap_clustered.png')
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Clustered cosine similarity heatmap saved to {heatmap_path}")
            
        except Exception as e:
            print(f"Clustering failed ({e}), creating simple heatmap...")
            # Fallback to simple heatmap
            plt.figure(figsize=(12, 10))
            sns.heatmap(cos_sim_matrix, 
                        xticklabels=motifs,
                        yticklabels=motifs,
                        cmap='viridis',
                        annot=False,
                        square=True,
                        cbar_kws={'label': 'Cosine Similarity'})
            
            plt.title('Motif Cosine Similarity Matrix')
            plt.xlabel('Motifs')
            plt.ylabel('Motifs')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            
            heatmap_path = os.path.join(output_dir, 'cosine_similarity_heatmap.png')
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Simple cosine similarity heatmap saved to {heatmap_path}")
    else:
        # For very few motifs, just create simple heatmap
        sns.heatmap(cos_sim_matrix, 
                    xticklabels=motifs,
                    yticklabels=motifs,
                    cmap='viridis',
                    annot=True,
                    square=True,
                    cbar_kws={'label': 'Cosine Similarity'})
        
        plt.title('Motif Cosine Similarity Matrix')
        plt.xlabel('Motifs')
        plt.ylabel('Motifs')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        
        heatmap_path = os.path.join(output_dir, 'cosine_similarity_heatmap.png')
        plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Cosine similarity heatmap saved to {heatmap_path}")

def count_motif_occurrences(df, fai, min_frac=0.3):
    """
    Record the length of each contig and count the total length of contigs 
    that have this motif (with frac > min_frac)
    
    Args:
        df: DataFrame with motif profile data (columns: motifString, contig, fraction)
        fai: Path to FASTA index file (.fai)
        min_frac: Minimum fraction threshold (default: 0.3)
    
    Returns:
        tuple: (motif_total_lengths dict, total_profile_length, profile_stats dict)
    """
    # Read contig lengths from .fai file
    contig_lengths = {}
    with open(fai, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_name = parts[0]
                contig_length = int(parts[1])
                contig_lengths[contig_name] = contig_length
    
    ### count the total length for all contigs that show up in the motif profile
    all_contigs_in_profile = set(df['contig'].unique())
    total_profile_length = 0
    profile_contigs_found = 0
    
    for contig in all_contigs_in_profile:
        if contig in contig_lengths:
            total_profile_length += contig_lengths[contig]
            profile_contigs_found += 1
    
    print(f"Total length of all contigs in motif profile: {total_profile_length:,} bp")
    print(f"Number of contigs in profile: {len(all_contigs_in_profile)}")
    print(f"Contigs found in .fai file: {profile_contigs_found}")

    
    # Count total length for each motif (vectorized approach)
    # Filter data once for fraction threshold
    filtered_df = df[df['fraction'] > min_frac].copy()
    
    # Add contig lengths column for efficient lookup
    filtered_df['contig_length'] = filtered_df['contig'].map(contig_lengths)
    
    # Remove rows where contig length is not found (NaN)
    filtered_df = filtered_df.dropna(subset=['contig_length'])
    
    # Group by motif and sum contig lengths
    motif_total_lengths = filtered_df.groupby('motifString')['contig_length'].sum().to_dict()
    
    # Prepare summary statistics
    profile_stats = {
        'total_profile_length': total_profile_length,
        'total_contigs_in_profile': len(all_contigs_in_profile),
        'contigs_found_in_fai': profile_contigs_found,
        'missing_contigs': len(all_contigs_in_profile) - profile_contigs_found
    }
    
    return motif_total_lengths, total_profile_length, profile_stats

def motif_cluster_worker(motif_file, fai, output_dir, min_frac=0.3, similarity_threshold=0.7):

    # Default behavior for testing
    # prefix = "cow_bioreactor_5"
    # motif_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/motif_profile.csv"
    # fai = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
    # output_dir = "../tmp/results2/"

    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read motifs
    print(f"Reading motif data from {motif_file}")
    motif_dict, samples = read_motifs(motif_file, min_frac)
    
    print(f"Found {len(motif_dict)} motifs across {len(samples)} samples")
    
    # Calculate cosine similarity
    print("Calculating cosine similarity matrix...")
    cos_sim_matrix, motifs = calculate_cosine_similarity_matrix(motif_dict)
    
    # Plot cosine similarity heatmap
    # plot_cosine_similarity_heatmap(cos_sim_matrix, motifs, output_dir)
    
    # Cluster motifs by sequence similarity
    print("Clustering motifs by sequence similarity...")
    cluster_labels = cluster_motifs(motifs, similarity_threshold)
    n_clusters = len(np.unique(cluster_labels))
    print(f"Found {n_clusters} sequence similarity clusters")
    ## output the elements in each cluster
    cluster_dict = defaultdict(list)
    for i, label in enumerate(cluster_labels):
        cluster_dict[label].append(motifs[i])
    cluster_output_path = os.path.join(output_dir, 'motif_clusters.txt')
    with open(cluster_output_path, 'w') as f:
        for cluster_id, motif_list in cluster_dict.items():
            f.write(f"Cluster {cluster_id}:\n")
            for motif in motif_list:
                f.write(f"  {motif}\n")
            f.write("\n")
    print(f"Motif clusters saved to {cluster_output_path}")
    
    # Count motif occurrences by contig length
    print("Counting motif occurrences by contig length...")
    # Create DataFrame more efficiently using vectorized operations
    motif_data = []
    for motif, values in motif_dict.items():
        for i, fraction in enumerate(values):
            if fraction > 0:  # Only include non-zero fractions
                motif_data.append([motif, samples[i], fraction])
    
    if motif_data:
        # Create DataFrame in one go
        motif_occurrence_df = pd.DataFrame(motif_data, columns=['motifString', 'contig', 'fraction'])
        motif_lengths, total_profile_length, profile_stats = count_motif_occurrences(motif_occurrence_df, fai, min_frac)
        
        # Save motif length statistics
        length_stats_path = os.path.join(output_dir, 'motif_length_stats.csv')
        length_df = pd.DataFrame(list(motif_lengths.items()), columns=['motifString', 'total_length'])
        length_df = length_df.sort_values('total_length', ascending=False)
        
        # Add percentage of total profile length
        length_df['percentage_of_profile'] = (length_df['total_length'] / total_profile_length * 100).round(2)
        
        # Add cluster labels to the dataframe
        motif_to_cluster = {motifs[i]: cluster_labels[i] for i in range(len(motifs))}
        length_df['cluster_id'] = length_df['motifString'].map(motif_to_cluster)
        
        # Reorder columns to have cluster_id near the beginning
        length_df = length_df[['motifString', 'cluster_id', 'total_length', 'percentage_of_profile']]
        
        length_df.to_csv(length_stats_path, index=False)
        print(f"Motif length statistics saved to {length_stats_path}")
        
        # Print top 10 motifs by total contig length
        print("\nTop 10 motifs by total contig length:")
        for i, (motif, cluster_id, length, percentage) in enumerate(length_df.head(10)[['motifString', 'cluster_id', 'total_length', 'percentage_of_profile']].values):
            print(f"{i+1:2d}. {motif} (Cluster {cluster_id}): {length:,} bp ({percentage}% of profile)")
        
        print(f"\nProfile Summary:")
        print(f"Total profile length: {total_profile_length:,} bp")
        print(f"Total contigs in profile: {profile_stats['total_contigs_in_profile']}")
        print(f"Contigs found in .fai: {profile_stats['contigs_found_in_fai']}")
        if profile_stats['missing_contigs'] > 0:
            print(f"Missing contigs: {profile_stats['missing_contigs']}")
            
    else:
        print("No motif occurrence data found for length analysis")
    return length_df
    

if __name__ == "__main__":
    motif_cluster_worker()