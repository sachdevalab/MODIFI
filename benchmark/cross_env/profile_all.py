import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import plotly.graph_objects as go
from collections import defaultdict
import matplotlib.patches as mpatches
import random
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
from Bio.Seq import Seq
import profile
from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys, os, re
import numpy as np
from collections import defaultdict
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.metrics import jaccard_score
from sklearn.metrics.pairwise import cosine_similarity
import pickle
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster


sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'specificity'))
from sample_object import get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa
from check_motif_change import get_motif_sites, get_modified_ratio, read_ipd_ratio, read_ref
from profile_good_ctgs import read_metadata



def collect_all_motifs(min_frac=0.4, min_sites = 100):
    all_motif_df = pd.DataFrame()
    
    # ctg_taxa_dict = get_ctg_taxa(all_dir)
    for prefix in os.listdir(all_dir):

        sample_obj = My_sample(prefix, all_dir)
        motif_df = sample_obj.simple_load_motifs()
        if motif_df is None:
            continue
        all_motif_df = pd.concat([all_motif_df, motif_df], axis=0)
    df_motif = all_motif_df[(all_motif_df['fraction'] >= min_frac) & (all_motif_df['nDetected'] >= min_sites)]
    ## rm redundant motifs which are reverse complement

    unique_motifs_df = pd.DataFrame()
    unique_motifs = []
    for index, row in df_motif.iterrows():
        if row['motifString'] not in unique_motifs and  str(Seq(row['motifString']).reverse_complement()) not in unique_motifs:
            unique_motifs.append(row['motifString'])
            unique_motifs_df = pd.concat([unique_motifs_df, pd.DataFrame([row])], axis=0)
    print (f"Total unique motifs collected: {unique_motifs_df.shape[0]}")
    print (unique_motifs_df)
    return unique_motifs_df
    
def collect_contigs(all_dir, sample_env_dict, ctg_taxa_dict):
    all_data = []
    all_base_data = []
    genome_data = []
    all_genome_list = []
    ctg_taxa_dict = get_ctg_taxa(all_dir)
    for my_dir in os.listdir(all_dir):
        prefix = my_dir

        if re.search("sludge", prefix):
            continue
        # print (f"Processing {prefix}...")


        sample_obj = My_sample(prefix, all_dir)
        ## skip if not exist depth file
        if not os.path.exists(sample_obj.depth_file):
            # print (f"[⚠️] Depth file not found for sample {prefix}, skipping...")
            continue
        sample_obj.read_depth()
        best_ctgs = sample_obj.get_final_best_ctg()
        # print ("best ctgs num:", len(best_ctgs))
        genome_list , contig_list= sample_obj.get_high_dp_ctg_list()
        environment = sample_env_dict.get(prefix, "unknown")
        for contig in contig_list:
            ctg_lineage = ctg_taxa_dict[contig] if contig in ctg_taxa_dict else "Unknown"
            ctg_phylum = classify_taxa(ctg_lineage, level='phylum')
            all_genome_list.append([prefix, contig, environment, ctg_phylum])
        # if len(all_genome_list) > 500:
        #     break
    print (f"Total genomes collected: {len(all_genome_list)}")
    return all_genome_list

def start_profile(unique_motifs_df, all_genome_list, datafile):
    data = []
    for prefix, contig, environment, ctg_phylum in all_genome_list:
        ctg_obj = My_contig(prefix, all_dir, contig, "meta")
        ## skip if files do not exist
        if not os.path.exists(ctg_obj.gff) or not os.path.exists(ctg_obj.ipd_ratio_file) \
            or not os.path.exists(ctg_obj.ctg_ref):
            print (f"[⚠️] Missing files for contig {contig} in sample {prefix}, \
                    {ctg_obj.gff}, {ctg_obj.ipd_ratio_file}, {ctg_obj.ctg_ref} skipping...")
            continue
        REF = read_ref(ctg_obj.ctg_ref)
        # print (REF)
        modified_loci = get_modified_ratio(ctg_obj.gff)
        # motifs = pd.read_csv(all_motifs)
        # ipd_ratio_dict = read_ipd_ratio(ctg_obj.ipd_ratio_file)

        # Process motifs sequentially
        for index, row in unique_motifs_df.iterrows():
            motif_new = row['motifString']
            exact_pos = row["centerPos"]
            motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
            data.append([contig, motif_new + "_" + str(exact_pos), motif_profile[-2], environment, ctg_phylum])
    df = pd.DataFrame(data, columns = ["contig", "motifString", "fraction", "environment", "phylum"])
    df.to_csv(datafile, index=False, sep='\t')
    
    print(f"Created motif profile data with shape: {df.shape}")
    print(f"Unique environments: {df['environment'].unique()}")
    print(f"Unique phylums: {df['phylum'].unique()}")
    return df
    
def plot_SNE(df, fig_dir):
    """
    Create t-SNE scatter plot with colors indicating environment
    """
    # Create pivot table with contigs as rows and motifs as columns
    df_pivot = df.pivot(index='contig', columns='motifString', values='fraction').fillna(0)
    
    if df_pivot.shape[0] < 2:
        print("Not enough data points for t-SNE analysis")
        return
    
    print(f"Data shape for t-SNE: {df_pivot.shape}")
    
    # Get environment and phylum info for each contig
    contig_info = df.groupby('contig').agg({
        'environment': 'first',
        'phylum': 'first'
    }).reset_index()
    
    # Prepare data matrix
    matrix = df_pivot.values
    
    # Replace zero values with small random values to avoid issues
    # mask = matrix == 0
    # matrix[mask] = np.random.uniform(-0.01, 0.01, mask.sum())
    
    try:
        # Apply t-SNE
        print("Running t-SNE...")
        tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(matrix)-1))
        X_embedded = tsne.fit_transform(matrix)
        
        # Create DataFrame with t-SNE results
        tsne_df = pd.DataFrame({
            'contig': df_pivot.index,
            'tsne1': X_embedded[:, 0],
            'tsne2': X_embedded[:, 1]
        })
        
        # Merge with environment and phylum info
        tsne_df = tsne_df.merge(contig_info, on='contig')
        
        # Create the plot
        plt.figure(figsize=(12, 10))
        
        # Get unique environments and create color palette
        environments = tsne_df['environment'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(environments)))
        env_colors = dict(zip(environments, colors))
        
        # Plot points colored by environment
        for env in environments:
            env_data = tsne_df[tsne_df['environment'] == env]
            plt.scatter(env_data['tsne1'], env_data['tsne2'], 
                       c=[env_colors[env]], label=env, alpha=0.7, s=60)
        
        plt.xlabel('t-SNE Component 1')
        plt.ylabel('t-SNE Component 2')
        plt.title('t-SNE Analysis of Methylation Motif Profiles\nColored by Environment')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        # Save the plot
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, 'tsne_environment.pdf'), bbox_inches='tight')
        plt.close()  # Close the figure to avoid display
        
        # Also create a plot colored by phylum
        plt.figure(figsize=(12, 10))
        phylums = tsne_df['phylum'].unique()
        phylum_colors = plt.cm.Set1(np.linspace(0, 1, len(phylums)))
        phylum_color_dict = dict(zip(phylums, phylum_colors))
        
        for phylum in phylums:
            phylum_data = tsne_df[tsne_df['phylum'] == phylum]
            plt.scatter(phylum_data['tsne1'], phylum_data['tsne2'], 
                       c=[phylum_color_dict[phylum]], label=phylum, alpha=0.7, s=60)
        
        plt.xlabel('t-SNE Component 1')
        plt.ylabel('t-SNE Component 2')
        plt.title('t-SNE Analysis of Methylation Motif Profiles\nColored by Phylum')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        plt.savefig(os.path.join(fig_dir, 'tsne_phylum.pdf'), bbox_inches='tight')
        plt.close()  # Close the figure to avoid display
        
        # Save the t-SNE results
        tsne_df.to_csv(os.path.join(fig_dir, 'tsne_results.csv'), index=False)
        
        print(f"t-SNE plots saved to {fig_dir}")
        print(f"Environment distribution: {tsne_df['environment'].value_counts()}")
        print(f"Phylum distribution: {tsne_df['phylum'].value_counts()}")
        
    except Exception as e:
        print(f"Failed to create t-SNE: {e}")
        import traceback
        traceback.print_exc()
        return

def create_combined_analysis(df, fig_dir):
    """
    Create a combined analysis with PCA and t-SNE side by side
    """
    # Create pivot table
    df_pivot = df.pivot(index='contig', columns='motifString', values='fraction').fillna(0)
    
    if df_pivot.shape[0] < 2:
        print("Not enough data points for combined analysis")
        return
    
    # Get environment info
    contig_info = df.groupby('contig').agg({
        'environment': 'first',
        'phylum': 'first'
    }).reset_index()
    
    matrix = df_pivot.values
    
    # Replace zeros with small random values
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.01, 0.01, mask.sum())
    
    try:
        # Apply PCA
        pca = PCA(n_components=2, random_state=42)
        X_pca = pca.fit_transform(matrix)
        
        # Apply t-SNE
        tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(matrix)-1))
        X_tsne = tsne.fit_transform(matrix)
        
        # Create subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        
        # Get unique environments and colors
        environments = contig_info['environment'].unique()
        colors = plt.cm.Set3(np.linspace(0, 1, len(environments)))
        env_colors = dict(zip(environments, colors))
        
        # PCA plot
        for i, env in enumerate(environments):
            env_indices = contig_info[contig_info['environment'] == env].index
            ax1.scatter(X_pca[env_indices, 0], X_pca[env_indices, 1], 
                       c=[env_colors[env]], label=env, alpha=0.7, s=60)
        
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%} variance)')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.2%} variance)')
        ax1.set_title('PCA Analysis\nColored by Environment')
        ax1.legend()
        
        # t-SNE plot
        for i, env in enumerate(environments):
            env_indices = contig_info[contig_info['environment'] == env].index
            ax2.scatter(X_tsne[env_indices, 0], X_tsne[env_indices, 1], 
                       c=[env_colors[env]], label=env, alpha=0.7, s=60)
        
        ax2.set_xlabel('t-SNE Component 1')
        ax2.set_ylabel('t-SNE Component 2')
        ax2.set_title('t-SNE Analysis\nColored by Environment')
        ax2.legend()
        
        plt.tight_layout()
        
        # Save combined plot
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, 'pca_tsne_combined.pdf'), bbox_inches='tight')
        plt.close()  # Close the figure to avoid display
        
        print(f"Combined PCA/t-SNE analysis saved to {fig_dir}")
        
    except Exception as e:
        print(f"Failed to create combined analysis: {e}")
        import traceback


if __name__ == "__main__":
    all_dir = "/home/shuaiw/borg/paper/run2/"
    meta_file = "/home/shuaiw/mGlu/assembly_pipe/prefix_table.tab"
    fig_dir = "../../tmp/figures/multi_env_linkage/"
    datafile = os.path.join(fig_dir, "all_motif_profile.tsv")
    sample_env_dict = read_metadata(meta_file)
    ctg_taxa_dict = get_ctg_taxa(all_dir)

    unique_motifs_df = collect_all_motifs(min_frac=0.8, min_sites=1000)
    all_genome_list = collect_contigs(all_dir, sample_env_dict, ctg_taxa_dict)
    ## show genome number per environment, and motif number
    print(f"Number of motifs: {unique_motifs_df.shape[0]}")
    df = start_profile(unique_motifs_df, all_genome_list, datafile)

    # df = pd.read_csv(datafile, sep='\t')

    print(f"Data loaded with shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")
    
    plot_SNE(df, fig_dir)
    create_combined_analysis(df, fig_dir)
