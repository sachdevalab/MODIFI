import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
import re
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import networkx as nx

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
from sample_object import get_unique_motifs, My_sample, get_ctg_taxa, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
from check_motif_change import given_species_drep, given_species_drep_fast
from find_borg import find_assembly

class Borg_Entry:
    def __init__(self, line):
        """Parse a single line from the BORG file"""
        fields = line.strip().split('\t')
        if len(fields) >= 9:
            self.borg_ref = fields[0]
            self.type = fields[1]
            self.seq_name = fields[2]
            self.identity = float(fields[3])
            self.query_coverage = float(fields[4])
            self.target_coverage = float(fields[5])
            self.alignment_length = int(fields[6])
            self.length = int(fields[7])
            self.ctg_depth = float(fields[9])
            self.sample_name = fields[8]
        else:
            raise ValueError(f"Invalid line format: {line}")

    def __str__(self):
        return f"BORG: {self.borg_ref} -> {self.seq_name} (type: {self.type}, identity: {self.identity}, query_cov: {self.query_coverage}, target_cov: {self.target_coverage}, align_len: {self.alignment_length}, length: {self.length}, depth: {self.ctg_depth})"

class My_Borg:
    def __init__(self, borg_file):
        self.borg_file = borg_file
        self.borg_entries = []
        self.load_borg_data()

    def load_borg_data(self):
        """Load BORG data from file into a list of Borg_Entry objects"""
        with open(self.borg_file, 'r') as f:
            # Skip header line
            header = f.readline()
            for line in f:
                if line.strip():  # Skip empty lines
                    try:
                        ## skip the header line
                        if line.startswith("borg_ref"):
                            continue
                        borg_entry = Borg_Entry(line)
                        self.borg_entries.append(borg_entry)
                    except ValueError as e:
                        print(f"Warning: {e}")

    def get_high_depth_borgs(self, min_depth=10.0):
        """Get BORG entries with contig depth above threshold"""
        return [entry for entry in self.borg_entries if entry.ctg_depth >= min_depth]

def personal_plot(profile_df):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist
    # Add color bar for BORG/Host indication
    from matplotlib.patches import Rectangle
    
    # Create pivot table for heatmap
    pivot_df = profile_df.pivot(index='contig', columns='motifString', values='fraction')
    pivot_df = pivot_df.fillna(0)  # Fill NaN values with 0

    # Create a mapping from contig to BORG_Ref BEFORE clustering
    contig_to_borg = profile_df.groupby('contig')['BORG_Ref'].first().to_dict()

    # Perform hierarchical clustering on both rows and columns

    # Cluster rows (contigs)
    row_linkage = linkage(pdist(pivot_df, metric='correlation'), method='ward')
    row_dendro = dendrogram(row_linkage, no_plot=True)
    row_order = row_dendro['leaves']

    # Cluster columns (motifs)
    col_linkage = linkage(pdist(pivot_df.T, metric='correlation'), method='ward')
    col_dendro = dendrogram(col_linkage, no_plot=True)
    col_order = col_dendro['leaves']

    # Reorder pivot table based on clustering
    pivot_df = pivot_df.iloc[row_order, col_order]

    
    # Create labels with BORG reference information AFTER reordering
    y_labels = []
    for contig in pivot_df.index:  # Use reordered index
        # Use the pre-created mapping instead of filtering the dataframe
        # borg_ref = contig_to_borg.get(contig, 'NA')
        if contig in contig_to_borg:
            borg_ref = contig_to_borg[contig]
        else:
            print (f"Warning: {contig} not found in contig_to_borg mapping.")
        label = f"{contig},{borg_ref}"
        y_labels.append(label)
    
    # Update the pivot_df index with the new labels before creating heatmap
    pivot_df.index = y_labels
    
    # Plot independent dendrogram tree
    fig_tree = plt.figure(figsize=(13, 13))
    ax_tree = fig_tree.add_subplot(111)
    dendro_tree = dendrogram(row_linkage, orientation='right', ax=ax_tree,
                            labels=y_labels, color_threshold=0, leaf_font_size=8,
                            no_plot=False, count_sort=False, distance_sort=False)
    
    # Get the leaf positions from dendrogram
    leaf_positions = dendro_tree['leaves']
    icoord = np.array(dendro_tree['icoord'])
    dcoord = np.array(dendro_tree['dcoord'])
    
    # Normalize all horizontal distances to create uniform branch lengths (cladogram style)
    # Redraw the dendrogram with uniform spacing
    ax_tree.clear()
    
    # Create cladogram by setting all branches to unit length
    def make_cladogram(linkage_matrix):
        """Convert linkage matrix to have uniform branch lengths"""
        Z = linkage_matrix.copy()
        # Set all distances to 1 for uniform branch lengths
        Z[:, 2] = 1
        return Z
    
    cladogram_linkage = make_cladogram(row_linkage)
    dendro_tree = dendrogram(cladogram_linkage, orientation='right', ax=ax_tree,
                            labels=y_labels, color_threshold=0, leaf_font_size=8)
    
    # Get the leaf positions from dendrogram
    leaf_positions = dendro_tree['leaves']
    
    # Draw horizontal lines from leaves to labels for better visibility
    for i, (leaf_idx, y_pos) in enumerate(zip(leaf_positions, range(5, len(leaf_positions)*10+5, 10))):
        # Get the x position of the leaf (rightmost point)
        leaf_x = 0
        ax_tree.plot([leaf_x, -2.0], [y_pos, y_pos], 'k--', linewidth=0.5, alpha=0.5, dashes=(5, 5))
    
    ax_tree.set_xlabel("Distance", fontsize=14)
    ax_tree.set_yticklabels(ax_tree.get_yticklabels(), fontsize=7, ha='right')
    ax_tree.set_title(f"{cluster_species} - Hierarchical Clustering", fontsize=16, pad=15)
    ax_tree.spines['top'].set_visible(False)
    ax_tree.spines['right'].set_visible(False)
    ax_tree.spines['left'].set_visible(False)
    plt.tight_layout()
    tree_plot_name = plot_name.replace('.pdf', '_tree.pdf')
    plt.savefig(tree_plot_name, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved tree plot to: {tree_plot_name}")
    
    # Create figure with subplots for dendrogram and heatmap
    fig = plt.figure(figsize=(22, 15))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 10], wspace=0.02)
    
    # Plot dendrogram on the left - with no_plot first to get the proper ordering
    ax_dendro = fig.add_subplot(gs[0])
    # Use the same linkage but plot it correctly oriented
    dendro_plot = dendrogram(row_linkage, orientation='left', ax=ax_dendro, 
                             no_labels=True, color_threshold=0)
    ax_dendro.set_xticks([])
    ax_dendro.set_yticks([])
    ax_dendro.spines['top'].set_visible(False)
    ax_dendro.spines['right'].set_visible(False)
    ax_dendro.spines['bottom'].set_visible(False)
    ax_dendro.spines['left'].set_visible(False)
    ax_dendro.invert_yaxis()  # Invert to match heatmap orientation
    
    # Plot heatmap on the right
    ax_heatmap = fig.add_subplot(gs[1])
    sns.heatmap(pivot_df, cmap="YlGnBu", annot=False, fmt=".2f", 
                cbar_kws={'label': 'Fraction'}, ax=ax_heatmap, yticklabels=True, xticklabels=True)
    
    # Set labels and formatting
    ax_heatmap.set_xlabel("Motif String", fontsize=12)
    ax_heatmap.set_ylabel("Contig", fontsize=12)
    ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=90, ha='right', fontsize=6)
    ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), rotation=0, fontsize=6)
    ax_heatmap.set_title(f"{cluster_species}", fontsize=14, pad=20)
    

    

    plt.tight_layout()
    plt.savefig(plot_name, dpi=300, bbox_inches='tight')
    plt.show()

def collect_host_genus(all_dir, fasta_dict):
    members  = []
    anno_dict = {}
    non_Mp_members = []
    non_Mp_anno_dict = {}
    ctg_taxa_dict = get_ctg_taxa(all_dir)
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if not re.search("soil", prefix):
            continue
        reference = fasta_dict[prefix] if prefix in fasta_dict else None
        print (f"Processing {prefix}...")
        sample_obj = My_sample(prefix, all_dir)
        sample_obj.reference_fasta = reference
        sample_obj.fai = sample_obj.reference_fasta + ".fai" if reference else None
        sample_obj.read_depth()
        sample_obj.get_len_dict()
        genome_list, contig_list =  sample_obj.get_high_dp_ctg_list( min_depth=5, min_len=100000)
        
        for ctg in contig_list:
            ctg_obj = My_contig(prefix, all_dir, ctg)
            ctg_lineage = ctg_taxa_dict[ctg] if ctg in ctg_taxa_dict else "Unknown"
            # print (ctg_lineage)
            ctg_phylum = classify_taxa(ctg_lineage, level='phylum')
            ctg_genus = classify_taxa(ctg_lineage, level='genus')
            # ctg_len = ctg_obj.get_ctg_len()
            taxon = get_detail_taxa_name(ctg_lineage)
            # print (ctg_genus)
            if ctg_genus == "g__Methanoperedens":
                # print (ctg, ctg_lineage)
                # os.system(f"cp {ctg_obj.ctg_ref} /home/shuaiw/borg/paper/borg_data/align/refs")
                members.append((ctg, prefix))
                anno_dict[ctg] = taxon
            elif ctg_genus != "g__" and not re.search("Unclassified", ctg_lineage):
                non_Mp_members.append((ctg, prefix))
                non_Mp_anno_dict[ctg] = taxon
    print (f"Total host contigs collected: {len(members)}")
    return members, anno_dict, non_Mp_members, non_Mp_anno_dict

def count_mod_freq(all_dir, borg_anno_dict, score_cutoff = 30):
    data = []
    for ctg in borg_anno_dict:
        prefix = "_".join(ctg.split('_')[:-2])
        print (f"Processing {ctg} from sample {prefix}...")
        ctg_obj = My_contig(prefix, all_dir, ctg)
        if not os.path.exists(ctg_obj.ctg_ref):
            print (f"[!] fasta file not found for {ctg}, skipping...")
            continue
        ctg_type, ctg_anno = borg_anno_dict[ctg]
        ctg_obj.get_mod_ratio(score_cutoff)
        data.append([ctg, prefix, ctg_type, ctg_anno, ctg_obj.ctg_len, ctg_obj.modified_num, ctg_obj.modified_motif_num, 
                     ctg_obj.modified_ratio, ctg_obj.modified_motif_ratio, ctg_obj.motif_ratio])
    df = pd.DataFrame(data, columns=['contig', 'sample', 'type', 'annotation', 'length', 'modified_num', 'modified_motif_num', 
                                     'modified_ratio', 'modified_motif_ratio', 'motif_ratio'])
    ## sort by modified_ratio descending
    df = df.sort_values(by='modified_ratio', ascending=False)
    print (df)
    return df

def get_unique_motif(motifs_to_keep):
    ## only keep one for each reverse complementary motif pair
    unique_motifs = set()
    unique_motif_ids = set()
    for motif_id in motifs_to_keep:
        motif_string = motif_id.split("_")[0]
        rev_comp = str(Seq(motif_string).reverse_complement())
        if rev_comp not in unique_motifs:
            unique_motifs.add(motif_string)
            unique_motif_ids.add(motif_id)
    return list(unique_motif_ids)

def extract_sample_name(contig_name):
    parts = str(contig_name).split('_')
    if len(parts) >= 8:
        return '_'.join(parts[:8])
    else:
        return contig_name

def umap_plot(profile_df, plot_name):


    # Add sample column extracted from contig name (first 8 underscore-separated parts)


    
    
    ## Randomly retain 5 Non-Mp contigs
    non_mp_df = profile_df[profile_df['Genome'] == 'Non-Mp']
    other_df = profile_df[profile_df['Genome'] != 'Non-Mp']
    
    if len(non_mp_df) > 0:
        # Get unique Non-Mp contigs and randomly sample 5
        unique_non_mp_contigs = non_mp_df['contig'].unique()
        np.random.seed(42)
        selected_contigs = np.random.choice(unique_non_mp_contigs, size=min(50, len(unique_non_mp_contigs)), replace=False)
        non_mp_retained = non_mp_df[non_mp_df['contig'].isin(selected_contigs)]
        profile_df = pd.concat([other_df, non_mp_retained], ignore_index=True)
        print(f"Retained {len(selected_contigs)} Non-Mp contigs out of {len(unique_non_mp_contigs)}")
    else:
        profile_df = other_df
        print("No Non-Mp contigs found")

    # Create pivot table for dimensionality reduction
    pivot_df = profile_df.pivot(index='contig', columns='motifString', values='fraction')
    pivot_df = pivot_df.fillna(0)  # Fill NaN values with 0
    


    print("\n=== Creating network graph based on Pearson correlation ===")

    # Calculate Pearson correlation matrix
    similarity_matrix = pivot_df.T.corr(method='pearson').values


    ## get a jaccard matrix as well with binary cutoff of 0.5
    from sklearn.metrics import pairwise_distances
    # Binarize the matrix with cutoff 0.5
    pivot_binary = (pivot_df >= 0.4).astype(int)
    # # # Calculate Jaccard distance and convert to similarity
    # jaccard_distances = pairwise_distances(pivot_binary.values, metric='jaccard')
    # similarity_matrix = 1 - jaccard_distances

    ## get  another matrix recording the number of shared motifs 
    # Calculate shared motif counts
    # similarity_matrix = np.dot(pivot_binary.values, pivot_binary.values.T)


    # for motif in pivot_df.columns:
    #     val1 = pivot_df.at['SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_155_C', motif]
    #     val2 = pivot_df.at['SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_METAMDBG_724567_L', motif]
    #     print (motif, round(val1,3), round(val2,3))


    # Get genome type for each contig
    if 'Genome' in profile_df.columns:
        genome_per_contig = profile_df.groupby('contig')['Genome'].first()
        contig_genome_map = genome_per_contig.reindex(pivot_df.index).to_dict()
    else:
        contig_to_borg = profile_df.groupby('contig')['BORG_Ref'].first()
        contig_genome_map = {c: ('BORG' if (pd.notna(v) and v != 'NA') else 'HOST') 
                            for c, v in contig_to_borg.reindex(pivot_df.index).items()}
    
    # Build graph
    G = nx.Graph()
    for i, contig1 in enumerate(pivot_df.index):
        # Add node with genome type attribute
        G.add_node(contig1, Genome=contig_genome_map.get(contig1, 'NA'))
        for j, contig2 in enumerate(pivot_df.index):
            if i < j:
                similarity = similarity_matrix[i, j]
                if similarity > 0.7:
                    G.add_edge(contig1, contig2, weight=similarity)

                    ## if the edge is between none-Mp contig and others, print it out
                    if (contig_genome_map.get(contig1, 'NA') == 'Non-Mp' or 
                        contig_genome_map.get(contig2, 'NA') == 'Non-Mp'):
                        if not (contig_genome_map.get(contig1, 'NA') == 'Non-Mp' and 
                                contig_genome_map.get(contig2, 'NA') == 'Non-Mp'):
                            print(f"Edge between Non-Mp contig {contig1} and {contig2} with similarity {similarity:.2f}")
                            ## print their shared motifs
                            shared_motifs = pivot_binary.columns[(pivot_binary.loc[contig1] == 1) & (pivot_binary.loc[contig2] == 1)].tolist()
                            print(f"Shared motifs: {shared_motifs}")
    
    print(f"Network has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    # Get unique genome types and create color palette
    unique_genomes = list(set(contig_genome_map.values()))
    n_genomes = len(unique_genomes)
    if n_genomes <= 10:
        genome_colors = sns.color_palette("tab10", n_genomes)
    elif n_genomes <= 20:
        genome_colors = sns.color_palette("tab20", n_genomes)
    else:
        genome_colors = sns.color_palette("hls", n_genomes)
    genome_palette = dict(zip(unique_genomes, genome_colors))
    
    # Plot the network
    fig, ax = plt.subplots(figsize=(16, 14))
    pos = nx.spring_layout(G, seed=42, k=0.5, iterations=50)
    
    # Draw edges with width proportional to weight
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=[w*3 for w in weights], 
                          alpha=0.3, edge_color='gray', ax=ax)
    
    # Draw nodes colored by genome type
    for genome_type in unique_genomes:
        nodelist = [node for node in G.nodes() if G.nodes[node]['Genome'] == genome_type]
        nx.draw_networkx_nodes(G, pos, nodelist=nodelist, 
                              node_color=[genome_palette[genome_type]], 
                              node_size=300, alpha=0.9, 
                              edgecolors='black', linewidths=1, 
                              label=genome_type, ax=ax)
    
    # Draw labels with smaller font
    nx.draw_networkx_labels(G, pos, font_size=6, font_color='black', ax=ax)
    
    ax.set_title('Network of Genomes Based on Motif Profile Pearson Correlation (>0.3)', 
                fontsize=16, pad=20)
    ax.legend(title='Genome Type', fontsize=10, loc='upper left', frameon=True)
    ax.axis('off')
    
    plt.tight_layout()
    network_file = plot_name.replace('.pdf', '_network.pdf')
    plt.savefig(network_file, dpi=300, bbox_inches='tight')
    print(f"Saved network plot to: {network_file}")
    
    ## save the graph in gml format - convert all attributes to strings
    gml_file = plot_name.replace('.pdf', '_network.gml')
    G_copy = G.copy()
    # Convert all edge weights to strings
    for u, v in G_copy.edges():
        G_copy[u][v]['weight'] = str(G_copy[u][v]['weight'])
    # Convert all node attributes to strings
    for node in G_copy.nodes():
        for key, value in G_copy.nodes[node].items():
            G_copy.nodes[node][key] = str(value)
    nx.write_gml(G_copy, gml_file)
    print(f"Saved network to GML format: {gml_file}")
    
    plt.show()

def load_ece_classification():
    ece_dict = {}
    with open("ECE_anno.csv", 'r') as f:
        header = f.readline()
        for line in f:
            fields = line.strip().split(',')
            if len(fields) >= 2:
                ctg_name = fields[0]
                ece_type = fields[1]
                ece_dict[ctg_name] = ece_type
    return ece_dict

if __name__ == "__main__":
    
    all_dir = "/home/shuaiw/borg/paper/gg_run2/"
    seq_dir = "/home/shuaiw/borg/paper/borg_data/profile3/"
    cluster = "profile"

    cluster_species = "borg"
    borg_file = f"all_{cluster_species}_contigs_summary.tsv"
    os.system(f"cat {all_dir}/*/borg/{cluster_species}_contigs_summary.tsv > {borg_file}")
    plot_name = os.path.join(seq_dir, f"{cluster_species}_motif_profile_all.pdf")
    
    # """
    # Load BORG data
    borg_data = My_Borg(borg_file)
    high_dp_borgs = borg_data.get_high_depth_borgs(min_depth=5.0)
    
    members = []
    borg_anno_dict = {}
    borg_indicator = {}
    ece_dict = load_ece_classification()
    for i, entry in enumerate(high_dp_borgs):
        # print(f"{i+1}. {entry}")
        members.append((entry.seq_name, entry.sample_name))
        borg_anno_dict[entry.seq_name] = [entry.type, entry.borg_ref]
        if entry.borg_ref in ece_dict:
            borg_indicator[entry.seq_name] = ece_dict[entry.borg_ref]
        else:
            ## report error if not in ece_dict, and stop the program
            sys.stderr.write(f"Error: {entry.borg_ref} not found in ECE classification dictionary.\n")
            
        # print (entry.borg_ref)
    # members = ["soil_1_1336_L", "soil_s4_1_109_C"]
    
    
    fasta_dict = find_assembly()
    all_members, all_anno_dict, non_Mp_members, non_Mp_anno_dict = collect_host_genus(all_dir, fasta_dict)
    ## if member in all members not in members, add to members, also update borg_anno_dict
    for member in all_members:
        if member not in members:
            members.append(member)
            borg_anno_dict[member[0]] = ['HOST', all_anno_dict[member[0]][3:]+"<GTDB>"]
            borg_indicator[member[0]] = 'Mp'
    # # print (borg_anno_dict)
    # count_mod_freq(all_dir, borg_anno_dict)

    manual_members = [("SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META_16008_L", "soil_2")]
    for member in manual_members:
        if member not in members:
            members.append(member)
            borg_anno_dict[member[0]] = ['HOST', "<Manual host>"]
            borg_indicator[member[0]] = 'Mp'


    ### add non-Mp members as well, randomly select 100 members
    import random
    from sklearn.metrics.pairwise import cosine_similarity
    random.seed(42)
    selected_non_Mp = random.sample(non_Mp_members, min(75, len(non_Mp_members)))
    for member in selected_non_Mp:
        if member not in members:
            members.append(member)
            borg_anno_dict[member[0]] = ['Non-Mp', non_Mp_anno_dict[member[0]][3:]+"<GTDB>"]
            borg_indicator[member[0]] = 'Non-Mp'

    print (len(members), "members in total")
    for m in members:
        print (m, borg_anno_dict[m[0]], borg_indicator[m[0]])
    cluster_obj = given_species_drep_fast(all_dir, members, seq_dir, cluster,
                                    seq_dir, seq_dir, min_frac=0.6, 
                                    min_sites=20, score_cutoff = 30, max_len=100000)
    # cluster_obj.plot_profile(cluster, plot_name, cluster_species)

    # cluster_obj = My_cluster(cluster, members) 
    # cluster_obj.load_df(seq_dir)


    ## add a column to indicate borg_ref
    cluster_obj.profile_df['BORG_Ref'] = cluster_obj.profile_df['contig'].apply(lambda x: borg_anno_dict[x][1] if x in borg_anno_dict else 'NA')
    ## add a column to indicate borg/host
    cluster_obj.profile_df['Genome'] = cluster_obj.profile_df['contig'].apply(lambda x: borg_indicator[x] if x in borg_indicator else 'NA')

    
    cluster_obj.profile_df.to_csv(f"{seq_dir}/{cluster}_profile_df.tsv", index=False)
    profile_df = cluster_obj.profile_df

    remove_motifs = ['GATCH_4','BATC_2','RGAYCY_3','YGATCB_3','BGATATC_5',"GGHCCD_5","HGGCC_5","YCDVH_2","YCTAARAR_2",\
                     "YGATCBB_3","GGNCCH_5","GGCCH_4","DCCWGG_3","CCCTGH_3","GGNDCC_5","RGATCT_5","RGATCY_5","RGAYCB_3",\
                       "DYCACGRND_3","YCDV_2", "YDCCGGHR_3","VSAB_3","VNNNNNNNNNNNNNGCAYNNNNNNHTNGC_17","HNNNNCAGNNNNNNGTAG_7",\
                        "HNDNNNBGATCHNV_11","GGACCANNNNNNNNNNND_3","DNNNNNNNNNNSAGCTSNNNNNNNNNNHB_15","DNNNDYCTAADR_7","DNNGAGNNNNNNNTTTG_5",\
                            "DGAGNNNNNGGC_3","DGADNNNNNNTCGC_3","CTAGNNNNNBH_1"]
    # profile_df = pd.read_csv(f"{seq_dir}/{cluster}_profile_df.tsv")
    profile_df = profile_df[~profile_df['motifString'].isin(remove_motifs)]
    
    # Remove motifs where all contigs have fraction < 0.2
    motif_max_fractions = profile_df.groupby('motifString')['fraction'].max()
    motifs_to_keep = motif_max_fractions[motif_max_fractions >= 0.4].index

    ## only keep one for each reverse complementary motif pair
    unique_motif_ids = get_unique_motif(motifs_to_keep)
    print ("###get unique motifs", len(motifs_to_keep), len(unique_motif_ids))
    motifs_to_keep = unique_motif_ids

    profile_df = profile_df[profile_df['motifString'].isin(motifs_to_keep)]
    print(f"Kept {len(motifs_to_keep)} motifs with max fraction >= 0.2")
    
    # Remove contigs where all motifs have fraction < 0.2
    contig_max_fractions = profile_df.groupby('contig')['fraction'].max()
    contigs_to_keep = contig_max_fractions[contig_max_fractions >= 0.4].index
    profile_df = profile_df[profile_df['contig'].isin(contigs_to_keep)]
    print(f"Kept {len(contigs_to_keep)} contigs with max fraction >= 0.2")

    profile_df.to_csv(f"{seq_dir}/{cluster}_profile_df_filtered.csv", index=False)

    # """
    ## store the profile_df in a CSV file
    
    profile_df = pd.read_csv(f"{seq_dir}/{cluster}_profile_df_filtered.csv")
    umap_plot(profile_df, plot_name)
    # personal_plot(profile_df)
    # profile_df['sample'] = profile_df['contig'].apply(extract_sample_name)
    # ## get a df for each sample
    # for sample in profile_df['sample'].unique():
    #     sample_df = profile_df[profile_df['sample'] == sample]
    #     sample_plot_name = plot_name.replace('.pdf', f'_{sample}_heatmap.pdf')
    #     umap_plot(sample_df, sample_plot_name)