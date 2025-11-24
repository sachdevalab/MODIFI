import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
from sample_object import get_unique_motifs, My_sample, get_ctg_taxa, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
from check_motif_change import given_species_drep

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
            self.ctg_depth = float(fields[8])
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
                        borg_entry = Borg_Entry(line)
                        self.borg_entries.append(borg_entry)
                    except ValueError as e:
                        print(f"Warning: {e}")

    def get_high_depth_borgs(self, min_depth=10.0):
        """Get BORG entries with contig depth above threshold"""
        return [entry for entry in self.borg_entries if entry.ctg_depth >= min_depth]

def personal_plot(cluster_obj):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    from matplotlib.patches import Patch
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist
    # Add color bar for BORG/Host indication
    from matplotlib.patches import Rectangle
    
    # Create pivot table for heatmap
    pivot_df = cluster_obj.profile_df.pivot(index='contig', columns='motifString', values='fraction')
    pivot_df = pivot_df.fillna(0)  # Fill NaN values with 0

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
        contig_data = cluster_obj.profile_df[cluster_obj.profile_df['contig'] == contig].iloc[0]

        borg_ref = contig_data['BORG_Ref']

        label = f"{contig};({borg_ref})"

        y_labels.append(label)
    
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # Create heatmap
    sns.heatmap(pivot_df, cmap="YlGnBu", annot=False, fmt=".2f", 
                cbar_kws={'label': 'Fraction'}, ax=ax)
    
    
    # Set labels and formatting
    ax.set_xlabel("Motif String", fontsize=12)
    ax.set_ylabel("Contig", fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(y_labels, rotation=0, fontsize=8)
    ax.set_title(f"{cluster_species}", fontsize=14, pad=20)
    

    

    plt.tight_layout()
    plt.savefig(plot_name, dpi=300, bbox_inches='tight')
    plt.show()

def collect_host_genus(all_dir):
    members  = []
    anno_dict = {}
    ctg_taxa_dict = get_ctg_taxa(all_dir)
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if not re.search("soil", prefix):
            continue
        print (f"Processing {prefix}...")
        sample_obj = My_sample(prefix, all_dir)
        sample_obj.read_depth()
        sample_obj.get_len_dict()
        genome_list, contig_list =  sample_obj.get_high_dp_ctg_list( min_depth=5, min_len=100000)
        
        for ctg in contig_list:
            ctg_obj = My_contig(prefix, all_dir, ctg)
            ctg_lineage = ctg_taxa_dict[ctg] if ctg in ctg_taxa_dict else "Unknown"
            ctg_phylum = classify_taxa(ctg_lineage, level='phylum')
            ctg_genus = classify_taxa(ctg_lineage, level='genus')
            ctg_len = ctg_obj.get_ctg_len()
            taxon = get_detail_taxa_name(ctg_lineage)
            # print (ctg_genus)
            if ctg_genus == "g__Methanoperedens":
                print (ctg, ctg_lineage, ctg_len)
                members.append(ctg)
                anno_dict[ctg] = taxon
    print (f"Total host contigs collected: {len(members)}")
    return members, anno_dict


if __name__ == "__main__":
    borg_file = "all_borg_contigs_summary.tsv"
    all_dir = "/home/shuaiw/borg/paper/run2/"

    
    
    # Load BORG data
    borg_data = My_Borg(borg_file)
    high_dp_borgs = borg_data.get_high_depth_borgs(min_depth=5.0)
    
    members = []
    borg_anno_dict = {}
    for i, entry in enumerate(high_dp_borgs):
        print(f"{i+1}. {entry}")
        members.append(entry.seq_name)
        borg_anno_dict[entry.seq_name] = [entry.type, entry.borg_ref]
    print (members)
    # members = ["soil_1_1336_L", "soil_s4_1_109_C"]

    all_members, all_anno_dict = collect_host_genus(all_dir)
    ## if member in all members not in members, add to members, also update borg_anno_dict
    for member in all_members:
        if member not in members:
            members.append(member)
            borg_anno_dict[member] = ['HOST', all_anno_dict[member]]
    print (borg_anno_dict)
    seq_dir = "/home/shuaiw/borg/paper/borg_data/profile/"
    cluster = "profile"
    plot_name = os.path.join(seq_dir, f"borg_motif_profile.pdf")
    cluster_species = "borg & hosts"

    # cluster_obj = given_species_drep(all_dir, members, seq_dir, cluster,
    #                                 seq_dir, seq_dir, min_frac=0.3, 
    #                                 min_sites=10, score_cutoff = 30)
    # cluster_obj.plot_profile(cluster, plot_name, cluster_species)

    cluster_obj = My_cluster(cluster, members) 
    cluster_obj.load_df(seq_dir)

    ## remove all rows with motifstring contains GATCH_4
    cluster_obj.profile_df = cluster_obj.profile_df[cluster_obj.profile_df['motifString'] != 'GATCH_4']

    ## add a column to indicate borg_ref
    cluster_obj.profile_df['BORG_Ref'] = cluster_obj.profile_df['contig'].apply(lambda x: borg_anno_dict[x][1] if x in borg_anno_dict else 'NA')
    print (cluster_obj.profile_df)
    personal_plot(cluster_obj)




