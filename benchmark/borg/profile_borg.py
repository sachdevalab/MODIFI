import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
import re
from Bio.Seq import Seq

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
from sample_object import get_unique_motifs, My_sample, get_ctg_taxa, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
from check_motif_change import given_species_drep
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
                print (ctg, ctg_lineage)
                # os.system(f"cp {ctg_obj.ctg_ref} /home/shuaiw/borg/paper/borg_data/align/refs")
                members.append((ctg, prefix))
                anno_dict[ctg] = taxon
    print (f"Total host contigs collected: {len(members)}")
    return members, anno_dict

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

def umap_plot(profile_df, plot_name):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import umap

    # Add sample column extracted from contig name (first 8 underscore-separated parts)
    def extract_sample_name(contig_name):
        parts = str(contig_name).split('_')
        if len(parts) >= 8:
            return '_'.join(parts[:8])
        else:
            return contig_name

    profile_df['sample'] = profile_df['contig'].apply(extract_sample_name)

    # Create pivot table for UMAP
    pivot_df = profile_df.pivot(index='contig', columns='motifString', values='fraction')
    pivot_df = pivot_df.fillna(0)  # Fill NaN values with 0

    # Perform UMAP dimensionality reduction
    reducer = umap.UMAP(n_neighbors=5, min_dist=0.3, metric='correlation', random_state=42)
    embedding = reducer.fit_transform(pivot_df)

    # Create a DataFrame for the embedding
    umap_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'], index=pivot_df.index)

    # Determine genome type per contig for marker shapes
    if 'Genome' in profile_df.columns:
        # derive one genome label per contig (first occurrence) and reindex to pivot_df order
        genome_per_contig = profile_df.groupby('contig')['Genome'].first()
        genome_series = genome_per_contig.reindex(pivot_df.index).fillna('NA')
    else:
        # fallback: infer from BORG_Ref presence
        contig_to_borg = profile_df.groupby('contig')['BORG_Ref'].first()
        genome_series = contig_to_borg.reindex(pivot_df.index).apply(lambda v: 'BORG' if (pd.notna(v) and v != 'NA') else 'HOST')

    umap_df['Genome'] = genome_series.values

    # Dynamically assign marker shapes to each unique genome type
    unique_genomes = umap_df['Genome'].unique().tolist()
    marker_shapes = ['o', 's', '^', 'D', 'P', 'X', '*', 'v', '<', '>', 'h', 'H', '8', 'p']  # extend as needed
    markers_map = {genome: marker_shapes[i % len(marker_shapes)] for i, genome in enumerate(unique_genomes)}

    # Color points by sample
    if 'sample' not in profile_df.columns:
        # ensure sample column exists (fallback to contig-based extraction)
        def extract_sample_name(contig_name):
            parts = str(contig_name).split('_')
            return '_'.join(parts[:8]) if len(parts) >= 8 else contig_name
        profile_df['sample'] = profile_df['contig'].apply(extract_sample_name)
        # re-create pivot and umap_df index to align (should already match)
        # note: umap_df index is pivot_df.index which matches profile_df contigs used above

    sample_names = list(umap_df.index.map(lambda c: profile_df.loc[profile_df['contig'] == c, 'sample'].iloc[0]))
    umap_df['sample'] = sample_names
    unique_samples = umap_df['sample'].unique().tolist()
    n_samples = len(unique_samples)
    import seaborn as sns
    if n_samples <= 20:
        colors = sns.color_palette("tab20", n_samples)
    else:
        colors = sns.color_palette("hls", n_samples)
    sample_palette = dict(zip(unique_samples, colors))

    # Plot UMAP using shapes to distinguish genome type and color by sample
    plt.figure(figsize=(15, 8))
    ax = sns.scatterplot(
        data=umap_df,
        x='UMAP1', y='UMAP2',
        hue='sample',
        style='Genome',
        markers=markers_map,
        palette=sample_palette,
        s=60,
        edgecolor='black',
        linewidth=0.3,
        alpha=0.6,
        legend='brief'
    )
    plt.title('UMAP Projection of Motif Profiles', fontsize=16)
    plt.xlabel('UMAP1', fontsize=14)
    plt.ylabel('UMAP2', fontsize=14)
    plt.grid(True)
    # Move legend to the right outside the plot
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='Sample', fontsize=8)
    # Leave space on the right for the legend
    plt.tight_layout(rect=[0, 0, 0.83, 1])
    plt.savefig(plot_name.replace('.pdf', '_umap.pdf'), dpi=300, bbox_inches='tight')
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
    seq_dir = "/home/shuaiw/borg/paper/borg_data/profile/"
    cluster = "profile"

    cluster_species = "borg"
    borg_file = f"all_{cluster_species}_contigs_summary.tsv"
    os.system(f"cat {all_dir}/*/borg/{cluster_species}_contigs_summary.tsv > {borg_file}")
    plot_name = os.path.join(seq_dir, f"{cluster_species}_motif_profile_all.pdf")
    
    
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
    
    # """
    fasta_dict = find_assembly()
    all_members, all_anno_dict = collect_host_genus(all_dir, fasta_dict)
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


    # cluster_obj = given_species_drep(all_dir, members, seq_dir, cluster,
    #                                 seq_dir, seq_dir, min_frac=0.3, 
    #                                 min_sites=50, score_cutoff = 20)
    # cluster_obj.plot_profile(cluster, plot_name, cluster_species)

    cluster_obj = My_cluster(cluster, members) 
    cluster_obj.load_df(seq_dir)

    ## remove all rows with motifstring contains GATCH_4

    ## add a column to indicate borg_ref
    cluster_obj.profile_df['BORG_Ref'] = cluster_obj.profile_df['contig'].apply(lambda x: borg_anno_dict[x][1] if x in borg_anno_dict else 'NA')
    ## add a column to indicate borg/host
    cluster_obj.profile_df['Genome'] = cluster_obj.profile_df['contig'].apply(lambda x: borg_indicator[x] if x in borg_indicator else 'NA')

    
    cluster_obj.profile_df.to_csv(f"{seq_dir}/{cluster}_profile_df.tsv", index=False)
    profile_df = cluster_obj.profile_df

    remove_motifs = ['GATCH_4','BATC_2','RGAYCY_3','YGATCB_3','BGATATC_5',"GGHCCD_5","HGGCC_5","YCDVH_2","YCTAARAR_2",\
                     "YGATCBB_3","GGNCCH_5","GGCCH_4","DCCWGG_3","CCCTGH_3","GGNDCC_5","RGATCT_5","RGATCY_5","RGAYCB_3",\
                       "DYCACGRND_3","YCDV_2", "YDCCGGHR_3","VSAB_3",""]
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
    # personal_plot(profile_df)
    umap_plot(profile_df, plot_name)






