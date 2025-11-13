"""
data object for the sample
"""

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

def get_unique_motifs(df_motif, min_frac=0.4, min_sites = 100):
    df_motif = df_motif[(df_motif['fraction'] >= min_frac) & (df_motif['nDetected'] >= min_sites)]
    ## rm redundant motifs which are reverse complement 
    unique_motifs = []
    for index, row in df_motif.iterrows():
        if row['motifString'] not in unique_motifs and  str(Seq(row['motifString']).reverse_complement()) not in unique_motifs:
            unique_motifs.append(row['motifString'])
    return len(unique_motifs), unique_motifs

def get_unique_motifs_simple(df_motif):
    unique_motifs = []
    for index, row in df_motif.iterrows():
        if row['motifString'] not in unique_motifs and  str(Seq(row['motifString']).reverse_complement()) not in unique_motifs:
            unique_motifs.append(row['motifString'])
    return unique_motifs

def count_N50_size(fai):
    """
    Count N50 size from a fasta index file.
    """
    total_length = 0
    contig_lengths = []
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            contig_lengths.append(length)
            total_length += length
    contig_lengths.sort(reverse=True)
    
    n50_size = 0
    half_length = total_length / 2
    current_length = 0
    for length in contig_lengths:
        current_length += length
        if current_length >= half_length:
            n50_size = length
            break
    return n50_size, total_length

def get_best_ctg(depth_file, fai, min_len = 1000000, min_dp = 10):
    """
    Get the best contig based on length from a fasta file.
    """
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    length_dict = {}
    for index, row in depth_df.iterrows():

        if row['depth'] >= min_dp:
            good_depth[row['contig']] = row['depth']
    best_ctgs = {}
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            length_dict[ctg] = length
            if ctg[-1] == "C" and length >= min_len and ctg in good_depth:
                best_ctgs[ctg] = length
    # print (f"Total {len(best_ctgs)} contigs with length >= {min_len} found.")
    return good_depth, length_dict

def classify_taxa(lineage, level="species"):
    ## given a lineage string, return the taxa at the given level
    levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    if lineage == "Unknown":
        return "Unknown"
    lineage = lineage.split(";")
    if len(lineage) == 0:
        return "Unknown"
    if len(lineage) < 7:
        return "Unknown"
    taxon = lineage[levels.index(level)]
    if len(taxon) == 3:
        return "Unknown"
    return taxon

def get_ctg_taxa(all_dir):
    ctg_taxa_dict = {}
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        sample_obj = My_sample(prefix, all_dir)
        sample_taxa_dict = sample_obj.read_meta_gtdb()
        ctg_taxa_dict.update(sample_taxa_dict)
    print (len(ctg_taxa_dict), "contig taxa info collected")
    return ctg_taxa_dict

class My_sample(object):
    def __init__(self, prefix, all_dir):
        self.prefix = prefix
        self.all_dir = all_dir

        self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation3"
        self.reference_fasta = f"{self.all_dir}/{self.prefix}/{self.prefix}.hifiasm.p_ctg.rename.fa"
        self.fai = f"{self.all_dir}/{self.prefix}/{self.prefix}.hifiasm.p_ctg.rename.fa.fai"
        self.map_sum = f"{self.all_dir}/{self.prefix}/{self.prefix}.align.count.csv"
        self.all_host_file = f"{self.all_dir}/{self.prefix}/all_host_ctgs.tsv"
        self.depth_file = os.path.join(self.work_dir, "mean_depth.csv")
        self.host_sum_file = os.path.join(self.work_dir, "host_summary.csv")
        self.orphan_file = os.path.join(self.work_dir, "regulatory_motif_enrichment.csv")
        self.motif_freq_file = os.path.join(self.work_dir, "motif_length_stats.csv")
        self.profile = os.path.join(self.work_dir, "motif_profile.csv")
        self.gtdb = os.path.join(self.work_dir, "../GTDB/gtdbtk.bac120.summary.tsv")
        self.checkm = os.path.join(self.work_dir, "../checkM2/quality_report.tsv")
        self.mge_file = f"{self.all_dir}/{self.prefix}/all_mge.tsv"
        self.all_motif_file = f"{self.work_dir}/all.motifs.csv"
        self.bin3c_cluster = f"{all_dir}/{prefix}/hic/bin3c_clust/clustering.mcl"
        self.contact_value_file = f"{all_dir}/{prefix}/hic/bin3c/contact_values.txt"
        self.spacer_linkage_file = f"{all_dir}/{prefix}/spacer/{prefix}_mge_spacer_hits.filter.tsv"
        
        self.mge_dict = None
        self.depth_dict = None
        self.length_dict = None
        self.depth_cutoff = 10
        self.length_cutoff = 5000
        self.specificity_cutoff = 0.01
        self.final_score_cutoff = 0.5

    def get_unique_motifs(self):
        if not os.path.exists(self.all_motif_file):
            print (f"[⚠️] Motif file not found: {self.all_motif_file}")
            return None, None
        df = pd.read_csv(self.all_motif_file)
        return get_unique_motifs(df)

    def read_MGE(self):
        if not os.path.exists(self.mge_file):
            print (f"[⚠️] MGE file not found: {self.mge_file}")
            return None
        self.mge_dict = {}
        f = open(self.mge_file, "r")
        for line in f:
            if line.startswith("seq_name"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            contig = fields[0]
            mge_type = fields[1]
            self.mge_dict[contig] = mge_type
        return self.mge_dict

    def get_MGE_bool(self):
        self.read_MGE()
        if self.mge_dict is None:
            return None
        else:
            if len(self.mge_dict) == 0:
                return 0
            else:
                return 1

    ## read length of each contig
    def read_fai(self):
        contig_len = {}
        if not os.path.exists(self.fai):
            print (f"[⚠️] FAI file not found: {self.fai}")
            return contig_len
        with open(self.fai, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 2:
                    continue
                contig = fields[0]
                length = int(fields[1])
                contig_len[contig] = length
        return contig_len

    # def read_host(self):
    #     line_num = sum(1 for line in open(self.host_sum_file) if line.strip())
    #     if line_num > 1:
    #         host_sum = pd.read_csv(self.host_sum_file)
    #         all_link_df = pd.concat([all_link_df, host_sum], axis=0)

    def read_host(self):
        # return 0
        print (f"Reading host summary file: {self.host_sum_file}")
        line_num = sum(1 for line in open(self.host_sum_file) if line.strip())
        if line_num < 2:
            print(f"Host summary file {self.host_sum_file} is empty or has only header.")
            return 0
        # Read the host summary file
        host_sum = pd.read_csv(self.host_sum_file)
        linkage_num = 0
        for index, row in host_sum.iterrows():
            # if row["pvalue"] >= 0.05:
            #     continue
            # if row['final_score'] <= 0.6:
            #     continue
            if row["specificity"] >= self.specificity_cutoff:
                continue
            if row['final_score'] <= self.final_score_cutoff:
                continue
            linkage_num += 1
        return linkage_num

    def read_linkage_dict(self, ctg2bin_dict={}):
        ## if self.host_sum_file not exists
        if not os.path.exists(self.host_sum_file):
            print(f"Host summary file {self.host_sum_file} not found.")
            return {}, {}
        line_num = sum(1 for line in open(self.host_sum_file) if line.strip())
        if line_num < 2:
            print(f"Host summary file {self.host_sum_file} is empty or has only header.")
            return {}, {}

        df = pd.read_csv(self.host_sum_file)
        df = df[df['specificity'] < self.specificity_cutoff]
        df = df[df['final_score'] > self.final_score_cutoff]
        our_linkages = defaultdict(list)
        our_ctg_linkages = {}
        for index, row in df.iterrows():
            if row['host'] not in ctg2bin_dict:
                bin_name = row['host']
                ## raise error
                # print (f"contig {row['host']} is not in ctg2bin_dict")
                # sys.exit(1)
            else:
                bin_name = ctg2bin_dict[row['host']]
            ## try plasmid is in header of row
            try:
                plasmid_name = row['plasmid']
            ## otherwise index with MGE
            except KeyError:
                plasmid_name = row['MGE']
            our_linkages[plasmid_name].append(bin_name)
            our_ctg_linkages[plasmid_name] = row['host']
            # our_linkages[row['plasmid']] = row['host']
        multiple_host_plasmid_num = 0
        for plasmid in our_linkages:
            if len(our_linkages[plasmid]) > 1:
                multiple_host_plasmid_num += 1
                # print (f"{plasmid} has multiple host: {our_linkages[plasmid]}")
        print (f"multiple host plasmid num: {multiple_host_plasmid_num} out of {len(our_linkages)}")
        return our_linkages, our_ctg_linkages

    def read_spacer(self, mismatch_allowed=0):
        spacer_linkage_dict = defaultdict(set)
        ## check if spacer_linkage_file empty
        if not os.path.exists(self.spacer_linkage_file) or os.path.getsize(self.spacer_linkage_file) == 0:
            return spacer_linkage_dict
        df = pd.read_csv(self.spacer_linkage_file, sep="\t")
        # print (df)
        for index, row in df.iterrows():
            if row['full_mismatch'] > mismatch_allowed:
                continue
            spacer_linkage_dict[row['target_id']].add(row['query_contig_id'])
        return spacer_linkage_dict

    def read_bin3c(self):
        bin3c_cluster_list = []
        with open(self.bin3c_cluster, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                bin3c_cluster_list.append(line.split())
        return bin3c_cluster_list


    def read_depth(self):
        self.depth_dict = {}
        self.length_dict = {}

        if not os.path.exists(self.depth_file):
            print(f"[⚠️] Depth file not found: {self.depth_file}")
            return self.depth_dict, self.length_dict
        
        try:
            # Read the CSV file
            df = pd.read_csv(self.depth_file)
            
            # Check if required columns exist
            required_cols = ['contig', 'depth']
            if 'contig' not in df.columns or 'depth' not in df.columns:
                print(f"[⚠️] Depth file missing required columns. Expected: contig, depth, Found: {list(df.columns)}")
                return self.depth_dict, self.length_dict
            
            # Check if length column exists, if not, we'll set length to 0
            has_length_col = 'length' in df.columns
            if not has_length_col:
                print(f"[⚠️] Length column not found in depth file, setting all lengths to 0")
            
            # Create dictionaries
            for _, row in df.iterrows():
                contig = row['contig']
                depth = float(row['depth'])
                if has_length_col:
                    length = int(row['length'])
                else:
                    length = 0

                self.depth_dict[contig] = depth
                self.length_dict[contig] = length

            # print(f"[✔] Successfully read {len(self.depth_dict)} contigs from depth file")
            # print(f"    Depth range: {min(self.depth_dict.values()):.2f} - {max(self.depth_dict.values()):.2f}")
            # print(f"    Length range: {min(self.length_dict.values())} - {max(self.length_dict.values())}")
            
        except Exception as e:
            print(f"[⚠️] Error reading depth file {self.depth_file}: {str(e)}")
            return {}, {}
        
        return self.depth_dict, self.length_dict

    def read_mapping(self):
        """
        map_read_count,raw_read_count,ratio
        985092,1043059,0.9444259624815087
        """
        if not os.path.exists(self.map_sum):
            return 0
        ## get the ratio
        f = open(self.map_sum, "r")
        for line in f:
            if line.startswith("#"):
                continue
            map_read_count, raw_read_count, ratio = line.strip().split(",")
            print(f"Map read count: {map_read_count}, Raw read count: {raw_read_count}, Ratio: {ratio}")
        f.close()
        return float(ratio)

    def read_orphan(self):
        regulate_motif_num = 0
        if not os.path.exists(self.orphan_file):
            return regulate_motif_num
        df = pd.read_csv(self.orphan_file)
        regulate_motif_set = set()
        for index, row in df.iterrows():
            if row["significant"] != True:
                continue
            motif_tag = f"{row['motif']}_{row['exact_pos']}"
            regulate_motif_set.add(motif_tag)
        regulate_motif_num = len(regulate_motif_set)
        return regulate_motif_num

    def get_N50_size(self):
        N50, genome_size = count_N50_size(self.fai)
        return N50, genome_size

    def get_final_best_ctg(self):
        good_depth, self.length_dict = get_best_ctg(self.depth_file, self.fai)
        """
        Read all host contigs from a file.
        """
        best_ctgs = {}
        with open(self.all_host_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                ctg, ctg, domain = line.strip().split("\t")
                if ctg not in good_depth:
                    continue
                best_ctgs[ctg] = domain
        return best_ctgs

    def get_high_dp_ctg_list(self, min_depth=10, min_len=100000):
        """
        Get high depth contig list.
        """
        genome_list = []
        for contig in self.depth_dict:
            if self.depth_dict[contig] >= min_depth and self.length_dict[contig] >= min_len:
                genome = self.work_dir + "/contigs/" + contig + ".fa"
                genome_list.append(genome)
        return genome_list

    def read_meta_gtdb(self):
        """
        Read the GTDB summary file and return a dictionary of contig to bin mapping.
        """
        if not os.path.exists(self.gtdb):
            print (f"[⚠️] GTDB file not found: {self.gtdb}")
            return {}
        gtdb_df = pd.read_csv(self.gtdb, sep='\t')
        isolation_taxa = {}
        for index, row in gtdb_df.iterrows():
            anno = row['classification']
            isolation_taxa[row['user_genome']] = anno
        return isolation_taxa


class Isolation_sample(My_sample):

    def __init__(self, prefix, all_dir):
        super().__init__(prefix, all_dir)
        self.phylum = None
        self.lineages = None
        self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation2"
        self.depth_file = os.path.join(self.work_dir, "mean_depth.csv")
        self.host_sum_file = os.path.join(self.work_dir, "host_summary.csv")
        self.orphan_file = os.path.join(self.work_dir, "regulatory_motif_enrichment.csv")
        self.motif_freq_file = os.path.join(self.work_dir, "motif_length_stats.csv")
        self.profile = os.path.join(self.work_dir, "motif_profile.csv")
        self.gtdb = os.path.join(self.work_dir, "../GTDB_2/gtdbtk.bac120.summary.tsv")
        self.all_gtdb = os.path.join(self.work_dir, "../GTDB/gtdbtk.bac120.summary.tsv")
        self.all_motif_file = f"{self.work_dir}/all.motifs.csv"

    def get_phylum(self):
        isolation_taxa = self.read_isolation_gtdb()
        self.lineage = isolation_taxa[list(isolation_taxa.keys())[0]]
        self.phylum = self.lineage.split(";")[1][3:] if ";" in self.lineage else "Unclassified"

    def read_isolation_gtdb(self):
        """
        Read the GTDB summary file and return a dictionary of contig to bin mapping.
        """
        gtdb_df = pd.read_csv(self.all_gtdb, sep='\t')
        isolation_taxa = {}
        for index, row in gtdb_df.iterrows():
            anno = row['classification']
            sra_id = row['user_genome']
            isolation_taxa[sra_id] = anno
            return isolation_taxa

    def check_pure2(self):
        ## read checkm file
        if not os.path.exists(self.checkm):
            print (f"[⚠️] CheckM file not found: {self.checkm}")
            return "unknown"
        df = pd.read_csv(self.checkm, sep="\t")
        if "Completeness" not in df.columns or "Contamination" not in df.columns:
            print (f"[⚠️] CheckM file does not have required columns.")
            return "unknown"
        completeness = df["Completeness"].values[0]
        contamination = df["Contamination"].values[0]
        if contamination <= 5:
            return "pure"
        else:
            return "mixed"

    def get_mge_specific_motif(self, MGE_ctg, host_ctg, min_frac=0.8, min_sites=20):
        MGE_motif_file = os.path.join(self.work_dir, "motifs/" + MGE_ctg + ".motifs.csv")
        df = pd.read_csv(MGE_motif_file)
        df = df[df["fraction"] >= min_frac]
        df = df[df["nDetected"] >= min_sites]

        host_df = pd.read_csv(os.path.join(self.work_dir, "motifs/" + host_ctg + ".motifs.csv"))
        host_motifs = set(host_df["motifString"].tolist())
        mge_specific_motifs = []
        for index, row in df.iterrows():
            if row["motifString"] not in host_motifs and str(Seq(row["motifString"]).reverse_complement()) not in host_motifs:
                mge_specific_motifs.append(row["motifString"])
        return mge_specific_motifs

    def explore_specific_motifs(self):    
        host_ctgs = []
        MGE_ctgs = []
        self.read_MGE()
        for contig in self.depth_dict:
            if self.depth_dict[contig] < self.depth_cutoff: continue
            if self.length_dict[contig] < self.length_cutoff: continue
            if contig not in self.mge_dict:
                host_ctgs.append(contig)
            else:
                MGE_ctgs.append(contig)
        for mge_ctg in MGE_ctgs:
            for host_ctg in host_ctgs:
                mge_specific_motifs = self.get_mge_specific_motif(mge_ctg, host_ctg)
                if len(mge_specific_motifs) > 0:
                    print (f"{self.prefix}\t{mge_ctg}\t{host_ctg}\t{len(mge_specific_motifs)}\t{mge_specific_motifs}")


class My_contig(My_sample):

    def __init__(self, prefix, all_dir, contig):
        super().__init__(prefix, all_dir)
        self.contig = contig
        self.ctg_ref = f"{all_dir}/{prefix}/{prefix}_methylation3/contigs/{contig}.fa"
        self.gff = f"{all_dir}/{prefix}/{prefix}_methylation3/gffs/{contig}.gff"
        self.ipd_ratio_file = f"{all_dir}/{prefix}/{prefix}_methylation3/ipd_ratio/{contig}.ipd3.csv"
        self.motif_file = f"{all_dir}/{prefix}/{prefix}_methylation3/motifs/{contig}.motifs.csv"

    def read_motif(self, min_frac=0.3, min_sites=30):
        ## check if file exists
        if not os.path.exists(self.motif_file):
            return None
        motif_df = pd.read_csv(self.motif_file)
        motif_df = motif_df[(motif_df['fraction'] >= min_frac) & (motif_df['nDetected'] >= min_sites)]
        return motif_df



class My_cluster(object):
    def __init__(self, cluster, members):
        self.cluster = cluster
        self.members = members  
        self.profile_df = None  # ["contig", "motifString", "fraction"]
        self.filtered_motif_file = "../motif_change/manual_filtered_motifs.csv"
    
    def get_profile(self, df, tmp_res_file):
        self.profile_df = df
        self.profile_df.to_csv(tmp_res_file)

    def plot_profile(self, closed_genome, plot_name, cluster_species):
        
        duplicates = self.profile_df.duplicated(subset=['contig', 'motifString'], keep=False)
        if duplicates.any():
            print(f"Warning: Found {duplicates.sum()} duplicate entries. Aggregating by taking mean.")
            self.profile_df = self.profile_df.groupby(['contig', 'motifString'])['fraction'].mean().reset_index()
        try:
            df_pivot = self.profile_df.pivot(index='contig', columns='motifString', values='fraction').fillna(0)
        except ValueError as e:
            print(f"Error creating pivot table: {e}")
            print("DataFrame info:")
            print(f"Shape: {self.profile_df.shape}")
            print(f"Columns: {self.profile_df.columns.tolist()}")
            print(f"Sample data:\n{self.profile_df.head()}")
            return
        
        if df_pivot.shape[0] < 2 or df_pivot.shape[1] < 2:
            print(f"Warning: Insufficient data for clustering (shape: {df_pivot.shape}). Plotting without clustering.")
            df_clustered = df_pivot
        else:
            try:
                # Cluster rows and columns
                row_linkage = linkage(df_pivot.values, method='average')
                col_linkage = linkage(df_pivot.values.T, method='average')
                row_order = leaves_list(row_linkage)
                col_order = leaves_list(col_linkage)
                df_clustered = df_pivot.iloc[row_order, col_order]
            except ValueError as e:
                print(f"Warning: Clustering failed ({e}). Plotting without clustering.")
                df_clustered = df_pivot
        ## check if df_clustered is empty
        if df_clustered.empty:
            print("No data available to plot.")
            return
        plt.figure(figsize=(10, 6))
        ax = sns.heatmap(df_clustered, cmap="YlGnBu", annot=False, fmt=".2f", cbar_kws={'label': 'Fraction'})
        ax.set_xlabel("Motif String")
        ax.set_ylabel("Contig")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ## set title as closed_genome
        ax.set_title(f"{cluster_species} - {closed_genome}")
        plt.tight_layout()
        plt.savefig(plot_name)
        
    
    def check_diff_motifs(self):
        ## if there is a contig has fraction < 0.1, regard lack a motif, so different motifs exist
        minus_frac = 0.1
        lack_motif_df = self.profile_df[self.profile_df['fraction'] < minus_frac]
        if not lack_motif_df.empty:
            return "variation"
        else:
            return "uniform"
        
    def pairwise_compare(self, bin_freq = 0.3):
        ## for every pair of contigs, compare their motif profiles using cosine similarity
        ## only retain unique motifs using get_unique_motifs
        unique_motifs = get_unique_motifs_simple(self.profile_df)
        ## filter the profile_df to only include these motifs
        filtered_df = self.profile_df[self.profile_df['motifString'].isin(unique_motifs)]

        df_pivot = filtered_df.pivot(index='contig', columns='motifString', values='fraction').fillna(0)
        contig_list = df_pivot.index.tolist()
        profile_matrix = df_pivot.values
        if len(df_pivot) < 2:
            print("Warning: Insufficient data for pairwise comparison.")
            return [], 0
        similarity_matrix = cosine_similarity(profile_matrix)
        ## also count Jaccard similarity, binary presence/absence with cutoff 0.3
        binary_df = df_pivot.applymap(lambda x: 1 if x >= bin_freq else 0)
        binary_matrix = binary_df.values

        n = len(contig_list)
        data = []
        for i in range(n):
            for j in range(i+1, n):
                jaccard = jaccard_score(binary_matrix[i], binary_matrix[j])
                data.append([similarity_matrix[i][j], jaccard])
        return data, self.uniq_clade_num(binary_matrix)

    def uniq_clade_num(self, binary_matrix):
        ### cluster the binary matrix rows and count unique clades, 
        # in each clade, the jaccard similarity between any two contigs is 1
        # use hierarchical clustering
        from scipy.cluster.hierarchy import fcluster
        linkage_matrix = linkage(binary_matrix, method='ward')
        # get the cluster labels with t=0 for perfect similarity (Jaccard=1)
        cluster_labels = fcluster(linkage_matrix, t=0.2, criterion='distance')
        return len(set(cluster_labels))


    def load_df(self, tmp_res):
        tmp_res_file = f"{tmp_res}/{self.cluster}.csv"
        if not os.path.exists(tmp_res_file):
            print(f"[⚠️] Profile file not found: {tmp_res_file}")
            return
        self.profile_df = pd.read_csv(tmp_res_file)

    def manual_filter_motifs(self):
        if not os.path.exists(self.filtered_motif_file):
            print(f"[⚠️] Filtered motif file not found: {self.filtered_motif_file}")
            return
        df_filter = pd.read_csv(self.filtered_motif_file)
        self.filtered_motifs = df_filter['motifString'].tolist()
        # filter self.profile_df to exclude these motifs
        self.profile_df = self.profile_df[~self.profile_df['motifString'].isin(self.filtered_motifs)]
       