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
import pickle

def get_unique_motifs(df_motif, min_frac=0.4, min_sites = 100):
    df_motif = df_motif[(df_motif['fraction'] >= min_frac) & (df_motif['nDetected'] >= min_sites)]
    ## rm redundant motifs which are reverse complement 
    unique_motifs = []
    unique_motifs_identifier = []
    for index, row in df_motif.iterrows():
        if row['motifString'] not in unique_motifs and  str(Seq(row['motifString']).reverse_complement()) not in unique_motifs:
            unique_motifs.append(row['motifString'])
            unique_motifs_identifier.append(row['motifString'] + "_" + str(row['centerPos']))
    return len(unique_motifs), unique_motifs, unique_motifs_identifier

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

def classify_taxa(lineage, level="species"):
    if level == "domain":
        if re.search(r'Unclassified Bacteria', lineage):
            return "d__Bacteria"
        if re.search(r'Unclassified Archaea', lineage):
            return "d__Archaea"
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

def get_detail_taxa_name(lineage):
    if lineage == "Unknown":
        detail_name = "Unknown"
    else:
        lineage = lineage.split(";")
        for level in reversed(lineage):
            if len(level.strip()) > 3:
                detail_name = level.strip()
                break
        else:
            detail_name = lineage[0].strip()
    return detail_name


class mge_obj:

    def __init__(self, name, type, length, method):
        self.name = name
        self.type = type
        self.length = length
        self.methods = [method]
        self.classification = None
    
    def add_method(self, method):
        if method not in self.methods:
            self.methods.append(method)

def get_ctg_taxa(all_dir, data_type="meta"):
    file_name = f"/home/shuaiw/borg/paper/gene_anno/{data_type}_ctg_taxa_dict.pkl"
    # if os.path.exists(file_name):
    #     with open(file_name, "rb") as f:
    #         return pickle.load(f)
    ctg_taxa_dict = {}
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if data_type == "meta":
            sample_obj = My_sample(prefix, all_dir)
        else:
            sample_obj = Isolation_sample(prefix, all_dir)
        sample_taxa_dict = sample_obj.read_meta_gtdb()
        ctg_taxa_dict.update(sample_taxa_dict)
    print (len(ctg_taxa_dict), "contig taxa info collected")
    ## save ctg_taxa_dict to a file
    
    with open(file_name, "wb") as f:
        pickle.dump(ctg_taxa_dict, f)
    return ctg_taxa_dict

class My_sample(object):
    def __init__(self, prefix, all_dir, data_type="meta", methy_v=4):
        self.prefix = prefix
        self.all_dir = all_dir

        if data_type == "meta":
            if methy_v == 3:
                self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation3"
            elif methy_v == 4:
                self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation4"
        else:
            self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation2"

        self.reference_fasta = f"{self.all_dir}/{self.prefix}/{self.prefix}.hifiasm.p_ctg.rename.fa"
        self.fai = f"{self.reference_fasta}.fai"
        self.map_sum = f"{self.all_dir}/{self.prefix}/{self.prefix}.align.count.csv"
        self.all_host_file = f"{self.all_dir}/{self.prefix}/all_host_ctgs.tsv"
        self.depth_file = os.path.join(self.work_dir, "mean_depth.csv")
        self.host_sum_file = os.path.join(self.work_dir, "host_summary.csv")
        self.orphan_file = os.path.join(self.work_dir, "regulatory_motif_enrichment.csv")
        self.motif_freq_file = os.path.join(self.work_dir, "motif_length_stats.csv")
        self.profile = os.path.join(self.work_dir, "motif_profile.csv")
        self.gtdb = os.path.join(self.work_dir, "../GTDB/gtdbtk.bac120.summary.tsv")
        self.arc_gtdb = os.path.join(self.work_dir, "../GTDB/gtdbtk.ar53.summary.tsv")
        self.checkm = os.path.join(self.work_dir, "../checkM2/quality_report.tsv")
        self.mge_file = f"{self.all_dir}/{self.prefix}/all_mge.tsv"
        self.all_motif_file = f"{self.work_dir}/all.motifs.csv"
        self.bin3c_cluster = f"{all_dir}/{prefix}/hic/bin3c_clust/clustering.mcl"
        self.contact_value_file = f"{all_dir}/{prefix}/hic/bin3c/contact_values.txt"
        self.spacer_linkage_file = f"{all_dir}/{prefix}/spacer/{prefix}_mge_spacer_hits.filter.tsv"
        self.genomad_plasmid = f"{self.all_dir}/{self.prefix}/Genomad/{self.prefix}.hifiasm.p_ctg.rename_summary/{self.prefix}.hifiasm.p_ctg.rename_plasmid_summary.tsv"
        self.genomad_virus = f"{self.all_dir}/{self.prefix}/Genomad/{self.prefix}.hifiasm.p_ctg.rename_summary/{self.prefix}.hifiasm.p_ctg.rename_virus_summary.tsv"
        
        self.mge_dict = None
        self.mge_genomad_dict = None
        self.depth_dict = None
        self.length_dict = None
        self.depth_cutoff = 10
        self.length_cutoff = 5000
        self.specificity_cutoff = 0.01
        self.final_score_cutoff = 0.5

    def get_unique_motifs(self, min_frac=0.4, min_sites = 100):
        if not os.path.exists(self.all_motif_file):
            print (f"[⚠️] Motif file not found: {self.all_motif_file}")
            return None, None
        df = pd.read_csv(self.all_motif_file)
        motif_num, unique_motifs, unique_motifs_identifier = get_unique_motifs(df, min_frac=min_frac, min_sites=min_sites)
        return motif_num, unique_motifs

    def simple_load_motifs(self):
        if not os.path.exists(self.all_motif_file):
            print (f"[⚠️] Motif file not found: {self.all_motif_file}")
            return None
        df = pd.read_csv(self.all_motif_file)
        return df

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
        linkage_info_list = []
        ## if self.host_sum_file not exists
        if not os.path.exists(self.host_sum_file):
            print(f"Host summary file {self.host_sum_file} not found.")
            return {}, {}, []
        line_num = sum(1 for line in open(self.host_sum_file) if line.strip())
        if line_num < 2:
            print(f"Host summary file {self.host_sum_file} is empty or has only header.")
            return {}, {}, []

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
            linkage_obj = Linkage_object()
            linkage_obj.load_from_row(row)
            linkage_info_list.append(linkage_obj)
            our_linkages[plasmid_name].append(bin_name)
            our_ctg_linkages[plasmid_name] = row['host']
            # our_linkages[row['plasmid']] = row['host']
        multiple_host_plasmid_num = 0
        for plasmid in our_linkages:
            if len(our_linkages[plasmid]) > 1:
                multiple_host_plasmid_num += 1
                # print (f"{plasmid} has multiple host: {our_linkages[plasmid]}")
        print (f"multiple host plasmid num: {multiple_host_plasmid_num} out of {len(our_linkages)}")
        return our_linkages, our_ctg_linkages, linkage_info_list

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

    def get_len_dict(self):
        with open(self.fai, "r") as f:
            for line in f:
                ctg, length, _, _, _ = line.strip().split("\t")
                length = int(length)
                self.length_dict[ctg] = length
        return self.length_dict

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
        self.length_dict = self.get_len_dict()
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

    def get_final_best_ctg(self, min_dp = 10):
        """
        Read all host contigs from a file.
        """
        best_ctgs = {}
        with open(self.all_host_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                ctg, ctg, domain = line.strip().split("\t")
                if self.depth_dict[ctg] < min_dp:
                    continue
                best_ctgs[ctg] = domain
        return best_ctgs
    
    def get_final_best_ctg2(self, min_dp = 10):
        """
        serve as high-quality contigs to show the landscape 
        of motif number distribution in envs
        """
        high_quality_ctgs = self.get_high_quality_ctgs()
        circular_ctgs = self.get_circular_ctgs()
        merged_ctgs = set(high_quality_ctgs) | set(circular_ctgs)
        final_selection = []
        genome_list = []
        taxa_dict = self.read_meta_gtdb()
        self.read_MGE()
        ## next examine if they have classified phylum
        for ctg in merged_ctgs:
            if ctg not in taxa_dict:
                continue
            if ctg in self.mge_dict:
                continue
            lineage = taxa_dict[ctg]
            phylum = classify_taxa(lineage, level="phylum")
            if phylum == "Unknown":
                continue
            if self.depth_dict[ctg] < min_dp:
                continue
            final_selection.append(ctg)
            genome = self.work_dir + "/contigs/" + ctg + ".fa"
            genome_list.append(genome)
        return genome_list, final_selection

    def get_high_dp_ctg_list(self, min_depth=10, min_len=100000):
        """
        Get high depth contig list.
        sever for drep, to explore motif variation 
        among strains
        """
        genome_list = []
        contig_list = []
        for contig in self.depth_dict:
            if self.depth_dict[contig] >= min_depth and self.length_dict[contig] >= min_len:
                genome = self.work_dir + "/contigs/" + contig + ".fa"
                genome_list.append(genome)
                contig_list.append(contig)
        return genome_list, contig_list

    def read_meta_gtdb(self):
        """
        Read the GTDB summary file and return a dictionary of contig to bin mapping.
        """
        isolation_taxa = {}
        if os.path.exists(self.gtdb):
            gtdb_df = pd.read_csv(self.gtdb, sep='\t')
            for index, row in gtdb_df.iterrows():
                anno = row['classification']
                isolation_taxa[row['user_genome']] = anno
        if os.path.exists(self.arc_gtdb):
            arc_gtdb_df = pd.read_csv(self.arc_gtdb, sep='\t')
            for index, row in arc_gtdb_df.iterrows():
                anno = row['classification']
                isolation_taxa[row['user_genome']] = anno
        return isolation_taxa

    def get_high_quality_ctgs(self):
        """
        Get the high quality contigs from a CheckM report.
        """
        high_quality_ctgs = []
        df = pd.read_csv(self.checkm, sep="\t", header=0)
        df = df[df['Completeness'] >= 50]
        df = df[df['Contamination'] <= 5]
        for ctg in df['Name']:
            high_quality_ctgs.append(ctg)

        print (f"Total {len(high_quality_ctgs)} high quality contigs found.")
        return high_quality_ctgs

    def get_circular_ctgs(self):
        """
        Get the best contig based on length from a fasta file.
        """
        circular_ctgs = []
        with open(self.fai, "r") as f:
            for line in f:
                ctg, length, _, _, _ = line.strip().split("\t")
                if ctg[-1] == "C":
                    circular_ctgs.append(ctg)
        return circular_ctgs

    def read_genomad(self, genomad_plasmid, mge_type="plasmid"):
        print (f"Reading {genomad_plasmid}...")
        genomad_dict = {}
        genomad = pd.read_csv(genomad_plasmid, sep = "\t")
        for i, row in genomad.iterrows():
            
            if re.search('\|provirus', row['seq_name']):
                continue
            if row['seq_name'] == 'seq_name':
                continue
            mge = mge_obj(
                name=row['seq_name'],
                type=mge_type,
                length=row['length'],
                method="genomad"
            )
            genomad_dict[row['seq_name']] = mge
        return genomad_dict

    def collect_all_mges(self):
        if os.path.exists(self.genomad_plasmid):
            genomad_plasmid = self.read_genomad(self.genomad_plasmid, mge_type="plasmid")
        else:
            print(f"Warning: {self.genomad_plasmid} does not exist. Skipping Genomad plasmid classification.")
            genomad_plasmid = {}
        if os.path.exists(self.genomad_virus):
            genomad_virus = self.read_genomad(self.genomad_virus, mge_type="virus")
        else:
            print(f"Warning: {self.genomad_virus} does not exist. Skipping Genomad virus classification.")
            genomad_virus = {}
        self.mge_genomad_dict = {**genomad_plasmid, **genomad_virus}
        return self.mge_genomad_dict


class Linkage_object(object):


    def __init__(self):
        self.mge = None
        self.host = None
        self.specificity = None
        self.final_score = None
        self.pvalue = None
        self.MGE_gc = None
        self.host_gc = None
        self.cos_sim = None
        self.MGE_cov = None
        self.host_cov = None
        self.mge_len = None

    def load_from_row(self, row):
        self.mge = row['MGE']
        self.host = row['host']
        self.specificity = row['specificity']
        self.final_score = row['final_score']
        self.pvalue = row['pvalue']
        self.MGE_gc = row['MGE_gc']
        self.host_gc = row['host_gc']
        self.cos_sim = row['cos_sim']
        self.MGE_cov = row['MGE_cov']
        self.host_cov = row['host_cov']
        self.mge_len = row['MGE_len']

class Isolation_sample(My_sample):

    def __init__(self, prefix, all_dir):
        super().__init__(prefix, all_dir)
        self.phylum = None
        self.species = None
        self.lineages = None
        self.average_dp = None
        self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation4"
        self.depth_file = os.path.join(self.work_dir, "mean_depth.csv")
        self.host_sum_file = os.path.join(self.work_dir, "host_summary.csv")
        self.orphan_file = os.path.join(self.work_dir, "regulatory_motif_enrichment.csv")
        self.motif_freq_file = os.path.join(self.work_dir, "motif_length_stats.csv")
        self.profile = os.path.join(self.work_dir, "motif_profile.csv")
        self.gtdb = os.path.join(self.work_dir, "../GTDB_2/gtdbtk.bac120.summary.tsv")
        self.all_gtdb_ark = os.path.join(self.work_dir, "../GTDB/gtdbtk.ar53.summary.tsv")
        self.all_gtdb = os.path.join(self.work_dir, "../GTDB/gtdbtk.bac120.summary.tsv")
        self.all_motif_file = f"{self.work_dir}/all.motifs.csv"
        self.checkm = os.path.join(self.work_dir, "../checkM2/quality_report.tsv")
        self.isolation_RM_file = f"{self.work_dir}/RM_systems/all_ctgs_RM.rm.genes.tsv"
        self.isolation_RM_dict = None

    def load_isolation_RM(self):
        self.isolation_RM_dict = defaultdict(list)
        if not os.path.exists(self.isolation_RM_file):
            print (f"[⚠️] RM file not found: {self.isolation_RM_file}")
            return self.isolation_RM_dict
        f = open(self.isolation_RM_file, "r")
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            if fields[3] == "RE" or fields[3] == "SP":
                continue
            system_name = fields[0]
            if not re.search("RM Operon #",system_name) and not re.search("Singleton #", system_name):
                continue
            gene_name = fields[1]
            self.isolation_RM_dict[system_name].append(gene_name)
        f.close()
        return self.isolation_RM_dict

    def get_phylum(self):
        isolation_taxa = self.read_isolation_gtdb()
        self.lineage = isolation_taxa[list(isolation_taxa.keys())[0]]
        self.phylum = self.lineage.split(";")[1][3:] if ";" in self.lineage else "Unclassified"
        self.species = self.lineage.split(";")[-1][3:] if ";" in self.lineage else "Unclassified"
    
    def read_isolation_gtdb(self):
        """
        Read the GTDB summary file and return a dictionary of contig to bin mapping.
        """
            
        isolation_taxa = {}
        if os.path.exists(self.all_gtdb):
            gtdb_df = pd.read_csv(self.all_gtdb, sep='\t')
            for index, row in gtdb_df.iterrows():
                anno = row['classification']
                sra_id = row['user_genome']
                isolation_taxa[sra_id] = anno

        if os.path.exists(self.all_gtdb_ark):
            arc_gtdb_df = pd.read_csv(self.all_gtdb_ark, sep='\t')
            for index, row in arc_gtdb_df.iterrows():
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

    def get_average_depth(self):
        ## depth is weighted by length
        total_depth = 0
        total_length = 0
        df = pd.read_csv(self.depth_file)

        for _, row in df.iterrows():
            contig = row['contig']
            depth = float(row['depth'])

            length = int(row['length'])
            total_depth += depth * length
            total_length += length
        if total_length == 0:
            return 0
        self.average_dp = total_depth / total_length
        return self.average_dp

    def search_megaP(self, length_cutoff=500000):
        potential_megaP = []
        for mge in self.mge_dict:
            mge_length = self.length_dict.get(mge, 0)
            if mge_length >= length_cutoff:
                potential_megaP.append((mge, mge_length, self.mge_dict[mge], self.lineage))
        return potential_megaP

    def get_iso_good_ctgs(self, min_depth=10, min_len=500000):
        """
        Get high depth contig list.
        sever for drep, to explore motif variation 
        among strains
        """
        genome_list = []
        contig_list = []
        ## sort self.length_dict by length descending
        sorted_contigs = sorted(self.length_dict.items(), key=lambda x: x[1], reverse=True)
        ## select one longest contig for each isolation sample
        for contig, length in sorted_contigs:
            if self.depth_dict.get(contig, 0) >= min_depth and length >= min_len\
                and contig not in self.mge_dict:
                genome = self.work_dir + "/contigs/" + contig + ".fa"
                genome_list.append(genome)
                contig_list.append(contig)
            else:
                print ("Skipping contig:", contig, "Depth:", self.depth_dict.get(contig, 0), "Length:", length)
            break
        return genome_list, contig_list

class My_contig(My_sample):

    def __init__(self, prefix, all_dir, contig, data_type="meta"):
        super().__init__(prefix, all_dir, data_type)
        self.contig = contig
        if data_type == "meta":
            self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation4"
        else:
            self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation4"
        self.ctg_ref = f"{self.work_dir}/contigs/{contig}.fa"
        self.ctg_ref_fai = f"{self.work_dir}/contigs/{contig}.fa.fai"
        self.gff = f"{self.work_dir}/gffs/{contig}.gff"
        self.reprocess_gff = f"{self.work_dir}/gffs/{contig}.reprocess.gff"
        self.ipd_ratio_file = f"{self.work_dir}/ipd_ratio/{contig}.ipd3.csv"
        self.motif_file = f"{self.work_dir}/motifs/{contig}.motifs.csv"
        self.RM_file = f"{self.work_dir}/RM_systems/{contig}.RM.csv"
        self.ctg_profile = f"{self.work_dir}/profiles/{contig}.motifs.profile.csv"
        self.ctg_len = None
        self.modified_ratio = None
        self.modified_motif_ratio = None
        self.motif_ratio = None
        self.modified_num = None
        self.modified_motif_num = None
        self.RM_dict = None
        self.modified_loci = None
        self.REF = None

    # def sample_name_corr(self):
    #     soil_dict = {
    #         "":"soil_1"
    #     }

    def read_motif(self, min_frac=0.3, min_sites=30):
        ## check if file exists
        if not os.path.exists(self.motif_file):
            return None
        motif_df = pd.read_csv(self.motif_file)
        motif_df = motif_df[(motif_df['fraction'] >= min_frac) & (motif_df['nDetected'] >= min_sites)]
        return motif_df

    def get_ctg_len(self):
        if self.ctg_len is not None:
            return self.ctg_len
        if not os.path.exists(self.ctg_ref_fai):
            print (f"[⚠️] FAI file not found: {self.ctg_ref_fai}")
            return 0
        with open(self.ctg_ref_fai, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[0] == self.contig:
                    self.ctg_len = int(fields[1])
                    return self.ctg_len
        return 0

    def get_mod_ratio(self, score_cutoff = 30):
        self.get_ctg_len()
        self.modified_num = 0
        self.modified_motif_num = 0
        if os.path.exists(self.reprocess_gff):
            for line in open(self.reprocess_gff, "r"):
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if int(fields[5]) < score_cutoff:
                    continue
                self.modified_num += 1
                if re.search(";motif=", fields[8]):
                    self.modified_motif_num += 1
        self.modified_ratio = self.modified_num / self.ctg_len
        self.modified_motif_ratio = self.modified_motif_num/ self.ctg_len
        self.motif_ratio =  self.modified_motif_num / self.modified_num if self.modified_num > 0 else 0
        return self.modified_num, self.modified_motif_num, self.modified_ratio, self.modified_motif_ratio, self.motif_ratio

    def get_borg_mod_ratio(self, max_len=float('inf'), score_cutoff = 30):
        # self.get_ctg_len()
        print (self.reprocess_gff)
        self.ctg_len = max_len
        self.modified_num = 0
        self.modified_motif_num = 0
        modified_sites = set()
        if os.path.exists(self.reprocess_gff):
            for line in open(self.reprocess_gff, "r"):
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if int(fields[3]) >= max_len:
                    break
                if int(fields[5]) < score_cutoff:
                    continue
                self.modified_num += 1
                if re.search(";motif=", fields[8]):
                    self.modified_motif_num += 1
                    modified_sites.add(fields[3])
        self.modified_ratio = self.modified_num / self.ctg_len
        self.modified_motif_ratio = self.modified_motif_num/ self.ctg_len
        self.motif_ratio =  self.modified_motif_num / self.modified_num if self.modified_num > 0 else 0
        return modified_sites

    def get_modified_loci(self, score_cutoff = 30):
        ## read the gff file
        f = open(self.reprocess_gff, "r")
        self.modified_loci = {}
        for line in f:
            if line[0] == "#":
                continue
            line = line.strip().split("\t")

            ref = line[0]
            pos = int(line[3]) 
            strand = line[6]
            score = int(line[5])
            if score <= score_cutoff:
                continue
            self.modified_loci[ref + ":" + str(pos) + strand] = score
        print ("no. of modified loci", len(self.modified_loci))
        return self.modified_loci

    def read_ref(self):
        self.REF = {}
        for record in SeqIO.parse(self.ctg_ref, "fasta"):
        #     seq_dict[record.id] = record.seq
        # return seq_dict
            self.REF[record.id] = record.seq
            # return str(record.seq), record.id
        return self.REF

    def count_mod_frac_in_motif(self, motif_new, exact_pos, start, end, score_cutoff=30):
        motif_len = len(motif_new)
        rev_exact_pos = motif_len - exact_pos + 1
        motif_sites = {}
        motif_loci_num = 0
        motif_modify_num = 0
        for_loci_num = 0
        rev_loci_num = 0
        for_modified_num = 0
        rev_modified_num = 0

        # motif_ipd_ratio = []
        record_modified_sites = {}

        for r, contig in self.REF.items():
            for site in nt_search(str(contig)[start:end], motif_new)[1:]:
                tag = r + ":" + str(site+exact_pos) + "+"
                record_modified_sites[tag] = motif_new

                motif_loci_num += 1
                for_loci_num += 1
                if tag in self.modified_loci:
                    motif_modify_num += 1
                    for_modified_num += 1


            for site in nt_search(str(contig)[start:end], Seq(motif_new).reverse_complement())[
                1:
            ]:
                tag = r + ":" + str(site+rev_exact_pos) + "-"
                record_modified_sites[tag] = motif_new
                # for i in range(site-100, site+100):
                #     record_modified_sites[r + ":" + str(i) + "+"] = motif_new
                motif_loci_num += 1
                rev_loci_num += 1
                if tag in self.modified_loci:
                    motif_modify_num += 1
                    rev_modified_num += 1

        if motif_loci_num == 0:
            ratio = 0
        else:
            ratio = motif_modify_num/motif_loci_num
        if for_loci_num == 0:
            for_ratio = 0
        else:
            for_ratio = for_modified_num/for_loci_num
        if rev_loci_num == 0:
            rev_ratio = 0
        else:
            rev_ratio = rev_modified_num/rev_loci_num


        return motif_loci_num, motif_modify_num, ratio

    def count_mod_dist(self, motif_new, exact_pos, start, end):

        motif_len = len(motif_new)
        rev_exact_pos = motif_len - exact_pos + 1

        for_loci_num = 0
        rev_loci_num = 0

        for r, contig in self.REF.items():
            for site in nt_search(str(contig)[start:end], motif_new)[1:]:
                tag = r + ":" + str(site+exact_pos) + "+"
                for_loci_num += 1

            for site in nt_search(str(contig)[start:end], Seq(motif_new).reverse_complement())[
                1:
            ]:
                tag = r + ":" + str(site+rev_exact_pos) + "-"
                rev_loci_num += 1
        total_loci = for_loci_num + rev_loci_num
        total_ratio = total_loci / (end - start)
        for_ratio = for_loci_num / (end - start)
        rev_ratio = rev_loci_num / (end - start)
        print (f"Motif: {motif_new}, Total loci: {total_loci} ({for_loci_num} forward, \
               {rev_loci_num} reverse) in range {start}-{end}\
                , Ratios: {total_ratio:.4f} ({for_ratio:.4f} forward, {rev_ratio:.4f} reverse)")


    def load_RM(self):
        self.RM_dict = defaultdict(list)
        if not os.path.exists(self.RM_file):
            print (f"[⚠️] RM file not found: {self.RM_file}")
            return self.RM_dict
        f = open(self.RM_file, "r")
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            if fields[3] == "RE" or fields[3] == "SP":
                continue
            system_name = fields[0]
            if not re.search("RM Operon #",system_name) and not re.search("Singleton #", system_name):
                continue
            gene_name = fields[1]
            self.RM_dict[system_name].append(gene_name)
        f.close()
        return self.RM_dict

    def read_ctg_profile(self):
        if not os.path.exists(self.ctg_profile):
            print (f"[⚠️] Contig profile file not found: {self.ctg_profile}")
            return None
        df = pd.read_csv(self.ctg_profile)
        motif_loci_num_dict = {}
        for index, row in df.iterrows():
            motif_loci_num_dict[row['motifString']] = row['motif_loci_num']
        return motif_loci_num_dict

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
        # contig_list = df_pivot.index.tolist()
        contig_list = self.members
        profile_matrix = df_pivot.values
        if len(df_pivot) < 2:
            print("Warning: Insufficient data for pairwise comparison.")
            return [], 0
        similarity_matrix = cosine_similarity(profile_matrix)
        jaccard_matrix = np.zeros((len(contig_list), len(contig_list)))
        ## also count Jaccard similarity, binary presence/absence with cutoff 0.3
        binary_df = df_pivot.applymap(lambda x: 1 if x >= bin_freq else 0)
        binary_matrix = binary_df.values

        n = len(contig_list)
        data = []
        for i in range(n):
            for j in range(i+1, n):
                jaccard = jaccard_score(binary_matrix[i], binary_matrix[j])
                jaccard_matrix[i][j] = jaccard
                jaccard_matrix[j][i] = jaccard
                data.append([similarity_matrix[i][j], jaccard])
        return data, self.uniq_clade_num(binary_matrix), jaccard_matrix

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

    def pairwise_edit_distance(self, seq_dir, edit_dir):
        dnadiff_list = []
        dnadiff_mat = np.zeros((len(self.members), len(self.members)))
        for i in range(len(self.members)):
            for j in range(i+1, len(self.members)):
                # compute edit distance between self.members[i] and self.members[j]
                
                genome_1_fasta = os.path.join(seq_dir, self.cluster, f"{self.members[i]}.fa")
                genome_2_fasta = os.path.join(seq_dir, self.cluster, f"{self.members[j]}.fa")
                # print (f"Edit distance between {self.members[i]} and {self.members[j]}: ")
                prefix = os.path.join(edit_dir, f"{self.members[i]}_{self.members[j]}")
                cmd = f"dnadiff -p {prefix} {genome_1_fasta} {genome_2_fasta}"
                
                dnadiff_report = f"{prefix}.report"
                if not os.path.exists(dnadiff_report):
                    # os.system(cmd)
                    dnadiff_list.append(cmd)
                    print ("## dnadiff report not found, adding command to list:", cmd)
            # else:
                total_snps, total_indels, edit_distance = read_dnadiff_report(dnadiff_report)
                dnadiff_mat[i][j] = edit_distance
                dnadiff_mat[j][i] = edit_distance
                # print (f"Total SNPs: {total_snps}, Total Indels: {total_indels}, Edit Distance: {edit_distance}")
                # ## remove tmp files except the report
                # for file in os.listdir(edit_dir):
                #     if file.startswith(f"{self.members[i]}_{self.members[j]}") and not file.endswith(".report"):
                #         os.remove(os.path.join(edit_dir, file))
        return dnadiff_list, dnadiff_mat

def read_dnadiff_report(dnadiff_report):
    ## read total TotalSNPs and TotalIndels from dnadiff report
    total_snps = 0
    total_indels = 0
    with open(dnadiff_report, 'r') as f:
        for line in f:
            if line.startswith("TotalSNPs"):
                total_snps = int(line.strip().split()[1])
            elif line.startswith("TotalIndels"):
                total_indels = int(line.strip().split()[1])
    edit_distance = total_snps + total_indels
    return total_snps, total_indels, edit_distance

class My_gene(object):

    def __init__(self, gff, genome, genome_file):
        self.gff = gff
        self.genome = genome
        self.genome_file = genome_file
        self.Gene_regulatory_regions = None
        self.cds_regions = None
        self.modified_regions = None
        self.genome_length = None

    def collect_regulation_region(self):
        import pandas as pd
        
        # Read GFF file efficiently using pandas
        try:
            # Define column names for GFF format
            gff_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
            
            # Read the file, manually filtering out comment lines (# and ##)
            valid_lines = []
            with open(self.gff, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        break
                    if line and not line.startswith('#'):
                        valid_lines.append(line.split('\t'))
            
            if not valid_lines:
                self.Gene_regulatory_regions = []
                self.cds_regions = []
                return self.Gene_regulatory_regions
            
            # Create DataFrame from filtered lines
            df = pd.DataFrame(valid_lines, columns=gff_columns)
            df['start'] = pd.to_numeric(df['start'], errors='coerce')
            df['end'] = pd.to_numeric(df['end'], errors='coerce')
            
            # Remove rows where start/end conversion failed
            df = df.dropna(subset=['start', 'end'])
            df['start'] = df['start'].astype(int)
            df['end'] = df['end'].astype(int)
            
            # Filter for the specific genome and CDS features
            cds_df = df[(df['seqname'] == self.genome) & (df['feature'] == 'CDS')].copy()
            
            if cds_df.empty:
                self.Gene_regulatory_regions = []
                self.cds_regions = []
                return self.Gene_regulatory_regions
            
            # Store CDS regions as list of tuples for compatibility
            self.cds_regions = list(zip(cds_df['start'], cds_df['end'], cds_df['strand']))
            
            # Calculate regulatory regions vectorized
            pos_strand = cds_df['strand'] == '+'
            neg_strand = cds_df['strand'] == '-'
            
            # For positive strand: 100bp upstream to 50bp downstream of start
            pos_reg_start = np.maximum(1, cds_df.loc[pos_strand, 'start'] - 100)
            pos_reg_end = cds_df.loc[pos_strand, 'start'] + 50
            
            # For negative strand: 50bp upstream to 100bp downstream of end
            neg_reg_start = np.maximum(1, cds_df.loc[neg_strand, 'end'] - 50)
            neg_reg_end = cds_df.loc[neg_strand, 'end'] + 100
            
            # Combine regulatory regions
            reg_regions = []
            if not pos_reg_start.empty:
                reg_regions.extend(list(zip(pos_reg_start, pos_reg_end, ['+'] * len(pos_reg_start))))
            if not neg_reg_start.empty:
                reg_regions.extend(list(zip(neg_reg_start, neg_reg_end, ['-'] * len(neg_reg_start))))
            
            self.Gene_regulatory_regions = reg_regions
            
        except Exception as e:
            print(f"Error reading GFF file {self.gff}: {e}")
            # Fallback to original method
            self.Gene_regulatory_regions = []
            self.cds_regions = []
            for line in open(self.gff, "r"):
                if line.startswith("#"):
                    continue
                if line.startswith(">"):
                    break
                line = line.strip().split("\t")
                if len(line) < 9:
                    continue
                if line[0] != self.genome:
                    continue
                if line[2] == "CDS":
                    start = int(line[3])
                    end = int(line[4])
                    strand = line[6]
                    self.cds_regions.append((start, end, strand))
                    if strand == "+":
                        reg_start = max(1, start - 100)
                        reg_end = start + 50
                        self.Gene_regulatory_regions.append((reg_start, reg_end, "+"))
                    elif strand == "-":
                        reg_start = max(1, end - 50)
                        reg_end = end + 100
                        self.Gene_regulatory_regions.append((reg_start, reg_end, "-"))
        
        return self.Gene_regulatory_regions

    def read_modified_gff(self, modified_gff, min_score=30):
        import pandas as pd
        
        try:
            # Define column names for GFF format
            gff_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
            
            # Read the file, manually filtering out comment lines (# and ##)
            valid_lines = []
            with open(modified_gff, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        valid_lines.append(line.split('\t'))
            
            if not valid_lines:
                return []
            
            # Create DataFrame from filtered lines
            df = pd.DataFrame(valid_lines, columns=gff_columns)
            df['start'] = pd.to_numeric(df['start'], errors='coerce')
            df['score'] = pd.to_numeric(df['score'], errors='coerce')
            
            # Remove rows where start/score conversion failed
            df = df.dropna(subset=['start', 'score'])
            df['start'] = df['start'].astype(int)
            df['score'] = df['score'].astype(float)
            
            # Filter in one operation
            filtered_df = df[
                (df['seqname'] == self.genome) & 
                (df['feature'] == 'modified_base') & 
                (df['score'] >= min_score)
            ]
            
            if filtered_df.empty:
                return []
            
            # Extract positions and strands as list of lists for compatibility
            self.modified_regions = filtered_df[['start', 'strand']].values.tolist()
            
        except Exception as e:
            print(f"Error reading modified GFF file {modified_gff}: {e}")
            # Fallback to original method
            modified_regions = []
            for line in open(modified_gff, "r"):
                if line.startswith("#"):
                    continue
                line = line.strip().split("\t")
                if len(line) < 9:
                    continue
                if line[0] != self.genome:
                    continue
                if int(line[5]) < min_score:
                    continue
                if line[2] == "modified_base":
                    pos = int(line[3])
                    strand = line[6]
                    self.modified_regions.append([pos, strand])

        return self.modified_regions

    def get_genome_len(self):
        ## count from fai file
        fai = open(self.genome_file + ".fai", "r")
        self.genome_length = 0
        for line in fai:
            line = line.strip().split("\t")
            self.genome_length += int(line[1])
        fai.close()
        return self.genome_length


    def intersect(self):

        ## usse dict to record all these info, no need for positions
        region_info = {"regulatory": {"count": 0},
                       "cds": {"count": 0},
                       "non_coding": {"count": 0}}
        genome_length = self.get_genome_len()
        for mod in self.modified_regions:
            pos = mod[0]
            in_regulatory = False
            in_cds = False
            for reg in self.Gene_regulatory_regions:
                if pos >= reg[0] and pos <= reg[1]:
                    region_info["regulatory"]["count"] += 1
                    in_regulatory = True
                    break

            for cds in self.cds_regions:
                if pos >= cds[0] and pos <= cds[1]:
                    region_info["cds"]["count"] += 1
                    in_cds = True
                    break
            if not in_cds:
                region_info["non_coding"]["count"] += 1
        # calculate frequencies
        non_coding_length, total_cds_length, total_regulatory_length = self.get_non_coding_length()
        region_info["regulatory"]["frequency"] = region_info["regulatory"]["count"] / total_regulatory_length if total_regulatory_length > 0 else 0
        region_info["cds"]["frequency"] = region_info["cds"]["count"] / total_cds_length if total_cds_length > 0 else 0
        region_info["non_coding"]["frequency"] = region_info["non_coding"]["count"] / non_coding_length if non_coding_length > 0 else 0
        ## transfer region_info to df, also add contig name, genome length
        region_info["genome"] = self.genome
        region_info["genome_length"] = genome_length
        mod_count = len(self.modified_regions)
        general_freq = mod_count / genome_length if genome_length > 0 else 0
        ## region_info to df
        data = [self.genome, region_info["genome_length"], region_info["regulatory"]["count"], 
                total_regulatory_length, region_info["regulatory"]["frequency"],
                region_info["cds"]["count"], total_cds_length, region_info["cds"]["frequency"],
                region_info["non_coding"]["count"], non_coding_length, region_info["non_coding"]["frequency"],
                mod_count, general_freq
                ]

        df = pd.DataFrame([data], columns=['genome', 'genome_length', 'regulatory_count', 'regulatory_length', 'regulatory_frequency',
                                           'cds_count', 'cds_length', 'cds_frequency',
                                           'non_coding_count', 'non_coding_length', 'non_coding_frequency',
                                           'modified_count', 'general_frequency'])
        return df

    def get_non_coding_length(self):
        ## as cds region might overlap with each other, first merge cds regions
        merged_cds = []
        sorted_cds = sorted(self.cds_regions, key=lambda x: x[0])
        for cds in sorted_cds:
            if not merged_cds:
                merged_cds.append(cds)
            else:
                last_cds = merged_cds[-1]
                if cds[0] <= last_cds[1]:
                    # overlap
                    merged_cds[-1] = (last_cds[0], max(last_cds[1], cds[1]), last_cds[2])
                else:
                    merged_cds.append(cds)
        total_cds_length = 0
        for cds in merged_cds:
            total_cds_length += (cds[1] - cds[0] + 1)
        non_coding_length = self.genome_length - total_cds_length
        ## also merge regulatory regions to get total regulatory length
        merged_regulatory = []
        sorted_regulatory = sorted(self.Gene_regulatory_regions, key=lambda x: x[0])
        for reg in sorted_regulatory:
            if not merged_regulatory:
                merged_regulatory.append(reg)
            else:
                last_reg = merged_regulatory[-1]
                if reg[0] <= last_reg[1]:
                    # overlap
                    merged_regulatory[-1] = (last_reg[0], max(last_reg[1], reg[1]), last_reg[2])
                else:
                    merged_regulatory.append(reg)
        total_regulatory_length = 0
        for reg in merged_regulatory:
            total_regulatory_length += (reg[1] - reg[0] + 1)
        return non_coding_length, total_cds_length, total_regulatory_length

    def intersect_gene(self):
        genome_length = self.get_genome_len()
        
        # Initialize counters for four regions
        region_counts = {
            "upstream_150": 0,      # CDS start - 150bp to CDS start
            "start_150": 0,         # CDS start to CDS start + 150bp
            "end_150": 0,           # CDS end - 150bp to CDS end
            "downstream_150": 0     # CDS end to CDS end + 150bp
        }
        
        # Assume all regions have the same length (150bp)
        # Total region lengths for frequency calculation
        num_cds = len(self.cds_regions)
        total_lengths = {
            "upstream_150": 150 * num_cds,
            "start_150": 150 * num_cds,
            "end_150": 150 * num_cds,
            "downstream_150": 150 * num_cds
        }
        
        # Count modifications in each region
        for mod in self.modified_regions:
            pos = mod[0]
            mod_strand = mod[1]
            
            for cds in self.cds_regions:
                start, end, strand = cds[0], cds[1], cds[2]
                
                # Only count modifications on the same strand as the gene
                if mod_strand != strand:
                    continue
                
                if strand == '+':
                    # For + strand genes (5' to 3' left to right)
                    # Upstream 150bp: [start - 150, start)
                    upstream_start = max(1, start - 150)
                    if upstream_start <= pos < start:
                        region_counts["upstream_150"] += 1
                        break
                    
                    # Start region: [start, start + 150]
                    start_region_end = min(end, start + 150)
                    if start <= pos <= start_region_end:
                        region_counts["start_150"] += 1
                        break
                    
                    # End region: [end - 150, end]
                    end_region_start = max(start, end - 150)
                    if end_region_start <= pos <= end:
                        region_counts["end_150"] += 1
                        break
                    
                    # Downstream 150bp: (end, end + 150]
                    downstream_end = end + 150
                    if end < pos <= downstream_end:
                        region_counts["downstream_150"] += 1
                        break
                
                else:  # strand == '-'
                    # For - strand genes (5' to 3' right to left)
                    # Upstream 150bp: (end, end + 150] (upstream of transcription start)
                    upstream_end = end + 150
                    if end < pos <= upstream_end:
                        region_counts["upstream_150"] += 1
                        break
                    
                    # Start region: [end - 150, end] (transcription start region)
                    start_region_start = max(start, end - 150)
                    if start_region_start <= pos <= end:
                        region_counts["start_150"] += 1
                        break
                    
                    # End region: [start, start + 150] (transcription end region)
                    end_region_end = min(end, start + 150)
                    if start <= pos <= end_region_end:
                        region_counts["end_150"] += 1
                        break
                    
                    # Downstream 150bp: [start - 150, start) (downstream of transcription end)
                    downstream_start = max(1, start - 150)
                    if downstream_start <= pos < start:
                        region_counts["downstream_150"] += 1
                        break
        
        # Calculate frequencies
        region_frequencies = {}
        for region in region_counts:
            if total_lengths[region] > 0:
                region_frequencies[region] = region_counts[region] / total_lengths[region]
            else:
                region_frequencies[region] = 0
        
        # Create DataFrame with results - only genome, genome_length, and frequencies
        data = [
            self.genome,
            genome_length,
            region_frequencies["upstream_150"],
            region_frequencies["start_150"],
            region_frequencies["end_150"],
            region_frequencies["downstream_150"]
        ]
        
        df = pd.DataFrame([data], columns=[
            'genome', 'genome_length',
            'upstream_150_frequency',
            'start_150_frequency',
            'end_150_frequency',
            'downstream_150_frequency'
        ])
        
        return df

if __name__ == "__main__":
    pass

    gff = "/home/shuaiw/borg/paper/gene_anno/meta/soil_1_129_C/prokka/soil_1_129_C.gff"
    modified_gff = "/home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/gffs/soil_1_129_C.reprocess.gff"
    genome_file = "/home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/contigs/soil_1_129_C.fa"
    count_dir = "/home/shuaiw/borg/paper/gene_anno/meta/soil_1_129_C/count/"
    count_file = os.path.join(count_dir, "soil_1_129_C_region_count.csv")   
    ## mkdir count_dir if not exists
    if not os.path.exists(count_dir):
        os.makedirs(count_dir)
    my_gene = My_gene(gff, "soil_1_129_C", genome_file)
    my_gene.collect_regulation_region()
    my_gene.read_modified_gff(modified_gff)
    region_info = my_gene.intersect()
    print(region_info)
    region_info.to_csv(count_file, index=False)