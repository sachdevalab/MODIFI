"""
data object for the sample
"""

import os
import sys
import pandas as pd
from Bio.Seq import Seq

def get_unique_motifs(df_motif):
    df_motif = df_motif[(df_motif['fraction'] >= 0.4) & (df_motif['nDetected'] >= 100)]
    ## rm redundant motifs which are reverse complement 
    unique_motifs = []
    for index, row in df_motif.iterrows():
        if row['motifString'] not in unique_motifs and  str(Seq(row['motifString']).reverse_complement()) not in unique_motifs:
            unique_motifs.append(row['motifString'])
    return len(unique_motifs), unique_motifs

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

def get_best_ctg(depth_file, fai, min_len = 1000000):
    """
    Get the best contig based on length from a fasta file.
    """
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    length_dict = {}
    for index, row in depth_df.iterrows():

        if row['depth'] >= 10:
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
        
        self.mge_dict = None
        self.depth_dict = None
        self.length_dict = None
        self.depth_cutoff = 10
        self.length_cutoff = 5000

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
            if row["specificity"] >= 0.01:
                continue
            if row['final_score'] <= 0.5:
                continue
            linkage_num += 1
        return linkage_num

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

class Isolation_sample(My_sample):

    def __init__(self, prefix, all_dir):
        super().__init__(prefix, all_dir)
        self.phylum = None
        self.lineages = None
        self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation2"
        
    def get_phylum(self):
        isolation_taxa = self.read_isolation_gtdb()
        self.lineage = isolation_taxa[list(isolation_taxa.keys())[0]]
        self.phylum = self.lineage.split(";")[1][3:] if ";" in self.lineage else "Unclassified"

    def read_isolation_gtdb(self):
        """
        Read the GTDB summary file and return a dictionary of contig to bin mapping.
        """
        gtdb_df = pd.read_csv(self.gtdb, sep='\t')
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