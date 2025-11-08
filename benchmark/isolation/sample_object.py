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


class My_sample(object):
    def __init__(self, prefix, all_dir):
        self.prefix = prefix
        self.all_dir = all_dir

        self.work_dir = f"{self.all_dir}/{self.prefix}/{self.prefix}_methylation2"
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

    def get_unique_motifs(self):
        df = pd.read_csv(self.all_motif_file)
        return get_unique_motifs(df)

    def read_MGE(self):
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
        line_num = sum(1 for line in open(self.host_sum_file) if line.strip())
        if line_num > 1:
            host_sum = pd.read_csv(self.host_sum_file)
            all_link_df = pd.concat([all_link_df, host_sum], axis=0)

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
            required_cols = ['contig', 'depth', 'length']
            if not all(col in df.columns for col in required_cols):
                print(f"[⚠️] Depth file missing required columns. Expected: {required_cols}, Found: {list(df.columns)}")
                return self.depth_dict, self.length_dict
            
            # Create dictionaries
            for _, row in df.iterrows():
                contig = row['contig']
                depth = float(row['depth'])
                length = int(row['length'])

                self.depth_dict[contig] = depth
                self.length_dict[contig] = length

            # print(f"[✔] Successfully read {len(self.depth_dict)} contigs from depth file")
            # print(f"    Depth range: {min(self.depth_dict.values()):.2f} - {max(self.depth_dict.values()):.2f}")
            # print(f"    Length range: {min(self.length_dict.values())} - {max(self.length_dict.values())}")
            
        except Exception as e:
            print(f"[⚠️] Error reading depth file {self.depth_file}: {str(e)}")
            return {}, {}
        
        return self.depth_dict, self.length_dict



class Isolation_sample(My_sample):

    def __init__(self, prefix, all_dir):
        super().__init__(prefix, all_dir)
        self.phylum = None
        self.lineages = None
        
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
            