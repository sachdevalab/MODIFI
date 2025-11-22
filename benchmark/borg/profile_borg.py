import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
from sample_object import My_sample, Isolation_sample
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


if __name__ == "__main__":
    borg_file = "all_borg_contigs_summary.tsv"
    all_dir = "/home/shuaiw/borg/paper/run2/"
    
    # Load BORG data
    borg_data = My_Borg(borg_file)
    high_dp_borgs = borg_data.get_high_depth_borgs(min_depth=5.0)
    
    members = []
    for i, entry in enumerate(high_dp_borgs):
        print(f"{i+1}. {entry}")
        members.append(entry.seq_name)
    print (members)
    # members = ["soil_1_1336_L", "soil_s4_1_109_C"]

    seq_dir = "/home/shuaiw/borg/paper/borg_data/profile/"
    cluster = "profile"
    plot_name = os.path.join(seq_dir, f"borg_motif_profile.pdf")
    cluster_species = "borg & hosts"
    cluster_obj = given_species_drep(all_dir, members, seq_dir, cluster,
                                    seq_dir, seq_dir, min_frac=0.3, 
                                    min_sites=10, score_cutoff = 30)
    cluster_obj.plot_profile(cluster, plot_name, cluster_species)