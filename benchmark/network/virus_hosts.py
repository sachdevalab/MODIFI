import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import plotly.graph_objects as go
from collections import defaultdict    
import re
    

def read_gtdb(gatk):
    """
    Read the GTDB summary file and return a dictionary of contig to bin mapping.
    """
    gtdb_df = pd.read_csv(gatk, sep='\t')
    bin2anno_dict = {}
    for index, row in gtdb_df.iterrows():
        if row['user_genome'] != 'user_genome':
            anno = row['classification']
            # print (f"anno: {anno}")
            # if re.search('Unclassified', anno):
            #     phylum = "Unclassified"
            # else:
            #     phylum = anno.split(';')[1].strip()
            bin2anno_dict[row['user_genome']] = anno
    return bin2anno_dict

def read_summary(host_sum_file, bin2anno_dict, MGE_type_dict, cutoff = 0.6):
    host_sum = pd.read_csv(host_sum_file)
    mge_hosts_dict = defaultdict(set)
    mge_contigs_dict = defaultdict(set)
    ## get node file
    for index, row in host_sum.iterrows():
        if row["pvalue"] >= 0.05:
            continue
        if row['final_score'] <= cutoff:
            continue
        mge = row['MGE'] 
        if MGE_type_dict[mge] != "virus":
            continue
        each_mge_file = f"{all_dir}/{prefix}/{prefix}_methylation3/hosts/{mge}.host_prediction.csv"
        each_df = pd.read_csv(each_mge_file)
        each_df = each_df[each_df['final_score'] > cutoff]
        for idx, r in each_df.iterrows():
            
            taxa = bin2anno_dict[r['host']] if r['host'] in bin2anno_dict else "Unclassified"
            if re.search("Unclassified", taxa):# or re.search("s__", taxa):
                continue
            mge_hosts_dict[mge].add(taxa)
            mge_contigs_dict[mge].add(r['host'])
        if len(mge_hosts_dict[mge]) > 1:
            print (f"Multiple hosts for {mge}: {mge_hosts_dict[mge]} {mge_contigs_dict[mge]}")

    return mge_hosts_dict, mge_contigs_dict

def read_mge_type(mge_file):
    MGE_type_dict = {}
    with open(mge_file, "r") as f:
        f.readline()  # skip header
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            mge_id = fields[0]
            mge_type = fields[1]
            MGE_type_dict[mge_id] = mge_type
    return MGE_type_dict

if __name__ == "__main__":  
    all_dir = "/home/shuaiw/borg/paper/run2/"
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        # print (f"Processing {prefix}...")
        # work_dir = f"{all_dir}/{prefix}/{prefix}_methylation2"

        print(f"Processing {prefix}...")
        work_dir = f"/home/shuaiw/borg/paper/run2/{prefix}/"
        mge_file = os.path.join(work_dir, "all_mge.tsv")
        ## check if the mge_file exists
        if not os.path.exists(mge_file):
            print(f"File {mge_file} does not exist. Skipping {prefix}.")
            continue
        gtdk_bac_file = os.path.join(work_dir, "GTDB/gtdbtk.bac120.summary.tsv")
        gtdk_arc_file = os.path.join(work_dir, "GTDB/gtdbtk.ar122.summary.tsv")
        gtdk_all_file = os.path.join(work_dir, "GTDB/gtdbtk.all.summary.tsv")
        host_summary_file = os.path.join(work_dir, f"{prefix}_methylation3/host_summary.csv")  
        ## check if the host_summary_file exists
        if not os.path.exists(host_summary_file):
            print(f"File {host_summary_file} does not exist. Skipping {prefix}.")
            continue
        ## skip if the host_summary_file is empty, count the number of lines
        ## host sum lines
        line_num = sum(1 for line in open(host_summary_file) if line.strip())
        if line_num < 2:
            print(f"Host summary file {host_summary_file} is empty or has only header. Skipping {prefix}.")
            continue
        bin2anno_dict = read_gtdb(gtdk_all_file)
        MGE_type_dict = read_mge_type(mge_file)
        read_summary(host_summary_file, bin2anno_dict, MGE_type_dict)