from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys


def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos, ipd_ratio_dict):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1
    motif_sites = {}
    motif_loci_num = 0
    motif_modify_num = 0

    for r, contig in REF.items():
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            tag = r + ":" + str(site+exact_pos) + "+"
            # print (tag)
            
            if tag in ipd_ratio_dict:
                motif_loci_num += 1
                if ipd_ratio_dict[tag] == 1:
                    motif_modify_num += 1

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            if tag in ipd_ratio_dict:
                motif_loci_num += 1
                if ipd_ratio_dict[tag] == 1:
                    motif_modify_num += 1
    modify_ratio = motif_modify_num / motif_loci_num
    print ("motif loci num", motif_loci_num)
    print ("motif modify num", motif_modify_num)
    print ("modify ratio", modify_ratio)
    return modify_ratio, motif_loci_num, motif_modify_num



def get_modified_ratio(gff):
    ## read the gff file
    f = open(gff, "r")
    modified_loci = {}
    for line in f:
        if line[0] == "#":
            continue
        line = line.strip().split("\t")

        ref = line[0]
        pos = int(line[3]) 
        strand = line[6]
        score = int(line[5])
        # if score <= score_cutoff:
        #     continue
        modified_loci[ref + ":" + str(pos) + strand] = score
    print ("no. of modified loci", len(modified_loci))
    return modified_loci

def read_ipd_ratio(csv_file, cutoff = 0.05, coverage = 10):
    """Deprecated - use read_ipd_ratio_from_df instead"""
    df = pd.read_csv(csv_file)
    return read_ipd_ratio_from_df(df, cutoff, coverage)

def read_ipd_ratio_from_df(df, cutoff = 0.05, coverage = 10):
    """Optimized version that works with pre-loaded DataFrame"""
    # Filter by coverage
    df = df[df['coverage'] >= coverage].copy()
    
    # Vectorized operations - much faster than iterrows()
    df['strand_string'] = df['strand'].map({1: "-", 0: "+"})
    df['tag'] = df['refName'] + ":" + (df['tpl'] + 1).astype(str) + df['strand_string']
    df['modified'] = (df['pvalue'] < cutoff).astype(int)
    
    # Convert to dictionary in one operation
    ipd_ratio_dict = df.set_index('tag')['modified'].to_dict()
    return ipd_ratio_dict

         
if __name__ == "__main__":
    motif_new = "CTGCAG"
    exact_pos = 5
    score_cutoff = 0
    # my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
    # gff = "/home/shuaiw/borg/bench/C227_native/gffs/CP011331.1.gff"

    my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"

    # native_csv = "/home/shuaiw/methylation/data/borg/bench/C227/native2/ipd_ratio/CP011331.1.ipd3.csv"
    # control_csv = "/home/shuaiw/methylation/data/borg/bench/C227/WGA2/ipd_ratio/CP011331.1.ipd3.csv"

    # native_csv = "/home/shuaiw/borg/allison/ecoli/native/ipd_ratio/CP011331.1.ipd3.csv"
    # control_csv = "/home/shuaiw/borg/allison/ecoli/WGA/ipd_ratio/CP011331.1.ipd3.csv"

    # gff = "/home/shuaiw/methylation/data/borg/bench/zymo2/gffs/E_coli_K12-MG1655_1.gff"
    # all_motifs = "/home/shuaiw/methylation/data/borg/bench/zymo2/all.motifs.csv"
    # profile = "/home/shuaiw/borg/test.csv"


    # native_csv = "/home/shuaiw/borg/paper/base/pure/native/ipd_ratio/CP011331.1.ipd3.csv"
    # control_csv = "/home/shuaiw/borg/paper/base/pure/control/ipd_ratio/CP011331.1.ipd3.csv"
    # output_csv = "../../tmp/figures/base_benchmark/ecoli_base_pure.csv"

    native_csv = "/home/shuaiw/borg/paper/base/meta/native/ipd_ratio/CP011331.1.ipd3.csv"
    control_csv = "/home/shuaiw/borg/paper/base/meta/control/ipd_ratio/CP011331.1.ipd3.csv"
    output_csv = "../../tmp/figures/base_benchmark/ecoli_base_meta.csv"

    REF = read_ref(my_ref)
    
    # Read CSV files once instead of 56 times
    print("Loading CSV files...")
    native_df = pd.read_csv(native_csv)
    control_df = pd.read_csv(control_csv)
    print("CSV files loaded.")
    
    data = []
    for coverage in [5, 10, 15, 20]:
        for cutoff in [0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2]:
            print ("coverage", coverage, "cutoff", cutoff)
            
            # Use pre-loaded DataFrames instead of reading files repeatedly
            native_ipd_ratio_dict = read_ipd_ratio_from_df(native_df, cutoff, coverage)
            recall, native_motif_loci_num, native_motif_modify_num = get_motif_sites(REF, motif_new, exact_pos, native_ipd_ratio_dict)
            
            control_ipd_ratio_dict = read_ipd_ratio_from_df(control_df, cutoff, coverage)
            FDR, motif_loci_num, motif_modify_num = get_motif_sites(REF, motif_new, exact_pos, control_ipd_ratio_dict)
            
            print ("recall", recall)
            print ("FDR", FDR)
            data.append([coverage, cutoff, recall, FDR, native_motif_loci_num, native_motif_modify_num , motif_loci_num, motif_modify_num])
    df = pd.DataFrame(data, columns=['coverage', 'cutoff', 'recall', 'FDR', 'native_motif_loci_num', 'native_motif_modify_num', 'motif_loci_num', 'motif_modify_num'])
    ## save df
    df.to_csv(output_csv, index = False)




