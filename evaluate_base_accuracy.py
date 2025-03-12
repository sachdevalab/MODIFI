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
    for_loci_num = 0
    rev_loci_num = 0
    for_modified_num = 0
    rev_modified_num = 0

    motif_ipd_ratio = []

    for r, contig in REF.items():
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            tag = r + ":" + str(site+exact_pos) + "+"
            # print (tag)
            motif_loci_num += 1
            for_loci_num += 1
            if tag in ipd_ratio_dict:
                if ipd_ratio_dict[tag] == 1:
                    motif_modify_num += 1

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            motif_loci_num += 1
            rev_loci_num += 1
            if tag in ipd_ratio_dict:
                if ipd_ratio_dict[tag] == 1:
                    motif_modify_num += 1
    modify_ratio = motif_modify_num / motif_loci_num
    print ("motif loci num", motif_loci_num)
    print ("motif modify num", motif_modify_num)
    print ("modify ratio", modify_ratio)
    return modify_ratio



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
        if score <= score_cutoff:
            continue
        modified_loci[ref + ":" + str(pos) + strand] = score
    print ("no. of modified loci", len(modified_loci))
    return modified_loci

def read_ipd_ratio(csv_file, cutoff = 0.05):
    ipd_ratio_dict = {}
    # df = pd.read_csv(csv_file, nrows = 100000)
    df = pd.read_csv(csv_file)
    for index, row in df.iterrows():
        # print (row['refName'], row['tpl'], row['strand'], row['coverage'], row['ipd_ratio'])
        if row['strand'] == 1:
            strand_string = "-"
        else:
            strand_string = "+"
        tag = row['refName'] + ":" + str(row['tpl']+1) + strand_string
        if row['pvalue'] < cutoff:
            ipd_ratio_dict[tag] = 1
        else:
            ipd_ratio_dict[tag] = 0
    return ipd_ratio_dict

         
if __name__ == "__main__":
    motif_new = "CTGCAG"
    exact_pos = 5
    score_cutoff = 0
    # my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
    # gff = "/home/shuaiw/borg/bench/C227_native/gffs/CP011331.1.gff"

    my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
    native_csv = "/home/shuaiw/methylation/data/borg/bench/C227/native2/ipd_ratio/CP011331.1.ipd3.csv"
    control_csv = "/home/shuaiw/methylation/data/borg/bench/C227/WGA2/ipd_ratio/CP011331.1.ipd3.csv"
    # gff = "/home/shuaiw/methylation/data/borg/bench/zymo2/gffs/E_coli_K12-MG1655_1.gff"
    # all_motifs = "/home/shuaiw/methylation/data/borg/bench/zymo2/all.motifs.csv"
    # profile = "/home/shuaiw/borg/test.csv"

    # my_ref = sys.argv[1]
    # gff = sys.argv[2]
    # all_motifs = sys.argv[3]
    # profile = sys.argv[4]

    REF = read_ref(my_ref)
    # print (REF)
    # modified_loci = get_modified_ratio(gff)
    # motifs = pd.read_csv(all_motifs)
    native_ipd_ratio_dict = read_ipd_ratio(native_csv)

    recall = get_motif_sites(REF, motif_new, exact_pos, native_ipd_ratio_dict)
    control_ipd_ratio_dict = read_ipd_ratio(control_csv)
    FDR = get_motif_sites(REF, motif_new, exact_pos, control_ipd_ratio_dict)
    print ("recall", recall)
    print ("FDR", FDR)



