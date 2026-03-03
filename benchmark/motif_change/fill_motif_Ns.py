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
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
from scipy.stats import fisher_exact
# ...existing code...
from scipy.stats import fisher_exact

def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos, max_len = 1000000000):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1


    record_modified_sites = {}
    record_modified_times = {}
    i = 0
    for r, contig in REF.items():
        contig = contig[:max_len]
        for site in nt_search(str(contig), motif_new)[1:]:
            ## priny the upstream and downstream 10 bp sequences
            upstream_site = max(0, site-20)
            downstream_site = site+20
            # print (motif_new, exact_pos, contig[site:site+motif_len], contig[upstream_site:downstream_site])
            record_modified_sites[contig[site:site+motif_len]] = [motif_new, exact_pos, contig[site:site+motif_len], contig[upstream_site:downstream_site]]
            record_modified_times[contig[site:site+motif_len]] = record_modified_times.get(contig[site:site+motif_len], 0) + 1
            i += 1
            # if i > 5:
            #     break

        # for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
        #     1:
        # ]:
        #     # tag = r + ":" + str(site+rev_exact_pos) + "-"
        #     # record_modified_sites[tag] = motif_new

        #     upstream_site = max(0, site-10)
        #     downstream_site = site+10
        #     # print (motif_new, exact_pos, contig[upstream_site:downstream_site])
    ##  sort the record_modified_sites by the times
    record_modified_sites = dict(sorted(record_modified_sites.items(), key=lambda item: item[1], reverse=True))
    ## print top 5
    for motif, info in list(record_modified_sites.items())[:5]:
        print(motif, record_modified_times[motif], info)


ref = "/home/shuaiw/borg/paper/E_faecalis/infant_14_31_C.fa"
REF = read_ref(ref)
motifs_list = [["TCAYNNNNNNTTG",3], ["CAANNNNNNRTGA", 3], ["CRTANNNNNNRTG", 4], ["CAYNNNNNNTAYG", 2]]
for motif, exact_pos in motifs_list:
    get_motif_sites(REF, motif, exact_pos)