import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
from sample_object import get_unique_motifs, My_sample, get_ctg_taxa, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
from check_motif_change import given_species_drep


motif_list = [["GAA", 3],['YCTB',2],['GATC',2]]

all_dir = "/home/shuaiw/borg/paper/borg_data/methy2/"
# prefix="soil_s1_2"
prefix="soil_1"
contig = "SRVP18_trench_6_60cm_scaf_214_117_86_FINAL"



ctg_obj = My_contig(prefix, all_dir, contig)
modified_sites = ctg_obj.get_borg_mod_ratio(max_len=852769)
print (ctg_obj.modified_num, ctg_obj.modified_motif_num, ctg_obj.modified_ratio, ctg_obj.modified_motif_ratio, ctg_obj.motif_ratio)
ctg_obj.get_modified_loci()
ctg_obj.read_ref()
for motif_new, exact_pos in motif_list:
    motif_loci_num, motif_modify_num, ratio = ctg_obj.count_mod_frac_in_motif(motif_new, exact_pos, start=0, end=852769)
    print (motif_new, exact_pos, motif_modify_num, motif_loci_num, ratio)


# ctg_obj = My_contig(prefix, all_dir, contig)
# ctg_obj.read_ref()
# print (ctg_obj.ctg_ref)
# ctg_obj.get_ctg_len()
# print (ctg_obj.ctg_len)
# for motif_new, exact_pos in motif_list:
#     # ctg_obj.count_mod_dist(motif_new, exact_pos, start=0, end=852769)
#     ctg_obj.count_mod_dist(motif_new, exact_pos, start=852769, end=918293)