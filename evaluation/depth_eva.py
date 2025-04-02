import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score



def compare_motif():
    df = pd.read_csv(bench_profile)
    df2 = pd.read_csv(depth_file)
    print (df)
    ## collect index names
    high_motifs = df['motifString'].tolist()
    low_motifs = df2['motifString'].tolist()
    print ("High motifs: ", high_motifs)
    print ("Low motifs: ", low_motifs)
    ## use high motifs as reference, calculate recall and precision
    recall_num  = 0
    precision_num = 0
    for motif in high_motifs:
        if motif in low_motifs:
            recall_num += 1
        else:
            print ("Motif not in low motifs: ", motif)
    ## calculate recall
    recall = recall_num / len(high_motifs)
    ## calculate precision
    for motif in low_motifs:
        if motif in high_motifs:
            precision_num += 1
    precision = precision_num / len(low_motifs)
    ## calculate F1 score
    f1 = 2 * (precision * recall) / (precision + recall)
    print ("Recall: ", recall)
    print ("Precision: ", precision)
    print ("F1 score: ", f1)

def read_gff(gff):
    modified_loci = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            locus_tag = cols[0] + "_" + str(cols[3]) + "_" + cols[6]
            score = int(cols[5])
            if score < 30:
                continue
            modified_loci[locus_tag] = score
    return modified_loci

def compare_gff():
    bench_loci = read_gff(bench_gff)
    depth_loci = read_gff(depth_gff)
    ## use bench loci as reference, calculate recall and precision
    recall_num  = 0
    precision_num = 0
    for locus in bench_loci:
        if locus in depth_loci:
            recall_num += 1
        # else:
        #     print ("Locus not in depth loci: ", locus)
    ## calculate recall
    recall = recall_num / len(bench_loci)
    ## calculate precision
    for locus in depth_loci:
        if locus in bench_loci:
            precision_num += 1
    precision = precision_num / len(depth_loci)
    ## calculate F1 score
    f1 = 2 * (precision * recall) / (precision + recall)
    print ("Recall: ", recall)
    print ("Precision: ", precision)
    print ("F1 score: ", f1)

bench_gff = "/home/shuaiw/borg/bench/zymo_new_ref_NM3/gffs/E_coli_H10407_2.gff"
depth_gff = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/gffs/E_coli_H10407_2.gff" 
compare_gff()





# bench_profile = "/home/shuaiw/borg/bench/zymo_new_ref_NM3/motif_profile.csv"
# depth_file = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/motif_profile.csv"
# # depth_file = "/home/shuaiw/borg/bench/zymo_new_ref/motif_profile.csv"
# compare_motif()
