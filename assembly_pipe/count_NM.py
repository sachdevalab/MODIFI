import pysam
import os
from multiprocessing import Pool
import sys
import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def handle_each_contig(bam, fig, prefix, NM_file, q=20):

    #  Create a new header that only includes the specific contig
    samfile = pysam.AlignmentFile(bam, "rb")
    f = open(NM_file, 'w')
    NM_list = []
    for read in samfile:
        if read.is_unmapped:
            continue
        if read.mapping_quality < q:
            continue
        if not read.has_tag('NM'):
            # print (f"Read {read.query_name} has no NM tag")
            continue
        NM= read.get_tag("NM")
        NM_list.append(NM)
        print(read.query_name, NM, file=f)
    f.close()
    samfile.close()
    plot(NM_list, fig, prefix)

def plot(NM_list, fig, prefix):
    ## calculate mean and median values if NM
    if NM_list:
        mean_NM = round(np.mean(NM_list),2)
        median_NM = round(np.median(NM_list),2)
        ## count the quatile value
        q1_NM = round(np.percentile(NM_list, 25),2)
        q3_NM = round(np.percentile(NM_list, 75),2)
        print(f"{prefix} Mean NM: {mean_NM}, Median NM: {median_NM}, Q1 NM: {q1_NM}, Q3 NM: {q3_NM}")
    else:
        print("No reads with NM tag found.")
        mean_NM, median_NM, q1_NM, q3_NM = 0, 0, 0, 0

    

    plt.figure(figsize=(10, 6))
    sns.histplot(NM_list, bins=30, kde=False)
    plt.xlabel('NM tag value')
    plt.ylabel('Read count')
    plt.title(f'mean: {mean_NM}, median: {median_NM}, Q1: {q1_NM}, Q3: {q3_NM}')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(fig)
    plt.close()




# prefix = "ERR12723529_mice"


all_dir = "/home/shuaiw/borg/paper/run2/"

for my_dir in os.listdir(all_dir):
    prefix = my_dir
    print (prefix)
    bam = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.bam"
    fig = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.NM.pdf"
    NM_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.NM.csv"

    ## check if bam exists
    if os.path.exists(bam):
        print(f"Processing {bam}")
        handle_each_contig(bam, fig, prefix,NM_file)