# -*- coding: utf-8 -*-
import pandas as pd
import os
import sys
from collections import defaultdict
import numpy as np
### plot the consistency rate line plot
## plot another subplot with the our link num
import matplotlib.pyplot as plt
import seaborn as sns

def read_hic(hic_file):
    hic_linkages = defaultdict(list)
    df = pd.read_csv(hic_file)
    for index, row in df.iterrows():
        hic_linkages[row['circular_element']].append(row['bin'])
    multiple_host_plasmid_num = 0
    for plasmid in hic_linkages:
        if len(hic_linkages[plasmid]) > 1:
            multiple_host_plasmid_num += 1
            # print (f"{plasmid} has multiple host: {hic_linkages[plasmid]}")
    print (f"multiple host plasmid num: {multiple_host_plasmid_num} out of {len(hic_linkages)}")
    return hic_linkages

def read_our(host_sum, ctg2bin_dict, score_cutoff = 0.45):
    df = pd.read_csv(host_sum)
    df = df[df['final_score'] > score_cutoff]
    our_linkages = defaultdict(list)
    our_ctg_linkages = {}
    for index, row in df.iterrows():
        if row['host'] not in ctg2bin_dict:
            bin_name = row['host']
            ## raise error
            # print (f"contig {row['host']} is not in ctg2bin_dict")
            # sys.exit(1)
        else:
            bin_name = ctg2bin_dict[row['host']]
        our_linkages[row['plasmid']].append(bin_name)
        our_ctg_linkages[row['plasmid']] = row['host']
        # our_linkages[row['plasmid']] = row['host']
    multiple_host_plasmid_num = 0
    for plasmid in our_linkages:
        if len(our_linkages[plasmid]) > 1:
            multiple_host_plasmid_num += 1
            # print (f"{plasmid} has multiple host: {our_linkages[plasmid]}")
    print (f"multiple host plasmid num: {multiple_host_plasmid_num} out of {len(our_linkages)}")
    return our_linkages, our_ctg_linkages

def read_our_multiple(host_sum, ctg2bin_dict, score_cutoff = 0.45):
    df = pd.read_csv(host_sum)
    df = df[df['final_score'] > score_cutoff]
    our_linkages = defaultdict(list)
    our_ctg_linkages = {}
    main_dir = "/".join(host_sum.split("/")[:-1])
    for index, row in df.iterrows():
        plasmid_host_file = os.path.join(main_dir, "hosts", row['plasmid'] + ".host_prediction.csv")
        ## read the host prediction file
        if not os.path.exists(plasmid_host_file):
            print (f"{plasmid_host_file} does not exist")
            continue
        plasmid_host_df = pd.read_csv(plasmid_host_file)
        plasmid_host_df = plasmid_host_df[plasmid_host_df['final_score'] > score_cutoff]
        for index1, row1 in plasmid_host_df.iterrows():
            if row1['host'] not in ctg2bin_dict:
                bin_name = row1['host']
                ## raise error
                # print (f"contig {row1['host']} is not in ctg2bin_dict")
                # sys.exit(1)
            else:
                bin_name = ctg2bin_dict[row1['host']]
            our_linkages[row['plasmid']].append(bin_name)
        our_ctg_linkages[row['plasmid']] = row['host']
        # our_linkages[row['plasmid']] = row['host']
    multiple_host_plasmid_num = 0
    for plasmid in our_linkages:
        if len(our_linkages[plasmid]) > 1:
            multiple_host_plasmid_num += 1
            # print (f"{plasmid} has multiple host: {our_linkages[plasmid]}")
    print (f"multiple host plasmid num: {multiple_host_plasmid_num} out of {len(our_linkages)}")
    return our_linkages, our_ctg_linkages

def compare_hic_our(hic_linkages, our_linkages, our_ctg_linkages, bin2ctg_dict):
    both_link = 0
    cosistency_num = 0
    for plasmid in our_linkages:
        if plasmid not in hic_linkages:
            # print(f"{plasmid} is not in HiC linkages")
            continue
        both_link += 1
        # if our_linkages[plasmid] in hic_linkages[plasmid]:
        #     cosistency_num += 1
        ## check if the two list has a same element
        if set(our_linkages[plasmid]) & set(hic_linkages[plasmid]):
            cosistency_num += 1
        else:
            print (f"{plasmid} is not consistent: {our_linkages[plasmid]} ({our_ctg_linkages[plasmid]}) vs Hi-C {hic_linkages[plasmid]} : ({bin2ctg_dict[hic_linkages[plasmid][0]][0]}, {len(bin2ctg_dict[hic_linkages[plasmid][0]])})")
    if both_link == 0:
        cosistency_rate = 0
    else:
        cosistency_rate = cosistency_num / both_link
    print(f"both linkages: {both_link}")
    print(f"cosistency linkages: {cosistency_num}")
    print (f"cosistency rate: {cosistency_rate}")
    print ("our link num", len(our_linkages))
    return both_link, cosistency_num, cosistency_rate, len(our_linkages)

def main():
    ctg2bin_dict, bin2ctg_dict = load_ctg2bin(ctg2bin)
    hic_linkages = read_hic(hic_file)
    ## output the MGE with multiple host

    # """
    data = []
    for cutoff in range(6, 18):
        my_cutoff = cutoff / 20
        
        our_linkages, our_ctg_linkages = read_our(host_sum, ctg2bin_dict, my_cutoff)
        # our_linkages, our_ctg_linkages = read_our_multiple(host_sum, ctg2bin_dict, my_cutoff)
        print (f"cutoff: {my_cutoff}")
        both_link, cosistency_num, consis_rate, our_num = compare_hic_our(hic_linkages, our_linkages, our_ctg_linkages, bin2ctg_dict)
        data.append([my_cutoff, both_link, cosistency_num, consis_rate, our_num])
    df = pd.DataFrame(data, columns=['cutoff', 'both_link', 'cosistency_num', 'consis_rate', 'our_num'])

    sns.set(style="whitegrid")
    plt.figure(figsize=(5, 9))  # Adjust the figure size to accommodate the third subplot

    ## plot three subplots
    plt.subplot(3, 1, 1)
    plt.plot(df['cutoff'], df['consis_rate'], marker='o')
    plt.title('Ratio of Consistency Linkages')
    plt.xlabel('Cutoff')
    plt.ylabel('Ratio of Consistency Linkages')

    plt.subplot(3, 1, 2)
    plt.plot(df['cutoff'], df['cosistency_num'], marker='o')
    plt.title('Number of Consistency Linkages')
    plt.xlabel('Cutoff')
    plt.ylabel('Number of Consistency Linkages')

    plt.subplot(3, 1, 3)
    plt.plot(df['cutoff'], df['our_num'], marker='o', color='green')  # Plot our_num
    plt.title('Number of Our Linkages')
    plt.xlabel('Cutoff')
    plt.ylabel('Number of Our Linkages')

    plt.tight_layout()
    plt.savefig('../tmp/both_linkages.png')
    # """

def contig2bin(bin_dir, ctg2bin):
    out = open(ctg2bin, "w")
    out.write("contig,bin\n")
    for bin_file in os.listdir(bin_dir):
        if bin_file.endswith(".fa"):
            bin_name = bin_file[:-3]
            with open(os.path.join(bin_dir, bin_file), "r") as f:
                for line in f:
                    if line.startswith(">"):
                        contig_name = line.strip().split(" ")[0][1:]
                        print(f"{contig_name},{bin_name}", file=out)
    out.close()
    print(f"contig2bin file is saved to {ctg2bin}")

def load_ctg2bin(ctg2bin):
    ctg2bin_dict = {}
    with open(ctg2bin, "r") as f:
        for line in f:
            if line.startswith("contig"):
                continue
            contig, bin_name = line.strip().split(",")
            ctg2bin_dict[contig] = bin_name
    bin2ctg_dict = {}
    for contig, bin_name in ctg2bin_dict.items():
        if bin_name not in bin2ctg_dict:
            bin2ctg_dict[bin_name] = []
        bin2ctg_dict[bin_name].append(contig)
    return ctg2bin_dict, bin2ctg_dict

def generate_bin_file(bin_file, fai):
    ctg2bin_dict, bin2ctg_dict = load_ctg2bin(ctg2bin)
    out = open(bin_file, "w")
    for line in open(fai, "r"):
        if line.startswith("#"):
            continue
        contig, length = line.strip().split("\t")[:2]
        if contig not in ctg2bin_dict:
            bin_name = contig
        else:
            bin_name = ctg2bin_dict[contig]
        print(f"{contig}\t{bin_name}", file=out)
    out.close()
    print(f"bin file is saved to {bin_file}")


if __name__ == "__main__":
    hic_file = "/home/shuaiw/borg/pengfan/hic_10mgs_linkages.csv"
    bin_file = "/home/shuaiw/borg/pengfan/10mgs_bins.tab"
    bin_dir = "/groups/diamond/projects/animal/rumen/RuReacBro20203/analysis/RuReacBro2023_CircCont/bins/Final_Genomes_qc_rmcirc/sequence"
    ctg2bin = "/home/shuaiw/borg/pengfan/ctg2bin.csv"
    fai = "/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa.fai"
    # host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_Comb_RF_LR_bin/host_summary.csv"
    host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/host_summary.csv"

    # host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_26_72h_NC_r4_LR_bin/host_summary.csv"
    # host_sum = "/home/shuaiw/methylation/data/borg/pengfan/total_summary.csv"

    main()
    # generate_bin_file(bin_file, fai)
    # contig2bin(bin_dir, ctg2bin)