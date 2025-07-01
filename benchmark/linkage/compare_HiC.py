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
from matplotlib_venn import venn2

def read_hic(individual_hic_file, sample, bin_depth_dict):
    filtered_linkage = {}
    df1 = pd.read_csv(hic_file)  ##  linkage exsits in >= 10 samples
    for index, row in df1.iterrows():
        tag = row['circular_element'] + "&" + row['bin']
        if tag not in filtered_linkage:
            filtered_linkage[tag] = 1
    hic_linkages = defaultdict(list)
    df = pd.read_csv(individual_hic_file, sep = "\t")
    for index, row in df.iterrows():
        if row['metagenome'] == sample:
            tag = row['circular_element'] + "&" + row['bin']
            if tag not in filtered_linkage:
                continue
            if row['bin'] not in bin_depth_dict:
                continue
            if bin_depth_dict[row['bin']] < 10:
                continue
            if row['circular_element'] not in bin_depth_dict:
                continue
            if bin_depth_dict[row['circular_element']] < 10:
                continue
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
        ## try plasmid is in header of row
        try:
            plasmid_name = row['plasmid']
        ## otherwise index with MGE
        except KeyError:
            plasmid_name = row['MGE']
        our_linkages[plasmid_name].append(bin_name)
        our_ctg_linkages[plasmid_name] = row['host']
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

def compare_hic_our(hic_linkages, our_linkages, our_ctg_linkages, bin2ctg_dict, bin_depth_dict,ctg_len_dict):
    both_link = 0
    cosistency_num = 0
    incosistency_num = 0
    both_linkage_dict = {}
    inconsistent_hic_set = set()
    hic_bin_set = {"consistent":[], "inconsistent":[]}
    consistent_hic_set = set()
    for plasmid in hic_linkages:
        if plasmid not in our_linkages:
            print (f"###  Hi-C {plasmid} {bin_depth_dict[plasmid]} : {hic_linkages[plasmid]} {bin_depth_dict[hic_linkages[plasmid][0]]} is not in our linkages")
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
            both_linkage_dict[plasmid] = our_linkages[plasmid]
            consistent_hic_set.add(bin2ctg_dict[hic_linkages[plasmid][0]][0])
            hic_bin_set["consistent"].append(hic_linkages[plasmid][0])
        else:
            print (f"{plasmid} is not consistent: {our_linkages[plasmid]} ({our_ctg_linkages[plasmid]}) vs Hi-C {hic_linkages[plasmid]} : ({bin2ctg_dict[hic_linkages[plasmid][0]][0]}, {len(bin2ctg_dict[hic_linkages[plasmid][0]])})\n")
            inconsistent_hic_set.add(bin2ctg_dict[hic_linkages[plasmid][0]][0])
            hic_bin_set["inconsistent"].append(hic_linkages[plasmid][0])
            incosistency_num += 1
    if both_link == 0:
        cosistency_rate = 0
    else:
        cosistency_rate = cosistency_num / both_link
    print(f"both linkages: {both_link}")
    print(f"cosistency linkages: {cosistency_num}")
    print (f"incosistency linkages: {incosistency_num}")
    print (f"cosistency rate: {cosistency_rate}")
    print ("our link num", len(our_linkages))
    print (len(inconsistent_hic_set), "inconsistent Hi-C contigs:", inconsistent_hic_set)
    print (len(consistent_hic_set), "consistent Hi-C contigs:", consistent_hic_set)

    for bin_name in set(hic_bin_set["consistent"]):
        if bin_name not in bin2ctg_dict:
            print (f"bin {bin_name} is not in bin2ctg_dict")
            continue
        len_list = []
        for ctg in bin2ctg_dict[bin_name]:
            if ctg not in ctg_len_dict:
                print (f"contig {ctg} is not in ctg_len_dict")
                continue
            len_list.append(ctg_len_dict[ctg])
        print ("consistent", bin_name, max(len_list), len(len_list), np.mean(len_list))

    for bin_name in set(hic_bin_set["inconsistent"]):
        len_list = []
        for ctg in bin2ctg_dict[bin_name]:
            len_list.append(ctg_len_dict[ctg])
        print ("inconsistent", bin_name, max(len_list), len(len_list), np.mean(len_list))

    return both_link, cosistency_num, cosistency_rate, len(our_linkages), both_linkage_dict

def main():
    my_cutoff = 0.5
    ctg_len_dict = get_ctg_len(fai)
    ctg2bin_dict, bin2ctg_dict = load_ctg2bin(ctg2bin)
    bin_depth_dict = read_depth(depth_file)
    # hic_linkages = read_hic(hic_file)
    hic_linkages = read_hic(individual_hic_file, sample, bin_depth_dict)
    data = []
    our_linkages, our_ctg_linkages = read_our(host_sum, ctg2bin_dict, my_cutoff)
    # our_linkages, our_ctg_linkages = read_our_multiple(host_sum, ctg2bin_dict, my_cutoff)
    print (f"cutoff: {my_cutoff}")
    both_link, cosistency_num, consis_rate, our_num, both_linkage_dict = compare_hic_our(hic_linkages, our_linkages, our_ctg_linkages, bin2ctg_dict, bin_depth_dict,ctg_len_dict)
    
    df2 = pd.read_csv(host_sum)
    ## add a new column for both
    df2['both_link'] = df2['MGE'].apply(lambda x: 1 if x in both_linkage_dict else 0)
    ## output df2 to a csv file
    df2.to_csv(host_sum_compare, index=False)
    data.append([my_cutoff, both_link, cosistency_num, consis_rate, our_num])
    print (len(hic_linkages), len(our_linkages), len(both_linkage_dict))
    venn2(subsets = (len(hic_linkages)-len(both_linkage_dict), len(our_linkages)-len(both_linkage_dict), len(both_linkage_dict)), set_labels = ('Hi-C', 'Methylation'))
    # plt.show()
    ## save the venn diagram
    plt.savefig(f"../../tmp/results/venn_{my_cutoff}.png")
    df = pd.DataFrame(data, columns=['cutoff', 'both_link', 'cosistency_num', 'consis_rate', 'our_num'])

    sns.set(style="whitegrid")
    plt.figure(figsize=(5, 9))  # Adjust the figure size to accommodate the third subplot

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
    ctg_len_dict = {}
    out = open(bin_file, "w")
    for line in open(fai, "r"):
        if line.startswith("#"):
            continue
        contig, length = line.strip().split("\t")[:2]
        ctg_len_dict[contig] = int(length)
        if contig not in ctg2bin_dict:
            bin_name = contig
        else:
            bin_name = ctg2bin_dict[contig]
        print(f"{contig}\t{bin_name}", file=out)
    out.close()
    print(f"bin file is saved to {bin_file}")
    return  ctg_len_dict

def get_ctg_len(fai):

    ctg_len_dict = {}
    for line in open(fai, "r"):
        if line.startswith("#"):
            continue
        contig, length = line.strip().split("\t")[:2]
        ctg_len_dict[contig] = int(length)
    return  ctg_len_dict

def read_depth(depth_file):
    depth_dict = {}
    with open(depth_file, "r") as f:
        for line in f:
            if line.startswith("contig"):
                continue
            contig, depth = line.strip().split(",")
            depth_dict[contig] = float(depth)
    ## get bin depth
    ctg2bin_dict = {}
    with open(bin_file, "r") as f:
        for line in f:
            contig, bin_name = line.strip().split("\t")
            ctg2bin_dict[contig] = bin_name
    bin_depth_dict = {}
    for contig, depth in depth_dict.items():
        bin_name = ctg2bin_dict.get(contig, contig)
        if bin_name not in bin_depth_dict:
            bin_depth_dict[bin_name] = []
        bin_depth_dict[bin_name].append(depth)
    for bin_name, depths in bin_depth_dict.items():
        bin_depth_dict[bin_name] = np.mean(depths)
    return bin_depth_dict

if __name__ == "__main__":
    hic_file = "/home/shuaiw/borg/pengfan/hic_10mgs_linkages.csv"
    individual_hic_file = "/home/shuaiw/borg/pengfan/bins_circular_elements_filt.tsv"
    bin_file = "/home/shuaiw/borg/pengfan/10mgs_bins.tab"
    bin_dir = "/groups/diamond/projects/animal/rumen/RuReacBro20203/analysis/RuReacBro2023_CircCont/bins/Final_Genomes_qc_rmcirc/sequence"
    ctg2bin = "/home/shuaiw/borg/pengfan/ctg2bin.csv"
    fai = "/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa.fai"
    # host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_Comb_RF_LR_bin/host_summary.csv"
    host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/host_summary.csv"
    depth_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/mean_depth.csv"
    # host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_12_72h_200ppm_r2_HMW_LR_bin/host_summary.csv"
    # host_sum =  "/home/shuaiw/methylation/data/borg/pengfan/total_summary.csv"
    host_sum_compare =  "/home/shuaiw/methylation/data/borg/pengfan/total_summary_compare.csv"
    sample = "RuReacBro_20230708_11_72h_200ppm_r1"

    # host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_26_72h_NC_r4_LR_bin/host_summary.csv"
    # host_sum = "/home/shuaiw/methylation/data/borg/pengfan/total_summary.csv"

    main()
    # generate_bin_file(bin_file, fai)
    # contig2bin(bin_dir, ctg2bin)

# metagenome
# RuReacBro_20230708_11_24h_200ppm_r1
# RuReacBro_20230708_11_40h_200ppm_r1
# RuReacBro_20230708_11_56h_200ppm_r1
# RuReacBro_20230708_11_72h_200ppm_r1
# RuReacBro_20230708_11_8h_200ppm_r1
# RuReacBro_20230708_12_24h_200ppm_r2
# RuReacBro_20230708_12_40h_200ppm_r2
# RuReacBro_20230708_12_56h_200ppm_r2
# RuReacBro_20230708_12_72h_200ppm_r2
# RuReacBro_20230708_12_8h_200ppm_r2
# RuReacBro_20230708_14_24h_200ppm_r3
# RuReacBro_20230708_14_40h_200ppm_r3
# RuReacBro_20230708_14_56h_200ppm_r3
# RuReacBro_20230708_14_72h_200ppm_r3
# RuReacBro_20230708_14_8h_200ppm_r3
# RuReacBro_20230708_19_24h_NC_r3
# RuReacBro_20230708_19_40h_NC_r3
# RuReacBro_20230708_19_56h_NC_r3
# RuReacBro_20230708_19_72h_NC_r3
# RuReacBro_20230708_19_8h_NC_r3
# RuReacBro_20230708_8_24h_NC_r1
# RuReacBro_20230708_8_40h_NC_r1
# RuReacBro_20230708_8_56h_NC_r1
# RuReacBro_20230708_8_72h_NC_r1
# RuReacBro_20230708_8_8h_NC_r1
# RuReacBro_20230708_9_24h_NC_r2
# RuReacBro_20230708_9_40h_NC_r2
# RuReacBro_20230708_9_56h_NC_r2
# RuReacBro_20230708_9_72h_NC_r2
# RuReacBro_20230708_9_8h_NC_r2
# RuReacBro_20230708_Comb_RF_WB
# RuReacBro_20230708_Cow1_RF
# RuReacBro_20230708_Cow2_RF
# RuReacBro_20230708_Cow3_RF