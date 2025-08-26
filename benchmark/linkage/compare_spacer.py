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
from scipy.stats import ttest_ind

from compare_bin3c import read_our

def read_spacer(spacer_linkage_file):
    spacer_linkage_dict = defaultdict(set)
    ## check if spacer_linkage_file empty
    if not os.path.exists(spacer_linkage_file) or os.path.getsize(spacer_linkage_file) == 0:
        return spacer_linkage_dict
    df = pd.read_csv(spacer_linkage_file, sep="\t")
    # print (df)
    for index, row in df.iterrows():
        spacer_linkage_dict[row['sseqid']].add(row['qseqid_ctg'])
    return spacer_linkage_dict

def compare_linkage(our_linkages, spacer_linkage_dict):
    consistent_num = 0
    for mge in our_linkages:
        if mge in spacer_linkage_dict:
            ## check if there is intersection
            if set(our_linkages[mge]) & set(spacer_linkage_dict[mge]):
                print(f"Intersection found for {mge}")
                consistent_num += 1
    spacer_linkage_num = len(spacer_linkage_dict)
    our_linkage_num = len(our_linkages)
    consistent_rate = consistent_num / our_linkage_num if our_linkage_num > 0 else 0
    print(f"Consistent linkages: {consistent_rate} (spacer: {spacer_linkage_num}, ours: {our_linkage_num})")
    return consistent_num, spacer_linkage_num, our_linkage_num


if __name__ == "__main__":
    data = []
    all_dir = "/home/shuaiw/borg/paper/run2/"
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        print (f"Processing {prefix}...")
        # prefix="cow_bioreactor_1"
        # prefix="ERR12723528_mice"
        bin3c_cluster = f"{all_dir}/{prefix}/hic/bin3c_clust/clustering.mcl"
        host_sum = f"{all_dir}/{prefix}/{prefix}_methylation2/host_summary.csv"
        contact_value_file = f"{all_dir}/{prefix}/hic/bin3c/contact_values.txt"
        spacer_linkage_file = f"{all_dir}/{prefix}/spacer/{prefix}_spacer_linkage.tsv"

        if not os.path.exists(host_sum):
            print(f"Skipping {prefix} as all_host_file does not exist.")
            continue


        our_linkages, our_ctg_linkages = read_our(host_sum)
        spacer_linkage_dict = read_spacer(spacer_linkage_file)
        consistent_num, spacer_linkage_num, our_linkage_num = compare_linkage(our_linkages, spacer_linkage_dict)
        data.append([prefix, consistent_num, spacer_linkage_num, our_linkage_num])

    df = pd.DataFrame(data, columns=["Sample", "Consistent Linkages", "Spacer Linkages", "Our Linkages"])
    df = df.sort_values(by='Sample')

    # Melt for grouped barplot
    df_melted = df.melt(id_vars=["Sample"], var_name="Linkage Type", value_name="Count")

    fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    # First subplot: grouped barplot
    sns.barplot(data=df_melted, x="Sample", y="Count", hue="Linkage Type", ax=axs[0])
    axs[0].set_title("Linkage Comparison", fontsize=14)
    axs[0].set_yscale('log')
    axs[0].tick_params(axis='x', rotation=90, labelsize=8)
    axs[0].tick_params(axis='y', labelsize=8)

    # Second subplot: only consistent linkages
    sns.barplot(data=df, x="Sample", y="Consistent Linkages", ax=axs[1], color="tab:blue")
    axs[1].set_title("Consistent Linkages Only", fontsize=14)
    axs[1].tick_params(axis='x', rotation=90, labelsize=8)
    axs[1].tick_params(axis='y', labelsize=8)

    plt.tight_layout()
    plt.savefig(f"../../tmp/results2/linkage_comparison_spacer.pdf")
