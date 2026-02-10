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

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'linkage'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))


from sample_object import get_unique_motifs, My_sample, Isolation_sample



def compare_linkage(our_linkages, spacer_linkage_dict, sample_obj):
    consistent_num = 0
    mge_type_consistent_dict = defaultdict(int)
    for mge in our_linkages:
        if mge in spacer_linkage_dict:
            ## check if there is intersection
            mge_type = sample_obj.mge_dict.get(mge, "unknown")
            if set(our_linkages[mge]) & set(spacer_linkage_dict[mge]):
                print(f"Intersection found for {mge}")
                consistent_num += 1
                mge_type_consistent_dict[mge_type] += 1
    spacer_linkage_num = len(spacer_linkage_dict)
    our_linkage_num = len(our_linkages)
    consistent_rate = consistent_num / our_linkage_num if our_linkage_num > 0 else 0
    print(f"Consistent linkages: {consistent_rate} (spacer: {spacer_linkage_num}, ours: {our_linkage_num})")

    return consistent_num, spacer_linkage_num, our_linkage_num, mge_type_consistent_dict

def plot_bar(df, fig_dir):
    # Melt for grouped barplot
    df_melted = df.melt(id_vars=["Sample"], var_name="Linkage Type", value_name="Count")

    # Define colors matching the provided colormap
    colors = {"Our Linkages": "#E57373", "Spacer Linkages": "#4DD0E1", "Consistent Linkages": "#4DD0E1"}
    
    fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    # First subplot: grouped barplot
    sns.barplot(data=df_melted, x="Sample", y="Count", hue="Linkage Type", ax=axs[0], palette=colors)
    axs[0].set_title("Linkage Comparison", fontsize=14)
    axs[0].set_yscale('log')
    axs[0].tick_params(axis='x', rotation=90, labelsize=8)
    axs[0].tick_params(axis='y', labelsize=8)

    # Second subplot: only consistent linkages
    sns.barplot(data=df, x="Sample", y="Consistent Linkages", ax=axs[1], color="#4DD0E1")
    axs[1].set_title("Consistent Linkages Only", fontsize=14)
    axs[1].tick_params(axis='x', rotation=90, labelsize=8)
    axs[1].tick_params(axis='y', labelsize=8)

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/linkage_comparison_spacer.pdf")

def plot_mge_type_bar(df, fig_dir):
    ## plot stacked barplot of consistent linkages by MGE type
    df_pivot = df.pivot(index='Sample', columns='MGE Type', values='Consistent Linkages')
    # Use custom colors based on the provided colormap
    ax = df_pivot.plot(kind='bar', stacked=True, figsize=(6, 6), color=['#F8766D', '#00BFC4'])
    plt.ylabel("Number of Consistent Linkages", fontsize=10)
    plt.xlabel("Sample", fontsize=10)
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=10)
    plt.gca().yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/mge_type_consistent_linkages.pdf")

if __name__ == "__main__":
    data = []
    all_dir = "/home/shuaiw/borg/paper/run2/"
    fig_dir = "../../tmp/figures/link_accuracy/"
    mge_type_consistent_dict_all = defaultdict(int)
    mge_type_df = []
    
    # Load linkages from CSV file
    print("Loading linkages from CSV file...")
    csv_file = "/home/shuaiw/mGlu/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
    linkage_df = pd.read_csv(csv_file)
    print(f"Loaded {len(linkage_df)} linkages from {csv_file}")
    
    # Group by sample
    for prefix, sample_group in linkage_df.groupby('sample'):
        print(f"Processing {prefix}...")
        sample_obj = My_sample(prefix, all_dir)
        
        # Read MGE types
        sample_obj.read_MGE()
        
        # Reconstruct our_linkages dictionary from CSV data
        our_linkages = defaultdict(list)
        for _, row in sample_group.iterrows():
            mge = row['MGE']
            host = row['host']
            our_linkages[mge].append(host)
        
        # Read spacer linkages for comparison
        spacer_linkage_dict = sample_obj.read_spacer(mismatch_allowed=0)

        consistent_num, spacer_linkage_num, our_linkage_num, mge_type_consistent_dict = compare_linkage(our_linkages, spacer_linkage_dict, sample_obj)
        mge_type_consistent_dict_all = {k: mge_type_consistent_dict_all.get(k, 0) + mge_type_consistent_dict.get(k, 0) for k in set(mge_type_consistent_dict_all) | set(mge_type_consistent_dict)}

        data.append([prefix, consistent_num, spacer_linkage_num, our_linkage_num])
        for mge_type in mge_type_consistent_dict:
            mge_type_df.append([prefix, mge_type, mge_type_consistent_dict[mge_type]])

    print ("Consistent linkages by MGE type:")
    for mge_type in mge_type_consistent_dict_all:
        print (f"{mge_type}\t{mge_type_consistent_dict_all[mge_type]}")
    df = pd.DataFrame(data, columns=["Sample", "Consistent Linkages", "Spacer Linkages", "Our Linkages"])
    df = df.sort_values(by='Sample')
    mge_type_df = pd.DataFrame(mge_type_df, columns=["Sample", "MGE Type", "Consistent Linkages"])
    mge_type_df = mge_type_df.sort_values(by=['Sample', 'MGE Type'])
    
    plot_bar(df, fig_dir)
    plot_mge_type_bar(mge_type_df, fig_dir)
    ## count the proportion of our linkages that are consistent in all
    ## sum all consistent linkages and our linkages
    total_consistent = df["Consistent Linkages"].sum()
    total_our = df["Our Linkages"].sum()
    overall_consistency = total_consistent / total_our if total_our > 0 else 0
    print(f"Overall consistency across all samples: {overall_consistency} ({total_consistent} consistent out of {total_our} our linkages)")




