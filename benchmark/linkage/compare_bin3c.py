
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


sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))


from sample_object import get_unique_motifs, My_sample, Isolation_sample



def assess_linage(bin3c_cluster_list, sample_obj):
    
    our_linkages, our_ctg_linkages,linkage_info_list = sample_obj.read_linkage_dict()
    consistent_num = 0
    inconsistent_num = 0
    for mge in our_linkages:
        ## check if mge is in bin3c_cluster_list
        found = False
        mge_bin3c_cluster = []
        for cluster in bin3c_cluster_list:
            if mge in cluster:
                found = True
                mge_bin3c_cluster = cluster
                break
        if not found:
            # print (f"{mge} not found in bin3c_cluster_list")
            continue
        if len(mge_bin3c_cluster) < 2:
            # print (f"{mge}'s bin has only one contig in bin3c_cluster_list: {mge_bin3c_cluster}")
            continue

        cosist_flag = False
        for bin_name in our_linkages[mge]:
            if bin_name in mge_bin3c_cluster:
                cosist_flag = True
                break
        if cosist_flag:
            consistent_num += 1
            # print (f"{mge} is consistent with bin3c_cluster_list")
        else:
            inconsistent_num += 1
            print (f"{mge} is inconsistent with bin3c_cluster_list, {our_linkages[mge]}")
            print (
                f"bin3c_cluster_list: {mge_bin3c_cluster}")

    print (f"consistent num: {consistent_num}, inconsistent num: {inconsistent_num}", f"rate: {consistent_num / (consistent_num + inconsistent_num) if (consistent_num + inconsistent_num) > 0 else 0}")
    return our_linkages


def read_contact_values(contact_value_file):
    contact_values = {}
    with open(contact_value_file, 'r') as f:
        f.readline()  # Skip header line
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            seq1, seq2, count = parts
            count = int(count)
            # Store both directions for easy lookup
            if seq1 not in contact_values:
                contact_values[seq1] = {}
            contact_values[seq1][seq2] = count
            if seq2 not in contact_values:
                contact_values[seq2] = {}
            contact_values[seq2][seq1] = count
    return contact_values

def plot_hic(fig_dir, df):
    plt.figure(figsize=(6, 6))
    ax = sns.boxplot(x="Sample", y="Contact Value", hue="Type", data=df, palette=["#FF6347", "#4682B4"])
    plt.yscale('log')
    plt.ylabel("Contact Values (log scale)", fontsize=12)
    plt.legend(title="Type", fontsize=12, title_fontsize=12)
    plt.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{fig_dir}/bin3c_boxplot.png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

def compare_contact(mge, bin_name, contact_values):
    skip = False
    if mge not in contact_values:
        print(f"{mge} not found in contact values")
        skip = True
        return 0, "", 0, [], skip
    ## sort contact_values[mge] by value, and return a tuple
    sorted_contact_values = sorted(contact_values[mge].items(), key=lambda x: x[1], reverse=True)
    ## if the mge only has contact value with itself, skip
    over_self_only = True
    total_contacts = 0 ## skip it self
    for a in contact_values[mge]:
        if a != mge and contact_values[mge][a] > 0:
            over_self_only = False
            total_contacts += contact_values[mge][a]
            # break
    if over_self_only:
        print(f">>>>>>>>>>>>>>>>{mge} only has contact value with itself, skipping")
        print (sorted_contact_values[:6], "\n\n")
        skip = True
        return 0, "", 0, [], skip
    if total_contacts < 20:
        print(f">>>>>>>>>>>>>>>>{mge} has very low total contact value {total_contacts}, skipping")
        print (sorted_contact_values[:6], "\n\n")
        skip = True
        return 0, "", 0, [], skip
    # Find the bin_name with highest contact value for contact_values[mge], excluding itself
    max_contact = 0
    max_bin_name = mge
    for a in contact_values[mge]:
        if a == mge:
            continue
        if contact_values[mge][a] >= max_contact:
            max_contact = contact_values[mge][a]
            max_bin_name = a


    if bin_name not in contact_values[mge]:
        contact_values[mge][bin_name] = 0
        print(f"******************{bin_name} not found in contact values for {mge}")
    return contact_values[mge][bin_name], max_bin_name, max_contact, sorted_contact_values, skip

def get_info_ctg(seq1):
    over_self_only = True
    total_contacts = 0 ## skip it self
    for a in contact_values[seq1]:
        if a != seq1 and contact_values[seq1][a] > 0:
            over_self_only = False
            total_contacts += contact_values[seq1][a]
    if over_self_only or total_contacts < 20:
        return False
    return True

def random_test(contact_values, iterations=1000):
    random_contact_values = []
    random_num = 0
    while random_num < iterations:
        seq1 = np.random.choice(list(contact_values.keys()))
        seq2 = np.random.choice(list(contact_values.keys()))
        if seq1 == seq2:
            continue

        if not get_info_ctg(seq1):
            continue


        if seq2 not in contact_values[seq1]:
            random_contact_values.append(0)
        else:
            random_contact_values.append(contact_values[seq1][seq2])
            # print (f"random contact value between {seq1} and {seq2}: {contact_values[seq1][seq2]}")
        random_num += 1
    return random_contact_values

if __name__ == "__main__":
    # prefix="cow_bioreactor_1"
    # prefix="cow_1"
    data = []
    all_dir = "/home/shuaiw/borg/paper/run2/"
    fig_dir = "../../tmp/figures/link_accuracy/"
    # for prefix in ["cow_bioreactor_1", "cow_bioreactor_4", "cow_bioreactor_5"]:
    for prefix in ["cow_bioreactor_1"]:
        sample_obj = My_sample(prefix, all_dir, methy_v=4)
        bin3c_cluster_list = sample_obj.read_bin3c()
        our_linkages = assess_linage(bin3c_cluster_list, sample_obj)
        contact_values = read_contact_values(sample_obj.contact_value_file)
        
        ## get the contact values for our linkages
        our_linkage_contact = []
        for mge in our_linkages:
            bin_name = our_linkages[mge][0]
            contact_value, max_bin_name, max_contact, sorted_contact_values, skip = compare_contact(mge, bin_name, contact_values)
            if skip:
                continue
            # print (contact_values[mge], max_bin_name)
            our_linkage_contact.append(contact_values[mge][bin_name])
            print (f"{mge} contact value with {bin_name}: {contact_values[mge][bin_name]}, highest contact value {max_bin_name} of {contact_values[mge][max_bin_name]}")
            print (sorted_contact_values[:6], "\n\n")

        print (f"our linkages contact values: {our_linkage_contact}")
        ## also collect the contact values for 1000 random seq pairs
        random_contact_values = random_test(contact_values, iterations=1000)

        # print (f"random contact values: {random_contact_values}")
        print (f"average contact value for our linkages: {np.mean(our_linkage_contact) if len(our_linkage_contact) > 0 else 0}")
        print (f"Median contact value for our linkages: {np.median(our_linkage_contact) if len(our_linkage_contact) > 0 else 0}")
        ## cont how many value are > 0, proportion
        positive_count = sum(1 for value in our_linkage_contact if value > 0)
        print(f"Proportion of positive contact values for our linkages: {positive_count}/{len(our_linkage_contact)} = {positive_count / len(our_linkage_contact) if len(our_linkage_contact) > 0 else 0}")
        print (f"average contact value for random pairs: {np.mean(random_contact_values) if len(random_contact_values) > 0 else 0}")
        for value in our_linkage_contact:
            data.append([prefix, "Our Linkages", value])
        for value in random_contact_values:
            data.append([prefix, "Random Pairs", value])

    df = pd.DataFrame(data, columns=["Sample", "Type", "Contact Value"])
    plot_hic(fig_dir, df)




