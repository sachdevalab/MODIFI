from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import adjusted_rand_score
import pandas as pd
import numpy as np  
from collections import defaultdict
import os

# Example: True labels vs. Predicted clusters
true_labels = [0, 0, 1, 1, 2, 2]
predicted_clusters = [1, 1, 0, 0, 2, 2]


def read_zymo_truth():
    fai = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged.fa.fai"
    spcies_index = {}
    index = 1

    contig_index_dict = {}
    for line in open(fai):
        contig = line.split("\t")[0]
        species_name = "_".join(contig.split("_")[:-1])
        if species_name not in spcies_index:
            spcies_index[species_name] = index
            index += 1
        contig_index_dict[contig] = spcies_index[species_name]
    return contig_index_dict

def read_all_break_truth():
    fai = "/home/shuaiw/methylation/data/borg/contigs/all_break.contigs.fa.fai"
    spcies_index = {}
    cluster_num_dict = defaultdict(int)
    index = 1

    contig_index_dict = {}
    for line in open(fai):
        contig = line.split("\t")[0]
        if contig[-1] == 'L' or contig[-1] == 'C':
            continue
        species_name = "_".join(contig.split("_")[:-2])
        if species_name not in spcies_index:
            spcies_index[species_name] = index
            index += 1
        contig_index_dict[contig] = spcies_index[species_name]
        cluster_num_dict[spcies_index[species_name]] += 1
    # print (cluster_num_dict)
    return contig_index_dict
    

def read_predicted_result(clster_out):
    df = pd.read_csv(clster_out)
    answer_label = {}
    for index, row in df.iterrows():
        # contig = row["motif"]
        contig = row["contigs"]
        answer_label[contig] = row['cluster']
    return answer_label

def read_metabat_result():


    # df = pd.read_csv("/home/shuaiw/borg/maxbat/zymo_bin.txt", header = None, sep = "\t")
    df = pd.read_csv("/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.metabat-bins--saveCls/bin", header = None, sep = "\t")
    # print (df)
    answer_label = {}
    for index, row in df.iterrows():
        # contig = row["motif"]
        contig = row[0]
        answer_label[contig] = row[1]
    return answer_label

def read_metabat_all_break():


    # df = pd.read_csv("/home/shuaiw/borg/maxbat/zymo_bin.txt", header = None, sep = "\t")
    df = pd.read_csv("/home/shuaiw/borg/contigs/all_break.contigs.fa.metabat-bins--saveCls/bin", header = None, sep = "\t")
    # print (df)
    answer_label = {}
    for index, row in df.iterrows():
        # contig = row["motif"]
        contig = row[0]
        answer_label[contig] = row[1]
    return answer_label

def for_zymo():
    # clster_out = "/home/shuaiw/borg/bench/zymo2/motif_cluster.h.csv"
    clster_out = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/motif_cluster.csv"
    # clster_out = "tmp/zymo.u.csv"
    answer_label = read_predicted_result(clster_out)
    contig_index_dict = read_zymo_truth()

    true_labels = []
    predicted_clusters = []
    for contig in answer_label:
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])
    print (len(true_labels))

    eva(true_labels, predicted_clusters)

def for_zymo_maxbat():

    answer_label = read_metabat_result()
    contig_index_dict = read_zymo_truth()

    true_labels = []
    predicted_clusters = []
    for contig in answer_label:
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])
    print (len(true_labels))

    eva(true_labels, predicted_clusters)

def compute_random_nmi(true_labels, num_iterations=1000, num_clusters=5):
    random_nmi_scores = []
    random_ari = []
    
    for _ in range(num_iterations):
        # Generate random cluster assignments
        random_clusters = np.random.randint(0, num_clusters, size=len(true_labels))
        # print (random_clusters)
        # Compute NMI between true labels and random clusters
        nmi_score = normalized_mutual_info_score(true_labels, random_clusters)
        ari = adjusted_rand_score(true_labels, random_clusters)
        random_nmi_scores.append(nmi_score)
        random_ari.append(ari)
    
    # Compute average NMI for random guessing
    return np.mean(random_nmi_scores), np.mean(random_ari)

def all_break():
    clster_out = "/home/shuaiw/borg/bench/all_break/motif_cluster.u.csv"
    clster_out = "tmp/all.u.csv"
    answer_label = read_predicted_result(clster_out)
    contig_index_dict = read_all_break_truth()
    print (len(contig_index_dict), "truth contig number")
    print (len(answer_label), "predicted contig number")

    true_labels = []
    predicted_clusters = []
    true_label_dict = defaultdict(int)
    for contig in answer_label:
        if contig not in contig_index_dict:
            continue
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])

    eva(true_labels, predicted_clusters)

def all_break_metabat():
    answer_label = read_metabat_all_break()
    contig_index_dict = read_all_break_truth()
    print (len(contig_index_dict), "truth contig number")
    print (len(answer_label), "predicted contig number")

    true_labels = []
    predicted_clusters = []
    true_label_dict = defaultdict(int)
    for contig in answer_label:
        if contig not in contig_index_dict:
            continue
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])

    eva(true_labels, predicted_clusters)

def read_all_binning():
    contigs = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    bin_anno =  "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.csv"
    f = open(bin_anno, 'w')
    ## use biopython to read the fasta
    from Bio import SeqIO
    import re
    for record in SeqIO.parse(contigs, "fasta"):
        contig_name = record.id
        # print (contig_name, record.description)
        field = str(record.description).split("bin=")
        bin_name = field[1][1:-1]
        # print (field[1], bin_name)
        print (contig_name, bin_name, sep = ",", file = f)
    f.close()

def get_binning_truth():
    bin_anno =  "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.csv"
    df = pd.read_csv(bin_anno, header = None)
    bin_index = 1
    bin_dict = {}
    for index, row in df.iterrows():
        # print (row)
        if row[1] == "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_UNK":
            continue
        if row[1] not in bin_dict:
            bin_dict[row[1]] = bin_index
            bin_index += 1
    contig_index_dict = {}
    for index, row in df.iterrows():
        # print (row)
        if row[1] == "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_UNK":
            continue
        contig_index_dict[row[0]] = bin_dict[row[1]]
    return contig_index_dict

def remove_single_cluster(true_labels, predicted_clusters):
    ## calculate the ## caluate the number of elements in each cluster
    cluster_num = defaultdict(int)
    for cluster in true_labels:
        cluster_num[cluster] += 1
    new_true, new_predict = [], []
    for i in range(len(true_labels)):
        if cluster_num[true_labels[i]] > 1:
            new_true.append(true_labels[i])
            new_predict.append(predicted_clusters[i])
    return new_true, new_predict 

def all_break2():
    # clster_out = "/home/shuaiw/borg/all_test_ccs3/motif_cluster.h.csv"
    # clster_out = "/home/shuaiw/borg/bench/all_subreads/motif_cluster.csv"
    # clster_out = "/home/shuaiw/borg/bench/s4/motif_cluster.h.csv"
    clster_out = "tmp/real.u.csv"
    answer_label = read_predicted_result(clster_out)
    contig_index_dict = get_binning_truth()
    print (len(contig_index_dict), "truth contig number")
    print (len(answer_label), "predicted contig number")

    true_labels = []
    predicted_clusters = []
    true_label_dict = defaultdict(int)
    for contig in answer_label:
        if contig not in contig_index_dict:
            continue
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])
    eva(true_labels, predicted_clusters)
   
def eva(true_labels, predicted_clusters):
    # true_labels, predicted_clusters  = remove_single_cluster(true_labels, predicted_clusters)
    # print (true_labels)
    print ("no of considered contigs:", len(true_labels))
    print (len(predicted_clusters))
    # print (contig_index_dict)
    # print (answer_label)

    nmi_score = normalized_mutual_info_score(true_labels, predicted_clusters)
    ARI = adjusted_rand_score(true_labels, predicted_clusters)
    print(f"Normalized Mutual Information: {nmi_score:.4f}", f"ARI: {ARI:.4f}")
    ## calculate random cluster nmi_score by shuffling the predicted_clusters
    # random_nmi_score, random_air = compute_random_nmi(true_labels, num_iterations=1000, num_clusters=max(true_labels))
    # print(f"Random Cluster NMI: {random_nmi_score:.4f}", f"random ARI: {random_air:.4f}")


def get_plasmid_dict(fai):
    index = 1
    plasmid_host_dict = {}
    host_contig_dict = defaultdict(list)

    contig_length_dict = {}
    for line in open(fai):
        contig = line.split("\t")[0]
        species_name = "_".join(contig.split("_")[:-1])
        length = int(line.split("\t")[1])
        index = int(contig.split("_")[-1])
        contig_length_dict[contig] = length
        if contig == "B_cepacia_UCB-717_4":
            continue
        if contig =="K_pneumoniae_BAA-2146_3":
            continue
        if contig in["S_aureus_USA300-TCH1516_3", "T_denticola_A_2"]:
            continue

        if species_name in ['V_parahaemolyticus_EB101', "B_cepacia_UCB-717", 'B_multivorans_249']:
            print (species_name)
            if index > 3:
                plasmid_host_dict[contig] = [species_name, length, []]
            else:
                host_contig_dict[species_name].append(contig)
        elif species_name in ['A_baumannii_AYE']:
            print (species_name)
            if index > 2:
                plasmid_host_dict[contig] = [species_name, length, []]
            else:
                host_contig_dict[species_name].append(contig)
        else:
            if index > 1:
                plasmid_host_dict[contig] = [species_name, length, []]
            else:
                host_contig_dict[species_name].append(contig)
    ## update host
    for contig in plasmid_host_dict:
        host = plasmid_host_dict[contig][0]
        plasmid_host_dict[contig][2] = host_contig_dict[host]
    # print (len(plasmid_host_dict), plasmid_host_dict)
    return plasmid_host_dict, contig_length_dict

def host_linkage_eva(): 
    fai = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai"
    plasmid_anno_file = fai + ".plasmid"
    plasmid_list = fai + ".plasmid.list"
    plasmid_host_dict, contig_length_dict = get_plasmid_dict(fai)
    f = open(plasmid_list, 'w')
    print ('seq_name', file = f)
    for contig in plasmid_host_dict:
        print (contig, file = f)
    f.close()
    output_host_linkage(plasmid_host_dict, contig_length_dict, plasmid_anno_file)


    # clster_out = "/home/shuaiw/borg/bench/zymo2/motif_cluster.h.csv"
    # clster_out = "/home/shuaiw/borg/bench/zymo6_NM200/motif_cluster.h.csv"
    # clster_out = "tmp/zymo.u.csv"
    # answer_label = read_predicted_result(clster_out)
    answer_label = read_metabat_result()
    plasmid_cluster_id = {}
    for contig in plasmid_host_dict:
        if contig in answer_label:
            plasmid_cluster_id[contig] = answer_label[contig]
        else:
            plasmid_cluster_id[contig] = -1
    cluster_dict = defaultdict(list)
    for contig in answer_label:
        cluster_dict[answer_label[contig]].append(contig)
    plasmid_cluster_dict = {}
    for plasmid in plasmid_cluster_id:
        plasmid_cluster_dict[plasmid] = cluster_dict[plasmid_cluster_id[plasmid]]

    recall = 0
    FP = 0
    report_num = 0
    total = 0
    for plasmid in plasmid_host_dict:
        host = plasmid_host_dict[plasmid]
        cluster = plasmid_cluster_dict[plasmid]
        valid_host = check_host(host[2], cluster, plasmid)
        if len(cluster) > 0:
            report_num += 1
        print (plasmid, host[2], cluster, valid_host, report_num)
        # print (valid_host)
        if valid_host == "link host" or valid_host == "FP":
            recall += 1
        if valid_host == "FP":
            FP += 1
        total += 1
    print (recall, total, "recall", recall/total, report_num, "FP", FP/report_num)

def cal_AUC(): 
    fai = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai"
    plasmid_anno_file = fai + ".plasmid"
    plasmid_list = fai + ".plasmid.list"
    plasmid_host_dict, contig_length_dict = get_plasmid_dict(fai)
    dir = "/home/shuaiw/borg/bench/zymo_new_ref_NM3/hosts/"
    dir = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30_filter/hosts/"
    # dir = "/home/shuaiw/borg/bench/zymo_new_ref_p0.1_cov1_s30/hosts/"
    cutoff = 0.4
    data = []
    for cutoff in range(0, 100, 5):
        cutoff = cutoff / 100
        guess_num = 0
        correct_num = 0
        
        for plasmid in plasmid_host_dict:
            host_file = dir + plasmid + ".host_prediction.csv"

            if not os.path.exists(host_file):
                print (host_file, "not exist")
                continue
            df = pd.read_csv(host_file)
            ## get first row of df
            for index, row in df.iterrows():
                best_host = row['host']
                best_score = row['final_score']
                break
            if best_score > cutoff:
                guess_num += 1
                if best_host in plasmid_host_dict[plasmid][2]:
                    correct_num += 1
        if guess_num == 0:
            continue
        recall = correct_num / len(plasmid_host_dict)
        precision = correct_num / guess_num
        print (recall, precision)
        data.append([cutoff, recall, precision, 1-precision])
    ## plot the AUC curve
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import seaborn as sns
    
    df = pd.DataFrame(data, columns=['cutoff', 'recall', 'precision', 'FP'])
    print (df)
    plt.figure(figsize=(8, 6))
    ## use seaborn
    sns.set(style="whitegrid")
    ## plot dot plot, with x as precision, y as recall
    ## sort df by recall
    df = df.sort_values(by=['cutoff'])
    plt.plot(df['FP'], df['recall'], marker='o', linestyle='-', color='b')
    plt.xlabel('FP')
    plt.ylabel('recall')
    plt.savefig("../tmp/precision_recall_curve.png")




def output_host_linkage(plasmid_host_dict, contig_length_dict, plasmid_anno_file):
    f = open(plasmid_anno_file, 'w')
    for plasmid in plasmid_host_dict:
        plasmid_length = contig_length_dict[plasmid]
        print (plasmid, plasmid_length, ",".join(plasmid_host_dict[plasmid][2]), sep="\t", file = f)
    f.close()


def check_host(host_list, cluster, plasmid):
    has_host_contig = False
    has_other_host = False
    host_species = "_".join(plasmid.split("_")[:-1])
    for contig in cluster:
        if contig in host_list:
            has_host_contig = True
        species_name = "_".join(contig.split("_")[:-1])
        if species_name != host_species:
            has_other_host = True
    flag = "no host"
    if has_host_contig:
        flag = "link host"
        if has_other_host:
            flag = "FP"
    return flag
        
    
    
    # if has_host_contig and not has_other_host:
    #     return "correct linkage"
    # elif has_host_contig and has_other_host:
    #     return "wrong linkage"
    # elif not has_host_contig and not has_other_host:
    #     return "no linkage"
    # elif not has_host_contig and has_other_host:
    #     print ("###########other")
    #     return "wrong linkage"



cal_AUC()
# host_linkage_eva()
# for_zymo()
# for_zymo_maxbat()
# all_break()
# all_break_metabat()

# all_break2()
