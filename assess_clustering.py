from sklearn.metrics import normalized_mutual_info_score
import pandas as pd
import numpy as np  
from collections import defaultdict

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
    print (cluster_num_dict)
    return contig_index_dict
    

def read_predicted_result(clster_out):


    df = pd.read_csv(clster_out)
    answer_label = {}
    for index, row in df.iterrows():
        # contig = row["motif"]
        contig = row["contigs"]
        answer_label[contig] = row['cluster']
    return answer_label


def for_zymo():
    clster_out = "/home/shuaiw/borg/bench/zymo2/motif_cluster.csv"
    answer_label = read_predicted_result(clster_out)
    contig_index_dict = read_zymo_truth()

    true_labels = []
    predicted_clusters = []
    for contig in answer_label:
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])
    print (len(true_labels))





    nmi_score = normalized_mutual_info_score(true_labels, predicted_clusters)
    print(f"Normalized Mutual Information: {nmi_score:.4f}")

def compute_random_nmi(true_labels, num_iterations=1000, num_clusters=5):
    random_nmi_scores = []
    
    for _ in range(num_iterations):
        # Generate random cluster assignments
        random_clusters = np.random.randint(0, num_clusters, size=len(true_labels))
        # print (random_clusters)
        # Compute NMI between true labels and random clusters
        nmi_score = normalized_mutual_info_score(true_labels, random_clusters)
        random_nmi_scores.append(nmi_score)
    
    # Compute average NMI for random guessing
    return np.mean(random_nmi_scores), np.std(random_nmi_scores)

def all_break():
    clster_out = "/home/shuaiw/borg/bench/all_break/motif_cluster.csv"
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

    print (len(true_labels))
    print (len(predicted_clusters))
    # print (contig_index_dict)
    # print (answer_label)

    nmi_score = normalized_mutual_info_score(true_labels, predicted_clusters)
    print(f"Normalized Mutual Information: {nmi_score:.4f}")
    ## calculate random cluster nmi_score by shuffling the predicted_clusters
    random_nmi_score, random_nmi_std = compute_random_nmi(true_labels, num_iterations=1000, num_clusters=max(true_labels))
    print(f"Random Cluster NMI: {random_nmi_score:.4f}")

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

    

# for_zymo()
# all_break()
read_all_binning()
