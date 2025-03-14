from sklearn.metrics import normalized_mutual_info_score
import pandas as pd

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
    clster_out = "/home/shuaiw/borg/bench/zymo2/motif_cluster.h.csv"
    answer_label = read_predicted_result(clster_out)
    contig_index_dict = read_zymo_truth()

    true_labels = []
    predicted_clusters = []
    for contig in answer_label:
        true_labels.append(contig_index_dict[contig])
        predicted_clusters.append(answer_label[contig])





    nmi_score = normalized_mutual_info_score(true_labels, predicted_clusters)
    print(f"Normalized Mutual Information: {nmi_score:.4f}")

def all_break():
    clster_out = "/home/shuaiw/borg/bench/all_break/motif_cluster.h.csv"
    answer_label = read_predicted_result(clster_out)
    contig_index_dict = read_all_break_truth()
    print (len(contig_index_dict), "truth contig number")
    print (len(answer_label), "predicted contig number")

    true_labels = []
    predicted_clusters = []
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

# for_zymo()
all_break()
