from sklearn.metrics import normalized_mutual_info_score
import pandas as pd

# Example: True labels vs. Predicted clusters
true_labels = [0, 0, 1, 1, 2, 2]
predicted_clusters = [1, 1, 0, 0, 2, 2]


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

clster_out = "/home/shuaiw/borg/bench/zymo2/motif_cluster.csv"
df = pd.read_csv(clster_out)
answer_label = {}
for index, row in df.iterrows():
    contig = row["contigs"]
    answer_label[contig] = row['cluster']

true_labels = []
predicted_clusters = []
for contig in answer_label:
    true_labels.append(contig_index_dict[contig])
    predicted_clusters.append(answer_label[contig])





nmi_score = normalized_mutual_info_score(true_labels, predicted_clusters)
print(f"Normalized Mutual Information: {nmi_score:.4f}")
