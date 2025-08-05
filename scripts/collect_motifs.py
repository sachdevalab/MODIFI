"""
Given the folder motifs/, enumerage all the motifs and collect them into a single file.
"""

import os
import sys
import pandas as pd
from Bio.Seq import Seq  # BioPython is used for reverse complement functionality


from derep_motifs import MotifFilter, uniq_similar_motifs

def collect_motifs(folder, all_motif, MIN_FRAC, MIN_detect):
    motifs = []
    for file in os.listdir(folder):
        if file.endswith(".csv"):
            motif = pd.read_csv(os.path.join(folder, file), sep = ",")
            ## add a new column for the motif indentifier, which is the combination of motif and pos
            motif = motif[(motif['fraction'] > MIN_FRAC) & (motif['nDetected'] > MIN_detect)]
            # print ("no of motifs", len(motif))
            if len(motif) == 0:
                continue
            motif['indentifier'] = motif['motifString'] + "_" + motif['centerPos'].astype(str)
            motifs.append(motif)
    if len(motifs) == 0:
        motifs = pd.DataFrame(columns=['motifString', 'centerPos', 'modificationType', 'fraction', 'nDetected', 'nGenome', 'indentifier'])
    else:
        motifs = pd.concat(motifs)
        ## only keep one row for each motif indentifier
        ## add the nDetected together for the same motif indentifier
        motifs = motifs.groupby('indentifier').agg({
            'motifString': 'first',
            'centerPos': 'first',
            'modificationType': 'first',
            'fraction': 'first',
            'nDetected': 'sum',
            "nGenome": 'sum',
        }).reset_index()
        ## fraction  == nDetected / nGenome
        motifs['fraction'] = round(motifs['nDetected'] / motifs['nGenome'],2)
    # motifs = motifs.drop_duplicates(subset=['indentifier'])
    motifs.to_csv(all_motif, index=False)
    print ("no of motifs from all contigs:", len(motifs))

    # drep_motifs(motifs, all_motif)
    print ("skip drep motifs")
    drep_motif_file = all_motif.replace(".csv", "_drep.csv")
    motifs.to_csv(drep_motif_file, index=False)




def drep_motifs(motifs, all_motif):
    drep_motif_file = all_motif.replace(".csv", "_drep.csv")
    motif_data = []
    for index, row in motifs.iterrows():
        motif_data.append({
            'motif': row['motifString'],
            'centerPos': row['centerPos'],
            'host_meth': row['nDetected'],
            'host_total': row['nGenome'],
            'indentifier': row['indentifier'],
        })
    motif_filter = MotifFilter(motif_data)
    motif_data,final_similarity_groups = motif_filter.filter()
    # print (len(motif_data))
    ## filter the motifs with identifier exsits in motif_data
    retained_motifs = {}
    for motif_obj in motif_data:
        retained_motifs[motif_obj['indentifier']] = 1
    motifs = motifs[motifs['indentifier'].isin(retained_motifs.keys())]
    print (len(motifs), "motifs after dereplication")
    motifs.to_csv(drep_motif_file, index=False)
    ## print the final similarity groups
    print("Final similarity groups:")
    for group in final_similarity_groups:
        print(group)

def test_drep_motifs():
    motifs = pd.read_csv(all_motif)
    drep_motifs(motifs, "/home/shuaiw/borg/paper/run2/cow_1/cow_1_methylation2/test.all.motifs.drep.csv")


if __name__ == "__main__":
    MIN_FRAC = 0.5
    MIN_detect = 100
    # folder = "/home/shuaiw/borg/all_test/motifs"
    # all_motif = "/home/shuaiw/borg/all_test/test_motifs.csv"

    folder = sys.argv[1]
    all_motif = sys.argv[2]
    MIN_FRAC = float(sys.argv[3])
    MIN_detect = int(sys.argv[4])
    collect_motifs(folder, all_motif, MIN_FRAC, MIN_detect)
    # test_drep_motifs()