"""
Given the folder motifs/, enumerage all the motifs and collect them into a single file.
"""

import os
import sys
import pandas as pd

def collect_motifs(folder, all_motif, MIN_FRAC, MIN_detect):
    motifs = []
    for file in os.listdir(folder):
        if file.endswith(".csv"):
            motif = pd.read_csv(os.path.join(folder, file), sep = ",")
            ## add a new column for the motif indentifier, which is the combination of motif and pos
            motif = motif[(motif['fraction'] > MIN_FRAC) & (motif['nDetected'] > MIN_detect)]
            # print ("no of motifs", len(motif))
            motif['indentifier'] = motif['motifString'] + "_" + motif['centerPos'].astype(str)
            motifs.append(motif)
    motifs = pd.concat(motifs)
    ## only keep one row for each motif indentifier
    motifs = motifs.drop_duplicates(subset=['indentifier'])
    motifs.to_csv(all_motif, index=False)
    print ("no of motifs", len(motifs))


if __name__ == "__main__":
    MIN_FRAC = 0.5
    MIN_detect = 100
    # folder = "/home/shuaiw/borg/all_test/motifs"
    # all_motif = "/home/shuaiw/borg/all_test/test_motifs.csv"

    folder = sys.argv[1]
    all_motif = sys.argv[2]
    collect_motifs(folder, all_motif, MIN_FRAC, MIN_detect)