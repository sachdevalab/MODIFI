import os
import sys
import pandas as pd

def collect_motifs(folder, MIN_FRAC, MIN_detect):
    data = []
    for file in os.listdir(folder):
        if file.endswith(".csv"):
            motif = pd.read_csv(os.path.join(folder, file), sep = ",")
            ## add a new column for the motif indentifier, which is the combination of motif and pos
            motif = motif[(motif['fraction'] > MIN_FRAC) & (motif['nDetected'] > MIN_detect)]
            # print ("no of motifs", len(motif))
            motif['indentifier'] = motif['motifString'] + "_" + motif['centerPos'].astype(str)
            data += list(motif['fraction'].values)
        # if len(data) > 100:
        #     break
    # print (data)
    ## plot the distribution of the data using seaborn
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.histplot(data, bins=100, kde=True)
    plt.xlabel('Fraction of methylated sites')
    plt.ylabel('Density')
    plt.title('Distribution of Methylation Fraction')
    plt.savefig(os.path.join("../../tmp/results", "methylation_fraction_distribution_cow.png"))




if __name__ == "__main__":
    MIN_FRAC = 0
    MIN_detect = 1000
    # folder = "/home/shuaiw/borg/all_test/motifs"
    # all_motif = "/home/shuaiw/borg/all_test/test_motifs.csv"

    # folder = "/home/shuaiw/borg/bench/soil/run1/motifs"
    folder = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/motifs"
    collect_motifs(folder, MIN_FRAC, MIN_detect)