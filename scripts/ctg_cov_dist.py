import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas




def get_cov_dist(dp_file):
    df = pandas.read_csv(dp_file, sep = "\t")
    ## to float 

    df['totalAvgDepth'] = df['totalAvgDepth'].astype(float)
    covs = df['totalAvgDepth'].values
    ## plot the distribution using seaborn
    sns.set(style="whitegrid")
    sns.histplot(covs, bins=200, kde=True)
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.title("SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META")
    
    ## set xlim to 100, it did not work
    plt.xlim(0, 100)
    print ("Coverage distribution plot saved.")
    plt.savefig("tmp/cov_dist.pdf")


dp_file = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.depth.txt"

get_cov_dist(dp_file)
