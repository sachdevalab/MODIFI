import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns

def find_linkage():
    col1 = 'SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_33175_L'
    datafile = "/home/shuaiw/methylation/data/borg/bench/all_subreads/motif_profile2.csv"
    # datafile = "/home/shuaiw/methylation/data/borg/all_test_ccs2/motif_profile.csv"

    df = pd.read_csv(datafile, index_col=0)
    df = df.loc[:, (df > 0.1).any(axis=0)]
    # print (df.shape)

    ## get all_column names
    euclidean_distance_dict = {}
    all_columns = df.columns
    if col1 not in all_columns:
        raise ValueError("column not in the dataframe")
    # print (all_columns)
    for col in all_columns:
        if col == col1:
            continue
        euclidean_distance = ((df[col1] - df[col]) ** 2).sum() ** 0.5
        euclidean_distance_dict[col] = euclidean_distance
    ## sort the dictionary by value
    euclidean_distance_tuple = sorted(euclidean_distance_dict.items(), key=lambda item: item[1])
    ## print top ten columns
    for i in range(10):
        print (euclidean_distance_tuple[i])
    ## filter the df with the top ten columns
    top_ten_columns = [item[0] for item in euclidean_distance_tuple[:10]] + [col1]
    df_top_ten = df[top_ten_columns]
    df = df_top_ten.loc[(df_top_ten > 0.2).any(axis=1)]
    print (df)
    df = df.T
    ## plot the heatmap using seaborn
    g= sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(10, 10), cbar_kws={'label': 'Intensity'}, yticklabels=1)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=5)  # Y-axis labels
    ## x-axis labels
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=5, rotation=90)
    plt.savefig("tmp/test_heatmap.pdf")


find_linkage()

