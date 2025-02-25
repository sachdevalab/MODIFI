import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def visu(df, figure_path):
    sns.set(style="whitegrid")
    fig, axs = plt.subplots(2, 2, figsize=(20, 10))
    ## first row is covergae, second is tMean, third is control, fourth is ipd_ratio
    ## plot for each strand separately
    sns.histplot(df, x="coverage", hue="strand", multiple="stack", ax=axs[0, 0])
    sns.histplot(df, x="tMean", hue="strand", multiple="stack", ax=axs[0, 1])
    sns.histplot(df, x="control", hue="strand", multiple="stack", ax=axs[1, 0])
    sns.histplot(df, x="ipd_ratio", hue="strand", multiple="stack", ax=axs[1, 1])
    ## save the plot

    ### cal the number with strand 1 and 0
    # print (df[df['strand'] == 1].shape[0])
    # print (df[df['strand'] == 0].shape[0])
    plt.savefig(figure_path)


df = pd.read_csv("/home/shuaiw/borg/bench/break/ipd_ratio/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_961_C_459827_919654.ipd3.csv")
figure_path = "tmp/visu.png"
visu(df, figure_path)