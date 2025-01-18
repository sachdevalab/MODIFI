import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict
import numpy as np

def read_log(log):
    f = open(log)
    f.readline()
    f.readline()
    f.readline()
    predicted_list = []
    real_list = []
    ipdsum_list = []
    for line in f:
        field = line.strip().split() 
        if len(field) < 6:
            continue
        predicted = float(field[3])
        real = float(field[4])
        ipdsum = float(field[5])
        predicted_list.append(predicted)
        real_list.append(real)
        ipdsum_list.append(ipdsum)
        if len(predicted_list) > 100000:
            break
    f.close()

    ### cal prearson correlation between predicted and real, and between predicted and ipdsum
    from scipy.stats import pearsonr
    print (pearsonr(predicted_list, real_list))
    print (pearsonr(predicted_list, ipdsum_list))

    plot(predicted_list, real_list)

def plot(predicted_list, real_list):
    import matplotlib.pyplot as plt
    plt.scatter(predicted_list, real_list)
    ## svae the plot
    plt.savefig("tmp/scatter.pdf")

def read_log_csv(csv):

    df = pd.read_csv(csv, nrows=100000000)
    # df = pd.read_csv(csv)
    print ("length of df", len(df))
    ## cal correlation between estimated column and real column
    corr = pearsonr(df['estimated'], df['real'])
    print ("mean correlation", corr)
    corr = pearsonr(df['median'], df['real'])
    print ("median correlation", corr)
    ## cal correlation between estimated column and ipdsum column
    corr = pearsonr(df['ipd_sum'], df['real'])
    print ("ipdsummary correlation", corr)

    ### CAL MAR between real and estimated, and real and ipdsum
    mae = np.mean(np.abs(df['real'] - df['estimated']))
    mae_ipd_sum = np.mean(np.abs(df['real'] - df['ipd_sum']))
    print ("MAE between real and estimated", mae)
    print ("MAE between real and ipdsum", mae_ipd_sum)

def count_log_csv(csv):

    # df = pd.read_csv(csv, nrows=1000000)
    df = pd.read_csv(csv)

    corr = pearsonr(df['estimated'], df['real'])
    print ("mean correlation", corr)
    corr = pearsonr(df['median'], df['real'])
    print ("median correlation", corr)
    ## cal correlation between estimated column and ipdsum column
    corr = pearsonr(df['ipd_sum'], df['real'])
    print ("ipdsummary correlation", corr)

    ### CAL MAR between real and estimated, and real and ipdsum
    mae = np.mean(np.abs(df['real'] - df['estimated']))
    mae_ipd_sum = np.mean(np.abs(df['real'] - df['ipd_sum']))
    print ("MAE between real and estimated", mae)
    print ("MAE between real and ipdsum", mae_ipd_sum)

    estimate_dict = defaultdict(list)
    ipdsum_dict = defaultdict(list)
    real_dict = defaultdict(list)
    ## enumerate the row
    for index, row in df.iterrows():
        # print(index, row['kmer'])
        estimate_dict[row['kmer']].append(row['estimated'])
        ipdsum_dict[row['kmer']].append(row['ipd_sum'])
        real_dict[row['kmer']].append(row['real'])
    ## sort the kmer with the length of value list
    sorted_estimate = sorted(estimate_dict.items(), key=lambda x: len(x[1]), reverse=True)
    for i in range(10):
        kmer = sorted_estimate[i][0]
        ## round two for each element in estimate_dict[kmer]
        estimate_dict[kmer] = [round(x, 2) for x in estimate_dict[kmer]]
        ## cal MAE between real and estimated
        mae = np.mean(np.abs(np.array(real_dict[kmer]) - np.array(estimate_dict[kmer])))
        mae_ipd_sum = np.mean(np.abs(np.array(real_dict[kmer]) - np.array(ipdsum_dict[kmer])))
        print ("kmer", kmer, mae, mae_ipd_sum, pearsonr(real_dict[kmer], estimate_dict[kmer]), pearsonr(real_dict[kmer], ipdsum_dict[kmer]))
        print ("real", real_dict[kmer][:10], np.mean(real_dict[kmer]))
        print ("estimated", estimate_dict[kmer][:10], np.mean(estimate_dict[kmer]))
        print ("ipdsum", ipdsum_dict[kmer][:10], np.mean(ipdsum_dict[kmer]))
        print ("")
        break


# log = "slurm-705929.out"
# read_log(log)
csv = "/home/shuaiw/methylation/data/borg/human/test_result6.csv"
read_log_csv(csv)
# count_log_csv(csv)