import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

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

def plot(predicted_list, real_list, file_name):
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr

    # Create a scatter plot using seaborn
    sns.scatterplot(x=predicted_list, y=real_list, s=10)

    # Add x label and y label
    plt.xlabel("Predicted IPD")
    plt.ylabel("Observed IPD")

    # Calculate the Pearson correlation coefficient
    corr, _ = pearsonr(predicted_list, real_list)
    print("Pearson correlation:", corr)

    # Add the correlation as text on the plot
    plt.text(0.05, 0.95, f'Pearson correlation: {corr:.2f}', transform=plt.gca().transAxes, 
             fontsize=8, verticalalignment='top')

    # # Set x and y limits
    # plt.xlim(0, 5)
    # plt.ylim(0, 5)

    # Make sure the x and y axis are the same
    plt.axis('equal')

    # Save the plot
    plt.savefig(file_name)

    # Close the plot
    plt.close()

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
    return df


def p_value_right_tail(x, mu, sigma):
    """
    Calculate the one-tailed p-value for x being in the right tail of a normal distribution.

    Parameters:
    x (float): The observed value
    mu (float): The mean of the normal distribution
    sigma (float): The standard deviation of the normal distribution

    Returns:
    float: The p-value for x being in the right tail
    """
    # Calculate the z-score
    z = (x - mu) / sigma
    
    # Calculate the right-tail p-value
    p_value = 1 - norm.cdf(z)
    
    return p_value[0]

def calculate_x_from_pvalue(p_value, mu, sigma, tail="right"):
    """
    Calculate the value(s) x that correspond to a given p-value in a normal distribution.

    Parameters:
    p_value (float): The p-value (probability)
    mu (float): The mean of the normal distribution
    sigma (float): The standard deviation of the normal distribution
    tail (str): The tail type. Options are "right", "left", or "two-sided".

    Returns:
    float or tuple: The value(s) x corresponding to the given p-value
    """
    if tail == "right":
        # Right-tail: Find the value where P(Z >= z) = p_value
        z = norm.ppf(1 - p_value)
        x = mu + z * sigma
        return x

def get_ipd_ratio(df):
    ## get ratio which is the ratio between estimated and ipdsum
    df['ipd_ratio'] = df['real'] / df['estimated'] 
    ## plot the distribution of ipd_ratio and save the plot in pdf
    ## use log scale for y axis
    import matplotlib.pyplot as plt
    ## add x label and y label
    plt.hist(df['ipd_ratio'], bins=100)
    plt.xlabel("ipd_ratio")
    plt.ylabel("count")
    plt.savefig("tmp/ipd_ratio.pdf")
    # Fit GMM model
    data = df['ipd_ratio'].values.reshape(-1, 1)
    gmm = GaussianMixture(n_components=1, covariance_type='full', random_state=42)
    gmm.fit(data)
    ## print std and mean
    # print(gmm.means_, np.sqrt(gmm.covariances_))
    # print (data[:10])
    mean = gmm.means_[0][0]
    std = np.sqrt(gmm.covariances_[0][0])
    ## ofr each value, calculate the probability belong to the model
    cutoff = calculate_x_from_pvalue(0.01, mean, std)
    print ("cutoff", cutoff)
    ## add pvalue column to the dataframe
    df['pvalue'] = df['ipd_ratio'].apply(lambda x: p_value_right_tail(x, mean, std))
    print (df)
    ## filter df with pvalue < 0.01
    df = df[df['pvalue'] < 0.01]
    ## svae the df to tmp/ipd_ratio_methy.csv
    df.to_csv("tmp/ipd_ratio_methy.csv", index=False)




def get_raw_ipd(df):
    ## plot the distribution of ipd_ratio and save the plot in pdf
    ## use log scale for y axis
    import matplotlib.pyplot as plt
    ## add x label and y label
    plt.hist(df['real'], bins=100)
    plt.xlabel("Observed_IPD")
    plt.ylabel("count")
    plt.savefig("tmp/raw_IPD.pdf")

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
    plot(df['estimated'], df['real'], "tmp/our_scatter.pdf")
    plot(df['ipd_sum'], df['real'], "tmp/ipd_sum_scatter.pdf")

    # estimate_dict = defaultdict(list)
    # ipdsum_dict = defaultdict(list)
    # real_dict = defaultdict(list)
    # ## enumerate the row
    # for index, row in df.iterrows():
    #     # print(index, row['kmer'])
    #     estimate_dict[row['kmer']].append(row['estimated'])
    #     ipdsum_dict[row['kmer']].append(row['ipd_sum'])
    #     real_dict[row['kmer']].append(row['real'])
    # ## sort the kmer with the length of value list
    # sorted_estimate = sorted(estimate_dict.items(), key=lambda x: len(x[1]), reverse=True)
    # for i in range(10):
    #     kmer = sorted_estimate[i][0]
    #     ## round two for each element in estimate_dict[kmer]
    #     estimate_dict[kmer] = [round(x, 2) for x in estimate_dict[kmer]]
    #     ## cal MAE between real and estimated
    #     mae = np.mean(np.abs(np.array(real_dict[kmer]) - np.array(estimate_dict[kmer])))
    #     mae_ipd_sum = np.mean(np.abs(np.array(real_dict[kmer]) - np.array(ipdsum_dict[kmer])))
    #     print ("kmer", kmer, mae, mae_ipd_sum, pearsonr(real_dict[kmer], estimate_dict[kmer]), pearsonr(real_dict[kmer], ipdsum_dict[kmer]))
    #     print ("real", real_dict[kmer][:10], np.mean(real_dict[kmer]))
    #     print ("estimated", estimate_dict[kmer][:10], np.mean(estimate_dict[kmer]))
    #     print ("ipdsum", ipdsum_dict[kmer][:10], np.mean(ipdsum_dict[kmer]))
    #     print ("")
    #     break


# log = "slurm-705929.out"
# read_log(log)
csv = "/home/shuaiw/methylation/data/borg/human/test_result5.csv"
df = read_log_csv(csv)
count_log_csv(csv)
# get_ipd_ratio(df)
# get_raw_ipd(df)