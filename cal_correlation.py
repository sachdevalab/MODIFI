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
    
    return p_value

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

class Benchmark:
    def __init__(self, ipd_summary, our):
        self.ipd_summary = ipd_summary
        self.our = our
        self.ipd_dict = {}
        self.ipd_dict_ctgs = defaultdict(dict)
        self.our_dict_ctgs = defaultdict(dict)
        self.our_dict = {}
        self.save_our_all = {}
        self.ratio_ana_dict = {}

    def read_ipd_summary(self):
        df = pd.read_csv(self.ipd_summary)
        ## keep the rows with strand == 1
        df = df[df['strand'] == 1]
        self.ratio_ana(df)
        ## keep the rows with score > 30
        df = df[df['score'] > 30]
        for index,row in df.iterrows():
            self.ipd_dict[row['tpl']] = row['score']

    def read_ipd_summary_all(self):
        df = pd.read_csv(self.ipd_summary)
        print ("read is done")
        ## keep the rows with strand == 1
        df = df[df['strand'] == 1]
        df = df[df['refName'] == contig]
        self.ratio_ana(df)
        ## keep the rows with score > 30
        df = df[df['score'] > 30]
        for index,row in df.iterrows():
            self.ipd_dict_ctgs[row['refName']][row['tpl']] = row['score']
            
    def ratio_ana(self, df):
        mean = df['ipdRatio'].mean()
        std = df['ipdRatio'].std()
        ## ofr each value, calculate the probability belong to the model
        ## add pvalue column to the dataframe
        df['pvalue'] = df['ipdRatio'].apply(lambda x: p_value_right_tail(x, mean, std))
        for index, row in df.iterrows():
            if row['pvalue'] < 0.05:
                self.ratio_ana_dict[row['tpl']] = row['pvalue']
        
    def read_our(self):
        df = pd.read_csv(self.our)
        ## only keep the rows with pvalue < 0.01
        for index, row in df.iterrows():
            if row['pvalue'] < 0.05:
                self.our_dict[row['tpl']+1] = row['ipd_ratio']

            self.save_our_all[row['tpl']+1] = row

    def read_our_all(self):
        df = pd.read_csv(self.our)
        df = df[df['strand'] == 1]
        ## only keep the rows with pvalue < 0.01
        for index, row in df.iterrows():
            if row['pvalue'] < 0.05:
                # self.our_dict[row['tpl']+1] = row['ipd_ratio']
                self.our_dict_ctgs[row['refName']][row['tpl']+1] = row['ipd_ratio']

            # self.save_our_all[row['tpl']+1] = row
    
    def compare_ctgs(self):
        
        for ctg in self.ipd_dict_ctgs:
            if ctg not in self.our_dict_ctgs:
                continue
            recall = 0
            for tpl in self.ipd_dict_ctgs[ctg]:
                if tpl in self.our_dict_ctgs[ctg]:
                    recall += 1
                # else:
                #     print (tpl, self.ipd_dict[tpl], "not in our")
                    # if tpl in self.save_our_all:
                    #     print (self.save_our_all[tpl])
            print (ctg)
            print (recall, len(self.ipd_dict_ctgs[ctg]), len(self.our_dict_ctgs[ctg]))
            print ("recall", recall / len(self.ipd_dict_ctgs[ctg]))
            print ("precision", recall / len(self.our_dict_ctgs[ctg]))

    
    def compare(self):
        recall = 0
        for tpl in self.ipd_dict:
            if tpl in self.our_dict:
                recall += 1
            # else:
            #     print (tpl, self.ipd_dict[tpl], "not in our")
                # if tpl in self.save_our_all:
                #     print (self.save_our_all[tpl])
        print (recall, len(self.ipd_dict), len(self.our_dict))
        print ("recall", recall / len(self.ipd_dict))
        print ("precision", recall / len(self.our_dict))

    def compare2(self):
        recall = 0
        for tpl in self.ratio_ana_dict:
            if tpl in self.our_dict:
                recall += 1
            # else:
            #     print (tpl, self.ratio_ana_dict[tpl], "not in our")
                # if tpl in self.save_our_all:
                #     print (self.save_our_all[tpl])
        print (recall, len(self.ratio_ana_dict), len(self.our_dict))
        print ("recall", recall / len(self.ratio_ana_dict))
        print ("precision", recall / len(self.our_dict))

    def compare3(self):
        recall = 0
        for tpl in self.ratio_ana_dict:
            if tpl in self.ratio_ana_dict:
                recall += 1
            else:
                print (tpl, self.ipd_dict[tpl], "not in our")
                # if tpl in self.save_our_all:
                #     print (self.save_our_all[tpl])
        print (recall, len(self.ipd_dict), len(self.ratio_ana_dict))
        print ("recall", recall / len(self.ipd_dict))
        print ("precision", recall / len(self.ratio_ana_dict))


if __name__ == "__main__":

    # infer = "/home/shuaiw/methylation/data/borg/b_contigs/ipds4/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C_0_852595.ipd2.csv"
    # ipd_summary = "/home/shuaiw/methylation/data/borg/b_contigs/1.csv"

    # contig = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_882_L"
    # ipd_summary = "/home/shuaiw/methylation/data/borg/test_100/test_100.csv"
    # infer = f"/home/shuaiw/methylation/data/borg/new_test7/ipd_ratio/{contig}.ipd3.csv"

    # bench = Benchmark(ipd_summary, infer)
    # bench.read_ipd_summary_all()
    # bench.read_our_all()
    # bench.compare_ctgs()


    # bench.read_ipd_summary()
    # bench.read_our()
    # bench.compare()
    # bench.compare2()
    # bench.compare3()

    # log = "slurm-705929.out"
    # read_log(log)
    # csv = "/home/shuaiw/methylation/data/borg/human/test_result5.csv"
    # df = read_log_csv(csv)
    # count_log_csv(csv)
    # get_ipd_ratio(df)
    # get_raw_ipd(df)

    df = pd.read_csv("/home/shuaiw/methylation/data/borg/bench/ecoli_control/ipd_ratio/CP064387.1.ipd3.csv")
    ## only keep the rows with pvalue < 0.01
    df = df[df['pvalue'] < 0.01]
    print (len(df))
    # for index, row in df.iterrows():
    #     if row['pvalue'] < 0.05:
    #         print (row)
    df = pd.read_csv("/home/shuaiw/methylation/data/borg/bench/ecoli_native/ipd_ratio/CP064387.1.ipd3.csv")
    ## only keep the rows with pvalue < 0.01
    df = df[df['pvalue'] < 0.01]
    print (len(df))