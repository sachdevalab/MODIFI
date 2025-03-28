import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import matplotlib.pyplot as plt



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

def truth_comp(our_tmeans, our_ipd_ratio, our_controls, true_controls):
    ### get ipd ratio for true controls
    true_ipd_ratio_dict = {}
    our_ipd_ratio_dict = {}
    for tpl in true_controls:
        if tpl in our_tmeans and tpl in our_ipd_ratio:
            if true_controls[tpl] == 0:
                continue
            true_ipd_ratio_dict[tpl] = our_tmeans[tpl] / true_controls[tpl]
            our_ipd_ratio_dict[tpl] = our_ipd_ratio[tpl]
    
    truth_ipd_ratio = list(true_ipd_ratio_dict.values())
    true_mean_ipd = np.mean(truth_ipd_ratio)
    true_std_ipd = np.std(truth_ipd_ratio)
    true_cutoff = calculate_x_from_pvalue(0.05, true_mean_ipd, true_std_ipd)
    true_modified = [x for x in true_ipd_ratio_dict if true_ipd_ratio_dict[x] > true_cutoff]

    our_ipd_ratio = list(our_ipd_ratio_dict.values())
    our_mean_ipd = np.mean(our_ipd_ratio)
    our_std_ipd = np.std(our_ipd_ratio)
    our_cutoff = calculate_x_from_pvalue(0.05, our_mean_ipd, our_std_ipd)
    our_modified = [x for x in our_ipd_ratio_dict if our_ipd_ratio_dict[x] > our_cutoff]

    ## calulate the recall and precision for our method
    recall = 0
    for tpl in true_modified:
        if tpl in our_modified:
            recall += 1
        else:
            print (tpl, true_ipd_ratio_dict[tpl], true_cutoff, our_ipd_ratio_dict[tpl], our_cutoff)
    print ("recall", recall/len(true_modified), len(true_modified), len(our_modified))
    print ("precision", recall/len(our_modified))

def corr_obs_cont(infer):
    df = pd.read_csv(infer)
    # df = df[df['kmer_count'] > 10]
    df = df[df['strand'] == 1]
    ## plot the distribution of ipd_ratio
    df['ipd_ratio'].plot.hist(bins=100)
    ## svae the plot
    plt.savefig("tmp/ipd_ratio.png")
    ## clean the plt
    plt.clf()
    ## plot the tMean
    df['tMean'].plot.hist(bins=100)
    ## svae the plot
    plt.savefig("tmp/tMean.png")
    plt.clf()
    ## cal the correlation between tMean and control
    corr = pearsonr(df['tMean'], df['control'])
    print ("tMean and control", corr)

def corr_ipd_sum(ipd_summary):
    df = pd.read_csv(ipd_summary, nrows=1705173*2)

    # df1 = df[df['strand'] == 1]
    # df2 = df[df['strand'] == 0]
    # ## align df1 and df2 with same tpl, remove the elements with tpl not in df2
    # df1 = df1[df1['tpl'].isin(df2['tpl'])]
    # df2 = df2[df2['tpl'].isin(df1['tpl'])]

    # print (pearsonr(df1['tMean'], df2['modelPrediction']))
    # print (pearsonr(df1['modelPrediction'], df2['tMean']))
    # print (pearsonr(df1['modelPrediction'], df1['tMean']))



    df = df[df['strand'] == 1]
    ## cal the correlation between tMean and control
    corr = pearsonr(df['tMean'], df['modelPrediction'])


    print ("tMean and control", corr)

def corr_depth(normal, p5):
    df = pd.read_csv(normal, nrows=1000)
    df2 = pd.read_csv(p5, nrows=1000)
    print (pearsonr(df['tMean'], df2['tMean']))
    print  (pearsonr(df['control'], df2['control']))
    print  (pearsonr(df['ipd_ratio'], df2['ipd_ratio']))

    ## print first 10 tMean 
    for i in range(100):
        print (df['tMean'][i],  df2['tMean'][i], df['ipd_ratio'][i], df2['ipd_ratio'][i])

if __name__ == "__main__":

    # infer = "/home/shuaiw/methylation/data/borg/b_contigs/ipds4/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C_0_852595.ipd2.csv"
    # ipd_summary = "/home/shuaiw/methylation/data/borg/b_contigs/1.csv"

    # contig = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_882_L"
    # ipd_summary = "/home/shuaiw/methylation/data/borg/test_100/test_100.csv"
    # infer = f"/home/shuaiw/methylation/data/borg/new_test7/ipd_ratio/{contig}.ipd3.csv"

    # contig = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C"
    # ipd_summary = "/home/shuaiw/borg/seven_contigs/seven_contigs.csv"
    # infer = f"/home/shuaiw/borg/seven_test/ipd_ratio/{contig}.ipd3.csv"

    # bench = Benchmark(ipd_summary, infer)
    # bench.read_ipd_summary_all()
    # bench.read_our_all()
    # bench.compare_ctgs()
    # infer = "/home/shuaiw/borg/new_test10/ipd_ratio/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd3.csv"
    # ipd_summary = "/home/shuaiw/methylation/data/borg/b_contigs/11.csv"
    # infer = "/home/shuaiw/borg/bench/ecoli_native2/ipd_ratio/CP064387.1.ipd3.csv"
    # corr_obs_cont(infer)
    # corr_ipd_sum(ipd_summary)

    # infer = "/home/shuaiw/borg/bench/break/test/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C_0_852595.ipd2.csv"
    # corr_obs_cont(infer)
    # ipd_summary = "/home/shuaiw/borg/break_contigs/break_contigs.csv"
    # corr_ipd_sum(ipd_summary)

    p5 = "/home/shuaiw/borg/bench/zymo_new_ref_p0.1/ipd_ratio/E_coli_H10407_1.ipd3.csv"
    normal = "/home/shuaiw/borg/bench/zymo_new_ref/ipd_ratio/E_coli_H10407_1.ipd3.csv"
    corr_depth(normal, p5)


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

    # strand = 1
    # # df = pd.read_csv("/home/shuaiw/methylation/data/borg/bench/ecoli_native/ipd_ratio/CP064387.1.ipd3.csv", nrows=200000)
    # df = pd.read_csv("/home/shuaiw/methylation/data/borg/bench/merge/ipd_ratio/CP064387.1.ipd3.csv")
    # ## only keep the rows with pvalue < 0.01
    # # df = df[df['pvalue'] < 0.01]
    
    # print (len(df))
    # our_controls = {}
    # our_ipd_ratio = {}
    # our_tmeans = {}
    # true_controls = {}
    # df = df[df['strand'] == strand]
    # df = df[df['coverage'] > 5]
    # for index, row in df.iterrows():
    #     our_controls[row['tpl']] = row['control']
    #     our_ipd_ratio[row['tpl']] = row['ipd_ratio']
    #     our_tmeans[row['tpl']] = row['tMean']

    # df = pd.read_csv("/home/shuaiw/methylation/data/borg/bench/ecoli_control/ipd_ratio/CP064387.1.ipd3.csv")
    # ## only keep the rows with pvalue < 0.01
    # # df = df[df['pvalue'] < 0.01]
    # print (len(df))
    # df = df[df['strand'] == strand]
    # df = df[df['coverage'] > 5]
    # for index, row in df.iterrows():
    #     true_controls[row['tpl']] = row['tMean']
    # print (len(our_controls), len(true_controls))

    
    # truth_comp(our_tmeans, our_ipd_ratio, our_controls, true_controls)

    # ## ipd summary result
    # ipd_result = "/home/shuaiw/borg/bench/ecoli_control/contigs/CP064387.1.csv"
    # df = pd.read_csv(ipd_result, sep = ";", nrows=4000000)
    # df = df[df['Strand'] == strand]
    # ipd_controls = {}
    # for index, row in df.iterrows():
    #     ipd_controls[row['Position']-1] = row['Prediction']


    # ipd_control_list = []
    # our_control_list = []
    # true_control_list = []
    # for tpl in our_controls:
    #     if tpl in true_controls:
    #         our_control_list.append(our_controls[tpl])
    #         true_control_list.append(true_controls[tpl])
    #         ipd_control_list.append(ipd_controls[tpl])
    #         # print (tpl, our_controls[tpl], true_controls[tpl])
    # print (len(our_controls), len(true_controls), len(ipd_controls), len(our_control_list), len(true_control_list), len(ipd_control_list))
    # print (len(our_control_list), pearsonr(our_control_list, true_control_list))
    
    # print (len(ipd_control_list), pearsonr(ipd_control_list, true_control_list))
    # # print (len(ipd_control_list), pearsonr(ipd_control_list, our_control_list))

