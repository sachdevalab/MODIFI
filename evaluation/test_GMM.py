import numpy as np
import warnings
from sklearn.mixture import GaussianMixture
from sklearn.exceptions import ConvergenceWarning
import pandas as pd
from scipy.stats import norm


# Fit a gaussian mixture model to the ipd of a single molecule, return the binarized values
def fitGaussian(ipds, initmean = False, zmw=-1):
    with np.errstate(invalid='ignore'):
        ridx = np.where(ipds > 0)[0]

    # lrnn = np.log(ipds[ridx]).reshape(-1,1)
    lrnn = ipds[ridx].reshape(-1,1)
    print ("len(lrnn)", len(lrnn))
    gmm2 = GaussianMixture(2, covariance_type='spherical')
    gmm2.fit(lrnn)
    print ("xxx", gmm2.means_, gmm2.weights_, gmm2.converged_)

    gmm1 = GaussianMixture(1, covariance_type='spherical')
    if initmean:
        gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                              weights_init=np.array([.85, .15]), tol=1e-6)
    else:
        gmm2 = GaussianMixture(2, covariance_type='spherical')
    
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        gmm1.fit(lrnn)
        if not gmm1.converged_:
            print('zmw #%d did not converge on gmm1' % (zmw))
        gmm2.fit(lrnn)
        print (gmm2.means_, gmm2.weights_, gmm2.converged_)
        convround = 0
        if not gmm2.converged_:
            gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                  n_init=10, init_params='random', tol=1e-5, max_iter=200)
            gmm2.fit(lrnn)
            convround = 1
            if not gmm2.converged_:
                gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                      n_init=20, init_params='random', tol=1e-4, max_iter=400)
                gmm2.fit(lrnn)
                convround = 2
                if not gmm2.converged_:
                    gmm2 = GaussianMixture(2, covariance_type='spherical', means_init=np.array([1.1, 3]).reshape((2,1)),
                                      n_init=40, init_params='random', tol=1e-3, max_iter=600)
                    gmm2.fit(lrnn)
                    convround = 3
                    if not gmm2.converged_:
                        convround = 9
                        print('zmw #%d did not converge on gmm2 even after extensions' % (zmw))
    print (gmm1.aic(lrnn), gmm2.aic(lrnn))
    aicdif = gmm1.aic(lrnn) - gmm2.aic(lrnn)
    mixmeans = gmm2.means_.flatten()
    mixweights = gmm2.weights_.flatten()
    elevstate = np.argmax(mixmeans)
    resp = gmm2.predict_proba(np.log(ipds[ridx].reshape(-1,1)))
    respfull = np.empty(len(ipds), dtype='float32')
    respfull.fill(np.nan)
    respfull[ridx] = resp[:,elevstate]
    convergInf = str(convround) + '.' + str(gmm2.n_iter_)
    print ("////")
    print (mixmeans, mixweights, aicdif, convergInf)
    ## print whether the model is 1 GMM or 2 GMM
    if aicdif >= 0:
        print ('2 GMM')
    else:
        print('1 GMM')
    
    return (respfull, np.array([mixmeans[1-elevstate], mixmeans[elevstate]]), 
            np.array([mixweights[1-elevstate], mixweights[elevstate]]), aicdif, convergInf)
def read_gff(gff):
    modified_loci = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            locus_tag = cols[0] + "_" + str(cols[3]) + "_" + cols[6]
            score = int(cols[5])
            # if score < 30:
            #     continue
            modified_loci[locus_tag] = score
    return modified_loci
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

# ## test the function
# file = "/home/shuaiw/borg/bench/C227/WGA2/ipd_ratio/CP011331.1.ipd3.csv"
file = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/ipd_ratio/E_coli_H10407_2.ipd3.csv" 
# file = "/home/shuaiw/borg/bench/breakdb/ipd_ratio/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_494_C_0_547023.ipd3.csv"
## read the ipd file using pd

df = pd.read_csv(file)
df = df[df['ipd_ratio'] != 0]

# X2 = df['ipd_ratio'].values.reshape(-1,1)
# gmm2 = GaussianMixture(2, weights_init=np.array([.99, .01]), means_init=np.array([1.1, 3]).reshape((2,1)))
# gmm2.fit(X2)
# print (gmm2.means_, gmm2.aic(X2), gmm2.weights_.flatten())
# ## add a pvalue column to the dataframe, indicating the p value of a ipd_ratio belong to the lower distribution
# df['pvalue'] = gmm2.predict_proba(X2)[:,0]
# print (df)


mean = df['ipd_ratio'].mean()
std = df['ipd_ratio'].std()
print (f"mean: {mean}, std: {std}, cal pvalue...")

df['pvalue'] = df['ipd_ratio'].apply(lambda x: p_value_right_tail(x, mean, std))

## count how much with pvalue < 0.05
print (len(df[df['pvalue'] < 0.05]))

my_loci = {}
for index, row in df.iterrows():
    if row['pvalue'] < 0.05:
        if row['strand'] == '0':
            row['strand'] = '-'
        else:
            row['strand'] = '+'
        one_based_tpl = row['tpl'] + 1
        tag = row['refName'] + "_" + str(one_based_tpl) + "_" + row['strand']
        # print (index, row['ipd_ratio'], row['pvalue'])
        my_loci[tag] = row['pvalue']

# print (my_loci)



bench_gff = "/home/shuaiw/borg/bench/zymo_new_ref_NM3/gffs/E_coli_H10407_2.gff"
# depth_gff = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/gffs/E_coli_H10407_1.gff" 



bench_loci = read_gff(bench_gff)
# print (bench_loci)
depth_loci = my_loci
# depth_loci = read_gff(depth_gff)
## use bench loci as reference, calculate recall and precision
recall_num  = 0
precision_num = 0
for locus in bench_loci:
    if locus in depth_loci:
        recall_num += 1
    # else:
    #     print ("Locus not in depth loci: ", locus)
## calculate recall
recall = recall_num / len(bench_loci)
## calculate precision
for locus in depth_loci:
    if locus in bench_loci:
        precision_num += 1
precision = precision_num / len(depth_loci)
## calculate F1 score
f1 = 2 * (precision * recall) / (precision + recall)
print ("Recall: ", recall)
print ("Precision: ", precision)
print ("F1 score: ", f1)









# ## generate random data with one or two gaussians
# np.random.seed(0)
# n_samples = 10000
# # generate random data with one gaussian
# # Generate samples from a normal distribution
# X1 = np.random.randn(n_samples) * 0.5 + 1
# X1 = X1.reshape(-1,1)

# X2 = np.random.randn(100) * 1 + 5
# X2 = X2.reshape(-1,1)

# # # generate random data with two gaussians
# X2 = np.vstack([X1, X2])
# print (X2)
# # ## fit the model
# respfull, mixmeans, mixweights, aicdif, convergInf = fitGaussian(X2.flatten(), True)
# # print (respfull, mixmeans, mixweights, aicdif, convergInf)
# print (aicdif)
# gmm2 = GaussianMixture(2, weights_init=np.array([.99, .01]), means_init=np.array([1.1, 3]).reshape((2,1)))
# gmm2.fit(X2)
# print (gmm2.means_, gmm2.aic(X2), gmm2.weights_.flatten())

# gmm1 = GaussianMixture(1)
# gmm1.fit(X2)
# print (gmm1.means_, gmm1.aic(X2))
# ## plot the distribution of X2
# import matplotlib.pyplot as plt
# import seaborn as sns
# sns.histplot(X2)
# ## save the plot
# plt.savefig('../tmp/test.png')