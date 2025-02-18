import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import sys



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

def get_ipd_ratio(csv, output):
    df = pd.read_csv(csv, sep = ",")
    mean = df['ipd_ratio'].mean()
    std = df['ipd_ratio'].std()
    df['pvalue'] = df['ipd_ratio'].apply(lambda x: p_value_right_tail(x, mean, std))

    ### round all the float values to 2 decimal places
    df['tpl'] = df['tpl'].astype(int)
    df['strand'] = df['strand'].astype(int)
    df['coverage'] = df['coverage'].astype(int)
    df['kmer_count'] = df['kmer_count'].astype(int)
    df = df.round(4)

    df.to_csv(output, index=False)


if __name__ == "__main__":
    csv = sys.argv[1]
    output = sys.argv[2]
    get_ipd_ratio(csv, output)


    # csv = "/home/shuaiw/methylation/data/borg/b_contigs/control/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd2.csv"
    # output = "/home/shuaiw/methylation/data/borg/b_contigs/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd3.csv"
    # get_ipd_ratio(csv, output)
    # get_raw_ipd(df)

