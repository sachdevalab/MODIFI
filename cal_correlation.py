import pandas as pd
from scipy.stats import pearsonr

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
    df = pd.read_csv(csv)
    print ("length of df", len(df))
    ## cal correlation between estimated column and real column
    corr = pearsonr(df['estimated'], df['real'])
    print ("correlation", corr)
    corr = pearsonr(df['median'], df['real'])
    print ("correlation", corr)
    ## cal correlation between estimated column and ipdsum column
    corr = pearsonr(df['ipd_sum'], df['real'])
    print ("correlation", corr)

# log = "slurm-705929.out"
# read_log(log)
csv = "/home/shuaiw/methylation/data/borg/human/test_result.csv"
read_log_csv(csv)