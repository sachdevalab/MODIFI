hic_file = "/home/shuaiw/borg/pengfan/hic_10mgs_linkages.csv"
# host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/host_summary.csv"
host_sum = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_Comb_RF_LR_bin/host_summary.csv"

import pandas as pd

def read_hic(hic_file):
    hic_linkages = {}
    df = pd.read_csv(hic_file)
    for index, row in df.iterrows():
        hic_linkages[row['circular_element']] = row['bin']
    return hic_linkages

def read_our(host_sum, score_cutoff = 0.45):
    df = pd.read_csv(host_sum)
    df = df[df['final_score'] > score_cutoff]
    our_linkages = {}
    for index, row in df.iterrows():
        our_linkages[row['plasmid']] = row['host']
    return our_linkages

def compare_hic_our(hic_linkages, our_linkages):
    both_link = 0
    cosistency_num = 0
    for plasmid in our_linkages:
        if plasmid not in hic_linkages:
            # print(f"{plasmid} is not in HiC linkages")
            continue
        both_link += 1
        if our_linkages[plasmid] == hic_linkages[plasmid]:
            cosistency_num += 1
        else:
            print (f"{plasmid} is not consistent: {our_linkages[plasmid]} vs {hic_linkages[plasmid]}")
    if both_link == 0:
        cosistency_rate = 0
    else:
        cosistency_rate = cosistency_num / both_link
    print(f"both linkages: {both_link}")
    print(f"cosistency linkages: {cosistency_num}")
    print (f"cosistency rate: {cosistency_rate}")
    print ("our link num", len(our_linkages))
    return both_link, cosistency_num, cosistency_rate, len(our_linkages)

data = []
for cutoff in range(0, 16):
    my_cutoff = cutoff / 20
    hic_linkages = read_hic(hic_file)
    our_linkages = read_our(host_sum, my_cutoff)
    print (f"cutoff: {my_cutoff}")
    both_link, cosistency_num, consis_rate, our_num = compare_hic_our(hic_linkages, our_linkages)
    data.append([my_cutoff, both_link, cosistency_num, consis_rate, our_num])
df = pd.DataFrame(data, columns=['cutoff', 'both_link', 'cosistency_num', 'consis_rate', 'our_num'])

### plot the consistency rate line plot
## plot another subplot with the our link num
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid")


plt.figure(figsize=(5, 9))  # Adjust the figure size to accommodate the third subplot

## plot three subplots
plt.subplot(3, 1, 1)
plt.plot(df['cutoff'], df['consis_rate'], marker='o')
plt.title('Ratio of Consistency Linkages')
plt.xlabel('Cutoff')
plt.ylabel('Ratio of Consistency Linkages')

plt.subplot(3, 1, 2)
plt.plot(df['cutoff'], df['cosistency_num'], marker='o')
plt.title('Number of Consistency Linkages')
plt.xlabel('Cutoff')
plt.ylabel('Number of Consistency Linkages')

plt.subplot(3, 1, 3)
plt.plot(df['cutoff'], df['our_num'], marker='o', color='green')  # Plot our_num
plt.title('Number of Our Linkages')
plt.xlabel('Cutoff')
plt.ylabel('Number of Our Linkages')

plt.tight_layout()
plt.savefig('../tmp/both_linkages.png')

# plt.figure(figsize=(5, 6))
# ## plot two subplots
# plt.subplot(2, 1, 1)
# plt.plot(df['cutoff'], df['consis_rate'], marker='o')
# plt.title('Ratio of Consistency Linkages')
# plt.xlabel('Cutoff')
# plt.ylabel('Ratio of Consistency Linkages')
# plt.subplot(2, 1, 2)
# plt.plot(df['cutoff'], df['cosistency_num'], marker='o')
# plt.title('Number of Consistency Linkages')
# plt.xlabel('Cutoff')
# plt.ylabel('Number of Consistency Linkages')

# plt.tight_layout()
# plt.savefig('../tmp/both_linkages.png')

# plt.plot(df['cutoff'], df['consis_rate'], marker='o')
# plt.title('Consistency Rate of HiC and Our Linkages')
# plt.xlabel('Cutoff')
# plt.ylabel('Consistency Rate')
# plt.savefig('../tmp/consistency_rate.png')