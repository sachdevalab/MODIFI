from collections import defaultdict
import pandas as pd

sample_occur_times = defaultdict(int)

f = open("rebase_bam_data.csv", )
f.readline()
for line in f:
    field = line.strip().split(',')
    sample_name = field[2]
    sample_occur_times[sample_name] += 1

## sort the dictionary by value
sorted_sample_occur_times = sorted(sample_occur_times.items(), key=lambda x: x[1], reverse=True)
for i in range(50):
    print (sorted_sample_occur_times[i])