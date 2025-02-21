"""
ipdtools predict the control IPD values for a fasta
Now we want to merge it with the loaded raw IPD values
"""


import pandas as pd
import numpy as np
import sys


#  ipdtools -f /home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa -o /home/shuaiw/methylation/data/borg/new_test10/test.csv --nproc 10 --indexing 0

# raw_ipd = "/home/shuaiw/methylation/data/borg/new_test10/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd1.csv"
# predicted_control = "/home/shuaiw/methylation/data/borg/new_test10/test.csv"
# output = "/home/shuaiw/methylation/data/borg/new_test10/test_merge.csv"

raw_ipd = sys.argv[1]
predicted_control = sys.argv[2]
output = sys.argv[3]

df_raw = pd.read_csv(raw_ipd)
df_control = pd.read_csv(predicted_control, sep = ";")

## if df_control["strand"] matched df_raw["strand"], and df_control["Position"] matched df_raw["tpl"], then add df_control["Prediction"] to df_raw["Prediction"]
# Perform a merge operation to match 'Strand' and 'Position'
df_merged = df_raw.merge(df_control[["Strand", "Position", "Prediction"]],
                         left_on=["strand", "tpl"], 
                         right_on=["Strand", "Position"], 
                         how="left")

# Replace the control column with merged Prediction values
df_raw["control"] = df_merged["Prediction"]
# Save the result
df_raw.to_csv(output, index=False)