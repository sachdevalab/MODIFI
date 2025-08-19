import pandas as pd
from Bio import SeqIO


linkage_file = "/home/shuaiw/borg/paper/run2/ocean_1/ocean_1_methylation2/host_summary.csv"
assembly = "/home/shuaiw/borg/paper/run2/ocean_1/ocean_1.hifiasm.p_ctg.rename.fa"
focus_ctg = "ocean_1_2090_L"
focus_ctgs_fa = "/home/shuaiw/borg/paper/run2/ocean_1/focus_ctgs.fa"


df = pd.read_csv(linkage_file)
focus_ctgs = set()

for index, row in df.iterrows():
    if row['MGE'] == focus_ctg or row['host'] == focus_ctg:
        focus_ctgs.add(row['host'])
        focus_ctgs.add(row['MGE'])

## extract focus_ctgs from assembly and output to focus_ctgs_fa
## using biopython

with open(focus_ctgs_fa, 'w') as out_f:
    for record in SeqIO.parse(assembly, "fasta"):
        if record.id in focus_ctgs:
            SeqIO.write(record, out_f, "fasta")
