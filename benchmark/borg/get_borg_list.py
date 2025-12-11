from Bio import SeqIO
import pandas as pd

# fasta = "/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.soil_1.fa"
# borg_list = "/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.soil_1.borg_list.txt"

fasta = "/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META.soil_2.fa"
borg_list = "/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META.soil_2.borg_list.txt"

## read the fasta using biopython, collect the sequence names with annotations contain any "Borg",
## output it to a list file

print(f"Reading FASTA file: {fasta}")
print(f"Output file: {borg_list}")

data = []

borg_count = 0
total_count = 0

for record in SeqIO.parse(fasta, "fasta"):
    total_count += 1
    if "Borg" in record.description or "Jumbo" in record.description or "Phage" in record.description:
        length = len(record.seq)
        data.append([record.id, record.description, length])

        borg_count += 1
        print(f"Found BORG: {record.id} - {record.description} (Length: {length})")
print(f"Total sequences: {total_count}")
print(f"BORG sequences found: {borg_count}")
print(f"BORG list saved to: {borg_list}")

## save the data to a csv file

df = pd.DataFrame(data, columns=["seq_name", "description", "length"])
df.to_csv(borg_list, index=False, sep="\t")