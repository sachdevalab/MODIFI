import pandas as pd
import re

motif_profile = "/home/shuaiw/borg/paper/linkage/pure2/m64004_210929_143746.p100/motif_profile.csv"
plasmid_list_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list"

# Read the motif profile
df = pd.read_csv(motif_profile)

# Read the plasmid list file
plasmid_df = pd.read_csv(plasmid_list_file)

# Create set of all plasmids (seq_name column)
plasmid_set = set(plasmid_df['seq_name'].tolist())

# Create set of all hosts - need to split the host column by semicolon and collect all
host_set = set()
for hosts in plasmid_df['host'].tolist():
    if pd.notna(hosts):
        host_set.update([h.strip() for h in str(hosts).split(';')])

print(f"Found {len(plasmid_set)} plasmids")
print(f"Found {len(host_set)} hosts")

# Transform from wide to long format
df_long = df.melt(id_vars=['motif_identifier'], 
                  var_name='sample', 
                  value_name='fraction')

# Extract strain name by removing the trailing _number pattern
df_long['strain_name'] = df_long['sample'].str.replace(r'_\d+$', '', regex=True)

# Add genome_type column
def assign_genome_type(sample_name):
    if sample_name in plasmid_set:
        return 'plasmid'
    elif sample_name in host_set:
        return 'host'
    else:
        return 'unknown'

df_long['genome_type'] = df_long['sample'].apply(assign_genome_type)
## remove the sample with genome_type unknown
df_long = df_long[df_long['genome_type'] != 'unknown']

## collect all the sample names for each strain_name, and check whether it has plasmid genome_type
strain_to_samples = df_long.groupby('strain_name')['sample'].apply(list)
strain_has_plasmid = {}
for strain, samples in strain_to_samples.items():
    genome_types = df_long[df_long['sample'].isin(samples)]['genome_type'].unique()
    strain_has_plasmid[strain] = 'plasmid' in genome_types

# Filter df_long to keep only strains that have plasmid
strains_with_plasmid = [strain for strain, has_plasmid in strain_has_plasmid.items() if has_plasmid]

df_long = df_long[df_long['strain_name'].isin(strains_with_plasmid)]

## remove motifs that have all < 0.3 values across all samples
motif_max = df_long.groupby('motif_identifier')['fraction'].max()
motifs_to_keep = motif_max[motif_max >= 0.3].index
df_long = df_long[df_long['motif_identifier'].isin(motifs_to_keep)]



# Display the result
print(df_long.head(20))
print(f"\nShape: {df_long.shape}")
print(f"\nUnique strains: {df_long['strain_name'].nunique()}")
print(f"\nGenome type distribution:")
print(df_long['genome_type'].value_counts())

# Save the long format data
output_file = motif_profile.replace('.csv', '_long.csv')
df_long.to_csv(output_file, index=False)
print(f"\nSaved to: {output_file}")





