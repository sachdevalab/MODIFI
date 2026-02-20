import pandas as pd

# Set pandas to display all rows
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

RM_file = "/home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/RM_systems/all_ctgs_RM.rm.genes.tsv"

# Read the file (tab-separated)
df = pd.read_csv(RM_file, sep='\t')

# Check if there's a contig column and filter if needed
if 'contig' in df.columns:
    # Filter for specific contig if needed
    df = df[df['contig'] == 'soil_1_332_C']
    print(f"Filtered for contig: soil_1_332_C\n")
elif 'Gene' in df.columns:
    # Try to extract contig from Gene column
    contig_name = 'soil_1_332_C'
    df = df[df['Gene'].str.startswith(contig_name)]
    print(f"Filtered for genes starting with: {contig_name}\n")

# Filter for MTases (Gene type = 'MT' or 'IIG')
mtases = df[df['Gene type'].isin(['MT', 'IIG'] )]
## unique by 'Operon'
mtases = mtases.drop_duplicates(subset=['Operon'])
## output all df
print("Full DataFrame:")
print(mtases)
# Count total MTases
total_mtases = len(mtases)

# Count MTases in operons (RM Operon) vs orphan (Singleton)
# Determine orphan status from the Operon column
operon_mtases = len(mtases[mtases['Operon'].str.contains('RM Operon', na=False)])
orphan_mtases = len(mtases[mtases['Operon'].str.contains('Singleton', na=False)])

# Calculate proportions
operon_proportion = operon_mtases / total_mtases if total_mtases > 0 else 0
orphan_proportion = orphan_mtases / total_mtases if total_mtases > 0 else 0

# Print results
print(f"Total MTases: {total_mtases}")
print(f"MTases in operons: {operon_mtases} ({operon_proportion:.2%})")
print(f"Orphan MTases: {orphan_mtases} ({orphan_proportion:.2%})")

# Show breakdown by operon vs singleton
print("\n--- Breakdown by Operon/Singleton ---")
is_singleton = mtases['Operon'].str.contains('Singleton', na=False)
print(f"In operons (RM Operon): {(~is_singleton).sum()}")
print(f"Singletons: {is_singleton.sum()}")

## count how many unique HMM are there
unique_hmms = mtases['Operon'].nunique()
print(f"\nUnique Operon among MTases: {unique_hmms}")
