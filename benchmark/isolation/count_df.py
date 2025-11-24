import pandas as pd

all_dir = "/home/shuaiw/borg/paper/isolation//batch2_results/"
fig_dir = "../../tmp/figures/motif_sharing/"  
    
df = pd.read_csv(f"{fig_dir}/jaccard_all_samples.csv")

## only retain the df with all_num >= 2
df = df[df['all_num'] >= 2]

print(f"Total rows in dataset: {len(df)}")
print(f"Unique relations: {df['relation'].unique()}")
print()

# Get top 10 rows with highest jaccard_similarity for each relation
for relation in df['relation'].unique():
    print(f"=== Top 10 highest Jaccard similarity for relation: {relation} ===")
    relation_df = df[df['relation'] == relation]
    top_10 = relation_df.nlargest(10, 'jaccard_similarity')
    
    print(f"Total entries for {relation}: {len(relation_df)}")
    print("Top 10:")
    for idx, row in top_10.iterrows():
        print(f"  {row['mge_contig']} -> {row['host_contig']} | "
              f"Jaccard: {row['jaccard_similarity']:.4f} | "
              f"Share: {row['share_num']}/{row['all_num']} | "
              f"Type: {row['mge_type']}")
    print()