import pandas as pd
import sys
from itertools import combinations

def count_shared_motifs(csv_file, fraction_cutoff=0.3):
    """
    Count shared motifs between contigs with a binary cutoff.
    
    Args:
        csv_file: Path to CSV file with columns: contig, motifString, fraction, BORG_Ref, Genome
        fraction_cutoff: Minimum fraction to consider a motif present (default: 0.3)
    
    Returns:
        DataFrame with contig pairs and their shared motif counts
    """
    # Read the data
    df = pd.read_csv(csv_file)
    
    # Filter motifs by fraction cutoff (binary: present if >= cutoff)
    df_present = df[df['fraction'] >= fraction_cutoff].copy()
    
    # Get unique motif sets for each contig
    contig_motifs = {}
    for contig in df_present['contig'].unique():
        contig_data = df_present[df_present['contig'] == contig]
        # Store unique motifs for this contig
        contig_motifs[contig] = set(contig_data['motifString'].unique())
    
    # Count shared motifs for all contig pairs
    results = []
    contigs = list(contig_motifs.keys())
    
    for contig1, contig2 in combinations(contigs, 2):
        shared_motifs = contig_motifs[contig1] & contig_motifs[contig2]
        shared_count = len(shared_motifs)
        
        results.append({
            'Contig1': contig1,
            'Contig2': contig2,
            'Shared_Motifs': shared_count,
            'Contig1_Total': len(contig_motifs[contig1]),
            'Contig2_Total': len(contig_motifs[contig2]),
            'Shared_Motifs_List': ','.join(sorted(shared_motifs)) if shared_motifs else ''
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Shared_Motifs', ascending=False)
    
    return results_df, contig_motifs

def count_shared_for_pair(csv_file, contig1, contig2, fraction_cutoff=0.3):
    """
    Count shared motifs between two specific contigs.
    """
    df = pd.read_csv(csv_file)
    
    # Filter by fraction cutoff
    df_present = df[df['fraction'] >= fraction_cutoff]
    
    # Get motifs for each contig
    motifs1 = set(df_present[df_present['contig'] == contig1]['motifString'].unique())
    motifs2 = set(df_present[df_present['contig'] == contig2]['motifString'].unique())
    
    # Find shared motifs
    shared = motifs1 & motifs2
    
    print(f"\n=== Shared Motifs between {contig1} and {contig2} ===")
    print(f"Fraction cutoff: >= {fraction_cutoff}")
    print(f"{contig1}: {len(motifs1)} motifs")
    print(f"{contig2}: {len(motifs2)} motifs")
    print(f"Shared: {len(shared)} motifs")
    
    if shared:
        print(f"\nShared motifs: {', '.join(sorted(shared))}")
    
    return len(shared), motifs1, motifs2, shared

if __name__ == "__main__":

    
    
    all_dir = "/home/shuaiw/borg/paper/gg_run3/"
    seq_dir = "/home/shuaiw/borg/paper/borg_data/profile5/"
    cluster = "profile"

    cluster_species = "borg"
    csv_file = f"{seq_dir}/{cluster}_profile_df_filtered.csv"
    
    # Read original data to get contig categories
    df_orig = pd.read_csv(csv_file)
    
    # Get contig categorization based on BORG_Ref or Genome column
    # Assuming BORG_Ref or similar column indicates whether contig is borg, Mp, or no-Mp
    contig_category = {}
    for contig in df_orig['contig'].unique():
        contig_data = df_orig[df_orig['contig'] == contig].iloc[0]
        # Determine category - adjust logic based on actual data structure
        if 'BORG_Ref' in df_orig.columns:
            if pd.notna(contig_data['BORG_Ref']) and 'borg' in str(contig_data['BORG_Ref']).lower():
                contig_category[contig] = 'borg'
            elif 'Genome' in df_orig.columns and pd.notna(contig_data['Genome']):
                # Check if this is Mp or no-Mp based on genome name or other criteria
                genome_name = str(contig_data['Genome'])
                if 'Mp' in genome_name or 'methylase' in genome_name.lower():
                    contig_category[contig] = 'Mp'
                else:
                    contig_category[contig] = 'no-Mp'
            else:
                contig_category[contig] = 'other'
        elif 'Genome' in df_orig.columns:
            genome_name = str(contig_data['Genome'])
            if 'borg' in genome_name.lower():
                contig_category[contig] = 'borg'
            elif 'Mp' in genome_name:
                contig_category[contig] = 'Mp'
            else:
                contig_category[contig] = 'no-Mp'
        else:
            # Fallback: try to infer from contig name itself
            if 'borg' in contig.lower():
                contig_category[contig] = 'borg'
            elif 'mp' in contig.lower():
                contig_category[contig] = 'Mp'
            else:
                contig_category[contig] = 'no-Mp'
    
    results_df, contig_motifs = count_shared_motifs(csv_file)
    
    # Add categories to results
    results_df['Contig1_Category'] = results_df['Contig1'].map(contig_category)
    results_df['Contig2_Category'] = results_df['Contig2'].map(contig_category)
    
    # Create pair type column
    def classify_pair(row):
        cat1, cat2 = sorted([row['Contig1_Category'], row['Contig2_Category']])
        if (cat1 == 'borg' and cat2 == 'Mp') or (cat1 == 'Mp' and cat2 == 'borg'):
            return 'borg-Mp'
        elif (cat1 == 'borg' and cat2 == 'no-Mp') or (cat1 == 'no-Mp' and cat2 == 'borg'):
            return 'borg-noMp'
        else:
            return 'other'
    
    results_df['Pair_Type'] = results_df.apply(classify_pair, axis=1)
    
    # Filter for borg-Mp and borg-noMp pairs
    filtered_df = results_df[results_df['Pair_Type'].isin(['borg-Mp', 'borg-noMp'])].copy()
    
    # Display results with all columns including shared motifs
    print("\n=== Shared Motifs Count for All Contig Pairs ===")
    print(f"Total pairs: {len(results_df)}")
    print(f"borg-Mp pairs: {len(filtered_df[filtered_df['Pair_Type'] == 'borg-Mp'])}")
    print(f"borg-noMp pairs: {len(filtered_df[filtered_df['Pair_Type'] == 'borg-noMp'])}")
    print(f"\nResults (sorted by shared motif count):\n")
    
    # Set pandas display options to show full content
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.width', None)
    
    # Save to CSV file
    output_file = f"{seq_dir}/{cluster}_shared_motifs.csv"
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    # Print summary statistics by pair type
    print("\n=== Summary Statistics by Pair Type ===")
    for pair_type in ['borg-Mp', 'borg-noMp']:
        subset = filtered_df[filtered_df['Pair_Type'] == pair_type]
        if len(subset) > 0:
            print(f"\n{pair_type}:")
            print(f"  Count: {len(subset)}")
            print(f"  Mean: {subset['Shared_Motifs'].mean():.2f}")
            print(f"  Median: {subset['Shared_Motifs'].median():.2f}")
            print(f"  Min: {subset['Shared_Motifs'].min()}")
            print(f"  Max: {subset['Shared_Motifs'].max()}")
    
    # Create box plot
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=filtered_df, x='Pair_Type', y='Shared_Motifs', palette='Set2')
    plt.xlabel('Contig Pair Type', fontsize=12)
    plt.ylabel('Number of Shared Motifs', fontsize=12)
    plt.title('Shared Motif Distribution: borg-Mp vs borg-noMp', fontsize=14)
    plt.tight_layout()
    
    plot_file = f"{seq_dir}/{cluster}_shared_motifs_boxplot.pdf"
    plt.savefig(plot_file, dpi=300)
    print(f"\nBox plot saved to: {plot_file}")