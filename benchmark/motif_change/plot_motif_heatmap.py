import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# ...existing code...




def heatmap(df, heat_map):
    # Prepare data for heatmap
    # Set 'anno' as index and extract motif columns
    motif_columns = [col for col in df.columns if col not in ['sample', 'type', 'anno', 'contig']]
    
    # Sort data by sample and type order (host, P1, P2)
    def get_sort_order(row):
        sample = row['sample']
        anno = row['anno']
        
        # Determine type order: host=0, P1=1, P2=2
        if 'host' in row['type']:
            type_order = 0
        elif 'P1' in anno:
            type_order = 1
        elif 'P2' in anno:
            type_order = 2
        else:
            type_order = 3
        
        return (sample, type_order)
    
    # Sort the dataframe
    df_sorted = df.copy()
    df_sorted['sort_key'] = df_sorted.apply(get_sort_order, axis=1)
    df_sorted = df_sorted.sort_values('sort_key')
    df_sorted = df_sorted.drop('sort_key', axis=1)
    
    # Create the heatmap data using 'anno' as index
    heatmap_data = df_sorted.set_index('anno')[motif_columns]
    
    # Create row colors based on type
    def get_base_type(type_str):
        if 'plasmid' in type_str:
            return 'plasmid'
        elif 'host' in type_str:
            return 'host'
        else:
            return 'unknown'
    
    df_indexed = df_sorted.set_index('anno')
    type_labels = df_indexed['type'].apply(get_base_type)
    sample_labels = df_indexed['sample']
    
    # Map types to colors
    type_palette = {'plasmid': '#1f77b4', 'host': '#ff7f0e', 'unknown': '#bbbbbb'}
    type_colors = type_labels.map(type_palette)
    
    # Map samples to colors using a colormap
    unique_samples = sample_labels.unique()
    sample_palette = plt.cm.Set3(range(len(unique_samples)))
    sample_color_map = dict(zip(unique_samples, sample_palette))
    sample_colors = sample_labels.map(sample_color_map)
    
    # Combine both color bars into a DataFrame
    row_colors = pd.DataFrame({
        'Type': type_colors,
        'Sample': sample_colors
    }, index=df_indexed.index)

    # Create custom legend patches for both type and sample
    type_handles = [
        mpatches.Patch(color='#1f77b4', label='plasmid'),
        mpatches.Patch(color='#ff7f0e', label='host'),
        # mpatches.Patch(color='#bbbbbb', label='unknown')
    ]
    
    # Create sample legend patches (limit to top samples to avoid clutter)
    sample_handles = []
    for sample in unique_samples[:10]:  # Show only top 10 samples
        sample_handles.append(mpatches.Patch(color=sample_color_map[sample], label=sample))
    
    if len(unique_samples) > 10:
        sample_handles.append(mpatches.Patch(color='lightgray', label=f'... +{len(unique_samples)-10} more samples'))
    
    # Plot clustered heatmap without row clustering to maintain sample grouping
    g = sns.clustermap(
        heatmap_data,
        row_colors=row_colors,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(6, 5),
        row_cluster=False,  # Don't cluster rows to keep sample grouping
        col_cluster=True,   # Still cluster columns (motifs)
        colors_ratio=0.03,  # Minimize the width of the row color bar
        yticklabels=True,   # Show all row labels (anno)
        xticklabels=True    # Show all column labels (motifs)
    )

    # Rotate x-axis labels for better readability
    g.ax_heatmap.tick_params(axis='x', rotation=90, labelsize=8)
    g.ax_heatmap.tick_params(axis='y', labelsize=10)

    # Add type legend outside the main figure on the left side
    # Type legend
    g.ax_heatmap.legend(
        handles=type_handles,
        title="Type",
        loc='center right',
        bbox_to_anchor=(-0.15, 0.5),  # Position to the left of the plot, centered vertically
        borderaxespad=0.,
        frameon=True,
        fancybox=False,
        shadow=False
    )

    plt.savefig(heat_map, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()



# motif_profile = "/home/shuaiw/borg/paper/linkage/m64004_210929_143746.p100/motif_profile.csv"
motif_profile = "180_4.csv"
heat_map = "../../tmp/results2/motif_change_heatmap.png"


df = pd.read_csv(motif_profile)
print("Data shape:", df.shape)
print("Columns:", df.columns.tolist())
print("Sample data:")
print(df.head())
heatmap(df, heat_map)