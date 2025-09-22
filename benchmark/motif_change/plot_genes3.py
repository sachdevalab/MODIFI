import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow





prefix_1 = "infant_14"
genes_1 = [
    {"name": "RE", "start": 1001, "end": 3988, "strand": "+"},
    {"name": "MT", "start": 4011, "end": 5603, "strand": "+"},
    {"name": "TRD2'", "start": 5600, "end": 6763, "strand": "+"},
    {"name": "TRD1'", "start": 7814, "end": 8971, "strand": "-"},
]

prefix_2 = "infant_26"
genes_2 = [
    {"name": "RE", "start": 1001, "end": 3988, "strand": "+"},
    {"name": "MT", "start": 4011, "end": 5603, "strand": "+"},
    {"name": "TRD1''", "start": 5600, "end": 6817, "strand": "+"},
    {"name": "TRD2''", "start": 7868, "end": 8971, "strand": "-"},
]


# Highlight inversion region
inversion_start = 6122
inversion_end = 8510


# Create figure with vertical layout
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 6), sharex=True, gridspec_kw={'hspace': 0.1})

# Calculate region bounds from both gene sets
all_genes = genes_1 + genes_2
region_start = min(g["start"] for g in all_genes) - 1000
region_end = max(g["end"] for g in all_genes) + 1000

# Plot colors for genes
colors = ["skyblue", "lightgreen", "salmon", "salmon"]

def plot_genes(ax, genes, title, y_center=1):
    # Draw genome axis
    ax.hlines(y_center, region_start, region_end, color="black", linewidth=2)
    
    # Plot each gene as arrow
    for i, gene in enumerate(genes):
        start, end, strand = gene["start"], gene["end"], gene["strand"]
        arrow_length = end - start
        
        # Adjust head_length to be proportional to arrow length
        head_len = min(300, arrow_length * 0.1)
        
        if strand == "+":
            # Forward strand: arrow from start to end
            ax.add_patch(FancyArrow(start, y_center, arrow_length, 0,
                                    width=0.15, head_width=0.3, head_length=head_len,
                                    length_includes_head=True, color=colors[i % len(colors)]))
        else:
            # Reverse strand: arrow from end to start (pointing backwards)
            ax.add_patch(FancyArrow(end, y_center, -arrow_length, 0,
                                    width=0.15, head_width=0.3, head_length=head_len,
                                    length_includes_head=True, color=colors[i % len(colors)]))
        
        # Label
        ax.text((start+end)//2, y_center + 0.4, gene["name"], ha="center", fontsize=11)
    
    # Formatting
    ax.set_xlim(region_start, region_end)
    ax.set_ylim(-0.3, 2.3)  # Extended range to accommodate longer lines
    ax.set_title(title, fontsize=12, pad=10)
    ax.axis("off")

# Plot both gene sets
plot_genes(ax1, genes_1, prefix_1)
plot_genes(ax2, genes_2, prefix_2)

# Add shared inversion region to both subplots
ax1.axvspan(inversion_start, inversion_end, color='orange', alpha=0.15)
ax2.axvspan(inversion_start, inversion_end, color='orange', alpha=0.15)

# # Add annotation for the shared inversion
# fig.text(0.5, 0.5, 'Shared Inversion Region', ha='center', va='center', 
#          fontsize=11, color='orange', fontweight='bold', rotation=90,
#          bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='orange', alpha=0.8))

plt.tight_layout()
plt.show()
## save the figure
fig.savefig(f"../../tmp/results2/gene_arrows_vertical_comparison.png", bbox_inches="tight", dpi=300)

