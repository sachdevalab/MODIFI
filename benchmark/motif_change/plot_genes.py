import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow



prefix = "infant_2"
genes = [
    {"name": "RE", "start": 1001, "end": 3988, "strand": "+"},
    {"name": "MT", "start": 4011, "end": 5603, "strand": "+"},
    {"name": "TRD1", "start": 5600, "end": 6787, "strand": "+"},
    {"name": "TRD2", "start": 7838, "end": 8971, "strand": "-"},
]

prefix = "infant_14"
genes = [
    {"name": "RE", "start": 1001, "end": 3988, "strand": "+"},
    {"name": "MT", "start": 4011, "end": 5603, "strand": "+"},
    {"name": "TRD2", "start": 5600, "end": 6763, "strand": "+"},
    {"name": "TRD1", "start": 7814, "end": 8971, "strand": "-"},
]


# Highlight inversion region
inversion_start = 5671
inversion_end = 8929

# # # Example genes from Streptococcus pneumoniae D39 (approximate coordinates)
# genes = [
#     {"name": "RE", "start": 2319186, "end": 2322173, "strand": "+"},
#     {"name": "MT", "start": 2322196, "end": 2323788, "strand": "+"},
#     {"name": "S1", "start": 2323785, "end": 2324948, "strand": "+"},
#     {"name": "S2", "start": 2325999, "end": 2327159, "strand": "+"},
# ]

# # Highlight inversion region
# inversion_start = 2323856
# inversion_end = 2327115


# ## infant 2 3_C
# genes = [
#     {"name": "RE", "start": 2889998, "end": 2892985, "strand": "+"},
#     {"name": "MT", "start": 2893008, "end": 2894600, "strand": "+"},
#     {"name": "S1", "start": 2894597, "end": 2895784, "strand": "+"},
#     {"name": "S2", "start": 2896835, "end": 2897902, "strand": "+"},
# ]

# # Highlight inversion region
# inversion_start = 2894667
# inversion_end = 2897926


# genes = [
#     {"name": "RE", "start": 2346891, "end": 2349878, "strand": "+"},
#     {"name": "MT", "start": 2349901, "end": 2351493, "strand": "+"},
#     {"name": "S1", "start": 2351490, "end": 2352707, "strand": "+"},
#     {"name": "S2", "start": 2353758, "end": 2354864, "strand": "+"},
# ]

# # Highlight inversion region
# inversion_start = 2352013
# inversion_end = 2354401

# Create figure
fig, ax = plt.subplots(figsize=(10, 2))

# Draw genome axis
region_start = min(g["start"] for g in genes) - 1000
region_end   = max(g["end"] for g in genes) + 1000
ax.hlines(1, region_start, region_end, color="black", linewidth=2)

# Plot each gene as arrow
colors = ["skyblue", "lightgreen", "salmon", "salmon"]  # different colors
y_positions = [1, 1, 1.2, 0.8]  # Different y positions to avoid overlap for S genes

for i, gene in enumerate(genes):
    start, end, strand = gene["start"], gene["end"], gene["strand"]
    arrow_length = end - start
    y_pos = y_positions[i]
    
    # Adjust head_length to be proportional to arrow length but not too large
    head_len = min(300, arrow_length * 0.1)
    
    if strand == "+":
        # Forward strand: arrow from start to end
        ax.add_patch(FancyArrow(start, y_pos, arrow_length, 0,
                                width=0.15, head_width=0.3, head_length=head_len,
                                length_includes_head=True, color=colors[i % len(colors)]))
    else:
        # Reverse strand: arrow from end to start (pointing backwards)
        ax.add_patch(FancyArrow(end, y_pos, -arrow_length, 0,
                                width=0.15, head_width=0.3, head_length=head_len,
                                length_includes_head=True, color=colors[i % len(colors)]))
    
    # Label
    ax.text((start+end)//2, y_pos + 0.4, gene["name"], ha="center", fontsize=12)


ax.axvspan(inversion_start, inversion_end, color='orange', alpha=0.3, label='Inversion region')
ax.annotate('Inversion', xy=((inversion_start+inversion_end)//2, 0.90), xycoords=('data', 'axes fraction'),
            ha='center', va='top', fontsize=10, color='orange', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='orange', alpha=0.5))

# Formatting
ax.set_xlim(region_start, region_end)
ax.set_ylim(0, 2.5)  # Increased to accommodate different y-positions
ax.axis("off")
handles, labels = ax.get_legend_handles_labels()
if 'Inversion region' not in labels:
    ax.legend(loc='upper right')

plt.show()
## save the figure
fig.savefig(f"../../tmp/results2/gene_arrows_{prefix}.png", bbox_inches="tight", dpi=300)

