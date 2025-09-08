import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

# # Example genes from Streptococcus pneumoniae D39 (approximate coordinates)
genes = [
    {"name": "RE", "start": 2319186, "end": 2322173, "strand": "+"},
    {"name": "MT", "start": 2322196, "end": 2323788, "strand": "+"},
    {"name": "S1", "start": 2323785, "end": 2324948, "strand": "+"},
    {"name": "S2", "start": 2325999, "end": 2327159, "strand": "+"},
]

# Highlight inversion region
inversion_start = 2323856
inversion_end = 2327115


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
ax.hlines(1, region_start, region_end, color="black")

# Plot each gene as arrow
colors = ["skyblue", "lightgreen", "salmon"]  # different colors
for i, gene in enumerate(genes):
    start, end, strand = gene["start"], gene["end"], gene["strand"]
    arrow_length = end - start
    direction = 1 if strand == "+" else -1
    ax.add_patch(FancyArrow(start, 1, arrow_length*direction, 0,
                            width=0.2, head_width=0.5, head_length=500,
                            length_includes_head=True, color=colors[i % len(colors)]))
    # Label
    ax.text((start+end)//2, 1.3, gene["name"], ha="center", fontsize=12)


# ax.axvspan(inversion_start, inversion_end, color='orange', alpha=0.3, label='Inversion region')
# ax.annotate('Inversion', xy=((inversion_start+inversion_end)//2, 0.95), xycoords=('data', 'axes fraction'),
#             ha='center', va='top', fontsize=10, color='orange', fontweight='bold',
#             bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='orange', alpha=0.5))

# Formatting
ax.set_xlim(region_start, region_end)
ax.set_ylim(0, 2)
ax.axis("off")
handles, labels = ax.get_legend_handles_labels()
if 'Inversion region' not in labels:
    ax.legend(loc='upper right')

plt.show()
## save the figure
fig.savefig("../../tmp/results2/gene_arrows_infant2.png", bbox_inches="tight", dpi=300)

