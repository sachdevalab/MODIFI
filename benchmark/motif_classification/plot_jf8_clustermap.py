#!/usr/bin/env python3
"""
JF8 (soil_1 control DB) motif x contig clustermap, matched to the paper's panel-b style:
 - columns grouped by species in the panel-b left-to-right order (no column dendrogram);
   within each species block contigs are sorted by descending total modification signal
 - species color strip (col_colors) + italic legend with one-letter genus, in panel-b legend order
 - each motif row label bolds its modified base (6mA A) via mathtext
 - viridis heatmap (kept from the previous python plot)
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns

DATADIR = "/home/shuaiw/borg/revision/motif_class"
FIGDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"
os.makedirs(FIGDIR, exist_ok=True)

# motif x contig matrix (rows = motif, cols = contig)
df = pd.read_csv(os.path.join(DATADIR, "jf8_heatmap_matrix.csv"), index_col=0)
ann = pd.read_csv(os.path.join(DATADIR, "jf8_heatmap_colann.csv"))
contig2sp = dict(zip(ann["contig"], ann["species"]))
cu = pd.read_csv(os.path.join(DATADIR, "jf8_contig_unitig.tsv"), sep="\t",
                 header=None, names=["contig", "unitig"])
contig2unitig = dict(zip(cu["contig"], cu["unitig"]))

# ---- species: column order, legend order, colors, one-letter-genus labels (panel b) ----
col_species_order = ["Ruminococcus gnavus", "Bacteroides vulgatus", "Collinsella aerofaciens",
                     "Bacteroides caccae", "Clostridium bolteae", "Escherichia coli",
                     "Bacteroides ovatus", "Bacteroides thetaiotaomicron"]
legend_order = ["Ruminococcus gnavus", "Escherichia coli", "Clostridium bolteae",
                "Collinsella aerofaciens", "Bacteroides vulgatus",
                "Bacteroides thetaiotaomicron", "Bacteroides ovatus", "Bacteroides caccae"]
sp_colors = {
    "Ruminococcus gnavus": "#E41A1C",           # red
    "Escherichia coli": "#F0C000",              # gold
    "Clostridium bolteae": "#7FE000",           # bright green
    "Collinsella aerofaciens": "#1B7837",       # dark green
    "Bacteroides vulgatus": "#00A0A0",          # teal
    "Bacteroides thetaiotaomicron": "#5B9BD5",  # light blue
    "Bacteroides ovatus": "#7030A0",            # purple
    "Bacteroides caccae": "#000000",            # black
}
sp_label = {
    "Ruminococcus gnavus": "R. gnavus", "Escherichia coli": "E. coli",
    "Clostridium bolteae": "C. bolteae", "Collinsella aerofaciens": "C. aerofaciens",
    "Bacteroides vulgatus": "B. vulgatus", "Bacteroides thetaiotaomicron": "B. theta",
    "Bacteroides ovatus": "B. ovatus", "Bacteroides caccae": "B. caccae",
}

# ---- modified-base (6mA A) position (1-based) per ground-truth motif, for bold labels ----
mod_pos = {
    "GATC": 2, "AGATCC": 3, "GGATCT": 3, "AGATCT": 3, "AATCC": 2,
    "CCANNNNNNCAT": 3, "ATGNNNNNNTGG": 1, "CAGNNNNNGGA": 2, "CCATC": 3, "GATGG": 2,
    "TCACNNNNNATG": 3, "GCACNNNNNNGTT": 3, "AACNNNNNNGTGC": 1, "GGAGC": 3, "CAGGAG": 5, "GAGC": 2,
}

def bold_label(motif):
    p = mod_pos.get(motif)
    parts = []
    for i, ch in enumerate(motif):
        parts.append(r"\mathbf{%s}" % ch if (p and i == p - 1) else r"\mathrm{%s}" % ch)
    return "$" + "".join(parts) + "$"

# ---- fixed motif (row) order matching the panel-b screenshot (GAGC at top ...) ----
# panel-b order = reverse of the J8_paper_motifs.table order
paper_order = ["GATC", "AGATCC", "GGATCT", "AGATCT", "AATCC", "CCANNNNNNCAT", "ATGNNNNNNTGG",
               "CAGNNNNNGGA", "CCATC", "GATGG", "TCACNNNNNATG", "GCACNNNNNNGTT",
               "AACNNNNNNGTGC", "GGAGC", "CAGGAG", "GAGC"]
motif_order = [m for m in reversed(paper_order) if m in df.index]
df = df.loc[motif_order]

# ---- order columns: group by species (panel-b order), sort by signal within each group ----
signal = df.sum(axis=0)                       # total modification signal per contig
ordered_cols = []
for sp in col_species_order:
    members = [c for c in df.columns if contig2sp.get(c) == sp]
    members.sort(key=lambda c: signal[c], reverse=True)
    ordered_cols.extend(members)
df = df[ordered_cols]

df_plot = df.copy()
df_plot.columns = [contig2unitig.get(c, c) for c in df.columns]
col_colors = pd.Series({contig2unitig.get(c, c): sp_colors[contig2sp[c]] for c in df.columns},
                       name="")

g = sns.clustermap(
    df_plot,
    row_cluster=True, col_cluster=True,       # cluster both motifs and contigs (with dendrograms)
    cmap="viridis",
    figsize=(12, 5.5),
    col_colors=col_colors,
    xticklabels=True, yticklabels=True,
    cbar_kws={"label": "Modification score"},
    dendrogram_ratio=(0.06, 0.10),
    colors_ratio=0.03,
)
g.ax_heatmap.set_xlabel("")
g.ax_heatmap.set_ylabel("")
g.ax_heatmap.tick_params(axis="y", labelsize=11)
g.ax_heatmap.tick_params(axis="x", labelsize=10)
for lab in g.ax_heatmap.get_xticklabels():
    lab.set_rotation(90)
# bold the modified base in each motif row label (clustered order)
g.ax_heatmap.set_yticklabels([bold_label(t.get_text()) for t in g.ax_heatmap.get_yticklabels()])

# open a right margin; shrink the column-aligned axes together so the species strip stays aligned
SHRINK = 0.86
axes = [g.ax_heatmap, g.ax_col_colors]
if getattr(g, "ax_col_dendrogram", None) is not None:
    axes.append(g.ax_col_dendrogram)
for ax in axes:
    b = ax.get_position()
    ax.set_position([b.x0, b.y0, b.width * SHRINK, b.height])

# species legend (panel-b order, one-letter genus, italic), right margin, upper
handles = [Patch(facecolor=sp_colors[s], edgecolor="none", label=sp_label[s]) for s in legend_order]
leg = g.figure.legend(handles=handles, title="Species",
                      bbox_to_anchor=(0.99, 0.90), loc="upper right",
                      frameon=False, fontsize=9, title_fontsize=10)
for t in leg.get_texts():
    t.set_fontstyle("italic")

# colorbar in the right margin, below the legend
g.cax.set_position([0.90, 0.34, 0.012, 0.18])
g.cax.yaxis.set_label_position("left")
g.cax.yaxis.set_ticks_position("right")

g.savefig(os.path.join(FIGDIR, "jf8_clustermap_soil1.pdf"), bbox_inches="tight")
g.savefig(os.path.join(FIGDIR, "jf8_clustermap_soil1.png"), dpi=200, bbox_inches="tight")

# source data next to the figure (columns in the plotted order)
df.to_csv(os.path.join(FIGDIR, "jf8_clustermap_soil1_sourcedata.csv"))
ann.to_csv(os.path.join(FIGDIR, "jf8_clustermap_soil1_species.csv"), index=False)
print("wrote", os.path.join(FIGDIR, "jf8_clustermap_soil1.pdf"), "and .png")
