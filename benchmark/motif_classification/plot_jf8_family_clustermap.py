#!/usr/bin/env python3
"""
JF8 (soil_1 control DB) clustermap over the 17 family-collapsed motifs (fraction>=0.4 &
nDetected>=100 filter). Rows = family representatives (profile = max over member motifs);
columns = contigs (>100kb, assigned). Same style as plot_jf8_clustermap.py:
species color strip + one-letter-genus italic legend, bold modified base, viridis, both clustered.
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

df = pd.read_csv(os.path.join(DATADIR, "jf8_family_matrix.csv"), index_col=0)
ann = pd.read_csv(os.path.join(DATADIR, "jf8_heatmap_colann.csv"))
contig2sp = dict(zip(ann["contig"], ann["species"]))
cu = pd.read_csv(os.path.join(DATADIR, "jf8_contig_unitig.tsv"), sep="\t",
                 header=None, names=["contig", "unitig"])
contig2unitig = dict(zip(cu["contig"], cu["unitig"]))
bp = pd.read_csv(os.path.join(DATADIR, "jf8_family_boldpos.csv"))
mod_pos = dict(zip(bp["motif"], bp["centerPos"]))

legend_order = ["Ruminococcus gnavus", "Escherichia coli", "Clostridium bolteae",
                "Collinsella aerofaciens", "Bacteroides vulgatus",
                "Bacteroides thetaiotaomicron", "Bacteroides ovatus", "Bacteroides caccae"]
sp_colors = {
    "Ruminococcus gnavus": "#E41A1C", "Escherichia coli": "#F0C000",
    "Clostridium bolteae": "#7FE000", "Collinsella aerofaciens": "#1B7837",
    "Bacteroides vulgatus": "#00A0A0", "Bacteroides thetaiotaomicron": "#5B9BD5",
    "Bacteroides ovatus": "#7030A0", "Bacteroides caccae": "#000000",
}
sp_label = {
    "Ruminococcus gnavus": "R. gnavus", "Escherichia coli": "E. coli",
    "Clostridium bolteae": "C. bolteae", "Collinsella aerofaciens": "C. aerofaciens",
    "Bacteroides vulgatus": "B. vulgatus", "Bacteroides thetaiotaomicron": "B. theta",
    "Bacteroides ovatus": "B. ovatus", "Bacteroides caccae": "B. caccae",
}

def bold_label(motif):
    p = mod_pos.get(motif)
    return "$" + "".join((r"\mathbf{%s}" if (p and i == p - 1) else r"\mathrm{%s}") % ch
                         for i, ch in enumerate(motif)) + "$"

df_plot = df.copy()
df_plot.columns = [contig2unitig.get(c, c) for c in df.columns]
col_colors = pd.Series({contig2unitig.get(c, c): sp_colors[contig2sp[c]] for c in df.columns},
                       name="")

g = sns.clustermap(
    df_plot, row_cluster=True, col_cluster=True, cmap="viridis", figsize=(12, 5.8),
    col_colors=col_colors, xticklabels=True, yticklabels=True,
    cbar_kws={"label": "Modification score"},
    dendrogram_ratio=(0.06, 0.10), colors_ratio=0.03,
)
g.ax_heatmap.set_xlabel(""); g.ax_heatmap.set_ylabel("")
g.ax_heatmap.tick_params(axis="y", labelsize=10)
g.ax_heatmap.tick_params(axis="x", labelsize=10)
for lab in g.ax_heatmap.get_xticklabels():
    lab.set_rotation(90)
g.ax_heatmap.set_yticklabels([bold_label(t.get_text()) for t in g.ax_heatmap.get_yticklabels()])

SHRINK = 0.86
axes = [g.ax_heatmap, g.ax_col_colors]
if getattr(g, "ax_col_dendrogram", None) is not None:
    axes.append(g.ax_col_dendrogram)
for ax in axes:
    b = ax.get_position()
    ax.set_position([b.x0, b.y0, b.width * SHRINK, b.height])

handles = [Patch(facecolor=sp_colors[s], edgecolor="none", label=sp_label[s]) for s in legend_order]
leg = g.figure.legend(handles=handles, title="Species", bbox_to_anchor=(0.99, 0.90),
                      loc="upper right", frameon=False, fontsize=9, title_fontsize=10)
for t in leg.get_texts():
    t.set_fontstyle("italic")

g.cax.set_position([0.90, 0.34, 0.012, 0.18])
g.cax.yaxis.set_label_position("left")
g.cax.yaxis.set_ticks_position("right")

g.savefig(os.path.join(FIGDIR, "jf8_family_clustermap_soil1.pdf"), bbox_inches="tight")
g.savefig(os.path.join(FIGDIR, "jf8_family_clustermap_soil1.png"), dpi=200, bbox_inches="tight")
df.to_csv(os.path.join(FIGDIR, "jf8_family_clustermap_soil1_sourcedata.csv"))
print("wrote", os.path.join(FIGDIR, "jf8_family_clustermap_soil1.png"))
