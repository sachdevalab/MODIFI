#!/usr/bin/env python3
"""Step 6: 4-panel linkage-validation figure on the REVISED high-confidence network:
 a) MGE vs host GC (color=habitat, shape=MGE_type, Pearson r/p, lm)
 b) MGE vs host coverage log-log (r/p, lm)
 c) cosine similarity by habitat (boxplot)
 d) CRISPR consistent linkages per sample (stacked bar by MGE type)
Reads network_99_revised/mge_host_gc_cov.csv + consistent_linkages_by_sample.csv."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

D = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised"
df = pd.read_csv(f"{D}/mge_host_gc_cov.csv")
for c in ("MGE_gc", "host_gc", "MGE_cov", "host_cov", "cos_sim"):
    df[c] = pd.to_numeric(df[c], errors="coerce")

habitats = sorted(df["environment"].dropna().unique())
cmap = plt.get_cmap("tab10")
hcol = {h: cmap(i % 10) for i, h in enumerate(habitats)}
shape = {"plasmid": "o", "virus": "x"}

fig, axes = plt.subplots(2, 2, figsize=(15, 13))

def scatter_panel(ax, xcol, ycol, letter, title, log=False):
    d = df.dropna(subset=[xcol, ycol]).copy()
    if log:
        d = d[(d[xcol] > 0) & (d[ycol] > 0)]
    for t in ("plasmid", "virus"):
        for h in habitats:
            s = d[(d.MGE_type == t) & (d.environment == h)]
            if len(s):
                ax.scatter(s[xcol], s[ycol], marker=shape[t], s=26, color=hcol[h],
                           alpha=0.75, linewidths=1.1 if t == "virus" else 0)
    x, y = d[xcol].values, d[ycol].values
    if log:
        r, p = pearsonr(np.log10(x), np.log10(y)); ax.set_xscale("log"); ax.set_yscale("log")
    else:
        r, p = pearsonr(x, y)
    # lm line
    lx = np.log10(x) if log else x; ly = np.log10(y) if log else y
    b, a = np.polyfit(lx, ly, 1)
    xs = np.linspace(lx.min(), lx.max(), 100)
    ys = a + b * xs
    ax.plot(10**xs if log else xs, 10**ys if log else ys, "k-", lw=1.5)
    ax.text(0.05, 0.93, f"r = {r:.2f}\np = {p:.1e}", transform=ax.transAxes, va="top", fontsize=11)
    ax.set_title(f"{letter}. {title}", loc="left", fontweight="bold")

scatter_panel(axes[0, 0], "MGE_gc", "host_gc", "a", "MGE vs host GC content")
axes[0, 0].set_xlabel("MGE GC content"); axes[0, 0].set_ylabel("Host GC content")
scatter_panel(axes[0, 1], "MGE_cov", "host_cov", "b", "MGE vs host coverage", log=True)
axes[0, 1].set_xlabel("MGE coverage (log)"); axes[0, 1].set_ylabel("Host coverage (log)")

# legend (habitat colors + shape) on panel a
hleg = [plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=hcol[h], markersize=9, label=h) for h in habitats]
sleg = [plt.Line2D([0], [0], marker=shape[t], color="k", linestyle="", markersize=8, label=t) for t in ("plasmid", "virus")]
axes[0, 0].legend(handles=hleg + sleg, fontsize=8, loc="lower right", framealpha=0.9)

# c) cosine by habitat
ax = axes[1, 0]; ax.set_title("c. Cosine similarity by habitat", loc="left", fontweight="bold")
d = df.dropna(subset=["cos_sim"])
order = d.groupby("environment")["cos_sim"].median().sort_values(ascending=False).index.tolist()
data = [d[d.environment == h]["cos_sim"].values for h in order]
bp = ax.boxplot(data, patch_artist=True, showfliers=True, flierprops=dict(marker=".", markersize=3))
for patch, h in zip(bp["boxes"], order):
    patch.set_facecolor(hcol[h])
for med in bp["medians"]:
    med.set_color("black")
ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order, rotation=45, ha="right")
ax.set_ylabel("Cosine similarity")

# d) CRISPR consistent linkages per sample, stacked by MGE type
ax = axes[1, 1]; ax.set_title("d. CRISPR-consistent linkages", loc="left", fontweight="bold")
cf = f"{D}/consistent_linkages_by_sample.csv"
if os.path.exists(cf):
    c = pd.read_csv(cf)
    piv = c.pivot_table(index="Sample", columns="MGE Type", values="Consistent Linkages",
                        aggfunc="sum", fill_value=0)
    for t in ("plasmid", "virus"):
        if t not in piv.columns:
            piv[t] = 0
    piv["tot"] = piv["plasmid"] + piv["virus"]
    piv = piv.sort_values("tot", ascending=False).drop(columns="tot")
    x = np.arange(len(piv))
    ax.bar(x, piv["plasmid"], color="#F8766D", label="plasmid")
    ax.bar(x, piv["virus"], bottom=piv["plasmid"], color="#00BFC4", label="virus")
    ax.set_xticks(x); ax.set_xticklabels(piv.index, rotation=90, fontsize=7)
    ax.set_ylabel("Number of consistent linkages"); ax.legend(title="MGE Type")
else:
    ax.text(0.5, 0.5, "consistent_linkages_by_sample.csv not found", ha="center")

fig.suptitle("Linkage validation, revised high-confidence set", fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.98])
for ext in ("pdf", "png"):
    fig.savefig(f"{D}/linkage_validation_revised.{ext}", bbox_inches="tight", dpi=200)
df.to_csv(f"{D}/linkage_validation_revised_sourcedata.csv", index=False)
print("wrote linkage_validation_revised.pdf/.png + sourcedata")
