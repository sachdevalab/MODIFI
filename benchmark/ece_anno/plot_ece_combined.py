#!/usr/bin/env python3
"""Combined 7-panel figure for the revised filter-passing ECE set (3,976), plasmid vs virus:
 a length dist  b depth dist  c environment dist  d linkage-score dist
 e length (score=0 vs >0)  f depth (score=0 vs >0)  g mod-density (score=0 vs >0)
Reads ece_master_revised.csv. Outputs fig + sourcedata to tmp/rev_figs/ece_anno/."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

A = "/home/shuaiw/borg/revision/ece_anno"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
os.makedirs(OUTDIR, exist_ok=True)
CP, CV = "#4C78A8", "#E45756"                 # plasmid, virus
TYPE_ORDER = ["plasmid", "virus"]
GRP_ORDER = ["score=0", "score>0"]

df = pd.read_csv(f"{A}/expanded/ece_master_revised.csv")
for c in ("length", "depth", "final_score", "mod_density_per_kb"):
    df[c] = pd.to_numeric(df[c], errors="coerce")
df = df[df["type"].isin(TYPE_ORDER)].copy()

sns.set_style("whitegrid")
plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold"})
fig, axes = plt.subplots(3, 3, figsize=(17, 15))
pal = {"plasmid": CP, "virus": CV}
grp_pal = {"score=0": "#BAB0AC", "score>0": "#59A14F"}

def panel(ax, letter, title):
    ax.set_title(f"{letter}. {title}", loc="left")

# a. length distribution (log x)
ax = axes[0, 0]; panel(ax, "a", "ECE length")
for t in TYPE_ORDER:
    v = df.loc[df.type == t, "length"].dropna()
    ax.hist(v, bins=np.logspace(np.log10(1000), np.log10(max(v.max(), 1e6)), 40),
            histtype="step", linewidth=2, color=pal[t], label=f"{t} (n={len(v)})")
ax.set_xscale("log"); ax.set_xlabel("length (bp)"); ax.set_ylabel("ECEs"); ax.legend()

# b. depth distribution (log x)
ax = axes[0, 1]; panel(ax, "b", "ECE read depth")
for t in TYPE_ORDER:
    v = df.loc[df.type == t, "depth"].dropna(); v = v[v > 0]
    ax.hist(v, bins=np.logspace(np.log10(max(v.min(), 1)), np.log10(v.max()), 40),
            histtype="step", linewidth=2, color=pal[t], label=t)
ax.set_xscale("log"); ax.set_xlabel("mean depth (x)"); ax.set_ylabel("ECEs"); ax.legend()

# c. environment distribution (grouped bars)
ax = axes[0, 2]; panel(ax, "c", "Environment")
ct = df.groupby(["environment", "type"]).size().unstack(fill_value=0)
ct = ct.reindex(columns=TYPE_ORDER, fill_value=0)
ct["tot"] = ct.sum(1); ct = ct.sort_values("tot", ascending=False).drop(columns="tot")
x = np.arange(len(ct)); w = 0.4
ax.bar(x - w/2, ct["plasmid"], w, color=CP, label="plasmid")
ax.bar(x + w/2, ct["virus"], w, color=CV, label="virus")
ax.set_xticks(x); ax.set_xticklabels(ct.index, rotation=45, ha="right")
ax.set_ylabel("ECEs"); ax.legend()

# d. linkage score distribution
ax = axes[1, 0]; panel(ax, "d", "Linkage score")
for t in TYPE_ORDER:
    v = df.loc[df.type == t, "final_score"].dropna()
    ax.hist(v, bins=np.linspace(0, 1, 26), histtype="step", linewidth=2, color=pal[t], label=t)
ax.set_xlabel("linkage final_score"); ax.set_ylabel("ECEs"); ax.set_yscale("log"); ax.legend()

# e/f/g: score=0 vs >0 boxplots, split by type
def grouped_box(ax, letter, title, col, logy, ylab):
    panel(ax, letter, title)
    d = df.dropna(subset=[col]).copy()
    sns.boxplot(data=d, x="type", y=col, hue="grp", order=TYPE_ORDER, hue_order=GRP_ORDER,
                palette=grp_pal, showfliers=False, ax=ax)
    if logy: ax.set_yscale("log")
    ax.set_xlabel(""); ax.set_ylabel(ylab); ax.legend(title="", loc="upper right")

grouped_box(axes[1, 1], "e", "Length: score=0 vs >0", "length", True, "length (bp)")
grouped_box(axes[1, 2], "f", "Depth: score=0 vs >0", "depth", True, "mean depth (x)")
grouped_box(axes[2, 0], "g", "Modification density: score=0 vs >0", "mod_density_per_kb", False,
            "modifications per kb")

# hide unused panels
for ax in (axes[2, 1], axes[2, 2]):
    ax.axis("off")

# summary text in the spare space
lines = [f"Revised filter-passing ECEs: {len(df)}  (plasmid {int((df.type=='plasmid').sum())}, "
         f"virus {int((df.type=='virus').sum())})",
         f"score>0: {int((df.grp=='score>0').sum())}   score=0: {int((df.grp=='score=0').sum())}"]
for t in TYPE_ORDER:
    for g in GRP_ORDER:
        s = df[(df.type == t) & (df.grp == g)]
        lines.append(f"{t} {g}: n={len(s)}  median len={s.length.median():.0f}  "
                     f"depth={s.depth.median():.1f}  moddens={s.mod_density_per_kb.median():.2f}/kb")
axes[2, 1].text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=10, family="monospace")

fig.suptitle("Metagenome ECEs (revised filter, 3,976): plasmid vs virus", fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.98])
for ext in ("pdf", "png"):
    fig.savefig(f"{OUTDIR}/ece_combined_panels.{ext}", bbox_inches="tight", dpi=200)
df.to_csv(f"{OUTDIR}/ece_combined_panels_sourcedata.csv", index=False)
print("wrote ece_combined_panels.pdf/.png + sourcedata")
print("\n".join(lines))
