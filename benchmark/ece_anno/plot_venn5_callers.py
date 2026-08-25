#!/usr/bin/env python3
"""5-set Venn of UNFILTERED ECE detections across geNomad, VirSorter1, VirSorter2, VIBRANT, PlasX,
from the combined metagenome ECE table (>=5x, >=1kb contigs; 64 metagenomes).

geNomad = genomad_virus OR genomad_plasmid. Each ECE is keyed sample|MGE (unique across samples).
Outputs (next to the figure, per source-data convention):
  tmp/rev_figs/ece_anno/venn5_callers.pdf / .png
  tmp/rev_figs/ece_anno/venn5_callers_sourcedata.csv   (count for every non-empty region)
"""
import os
from itertools import combinations
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from venn import venn

C = "/home/shuaiw/borg/revision/ece_anno/expanded/combined_metagenome_eces.csv"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(C)
df["key"] = df["sample"] + "|" + df["MGE"].astype(str)
m = df["methods"].fillna("")

sets = {
    "geNomad":     set(df.loc[m.str.contains("genomad"), "key"]),
    "VirSorter1":  set(df.loc[m.str.contains("virsorter1"), "key"]),
    "VirSorter2":  set(df.loc[m.str.contains("virsorter2"), "key"]),
    "VIBRANT":     set(df.loc[m.str.contains("vibrant"), "key"]),
    "PlasX":       set(df.loc[m.str.contains("plasx"), "key"]),
}
for k, v in sets.items():
    print(f"{k:12s} {len(v):6d}")
print(f"{'UNION':12s} {len(set().union(*sets.values())):6d}  (table rows {len(df)})")

# --- source data: exact count for every non-empty region of the 5-set partition ---
names = list(sets)
rows = []
for r in range(1, len(names) + 1):
    for combo in combinations(names, r):
        inc = set.intersection(*[sets[n] for n in combo])
        exc = set().union(*[sets[n] for n in names if n not in combo]) if r < len(names) else set()
        only = inc - exc
        if only:
            rows.append({"region": " & ".join(combo), "n_callers": r, "count": len(only)})
sd = pd.DataFrame(rows).sort_values(["n_callers", "count"], ascending=[True, False])
sd.to_csv(f"{OUTDIR}/venn5_callers_sourcedata.csv", index=False)
print(f"\nwrote {OUTDIR}/venn5_callers_sourcedata.csv ({len(sd)} non-empty regions)")

# --- figure ---
cmap = ["#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"]  # geNomad, VS1, VS2, VIBRANT, PlasX
fig, ax = plt.subplots(figsize=(11, 10))
venn(sets, cmap=cmap, fontsize=9, legend_loc="upper left", ax=ax)
ax.set_title("ECE detections across callers (unfiltered; 64 metagenomes, >=5x contigs)\n"
             f"{len(df)} distinct ECEs", fontsize=12)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUTDIR}/venn5_callers.{ext}", bbox_inches="tight", dpi=200)
print(f"wrote {OUTDIR}/venn5_callers.pdf / .png")
