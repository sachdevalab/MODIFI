#!/usr/bin/env python3
"""5-set Venn of the REVISED filter-passing ECE set (3,976) across geNomad, VirSorter1/2,
VIBRANT, PlasX. Reads filterpass_revised_final.csv (methods column). Outputs:
  tmp/rev_figs/ece_anno/venn5_callers_filtered.pdf / .png / _sourcedata.csv
"""
import os
from itertools import combinations
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from venn import venn

D = "/home/shuaiw/borg/revision/ece_anno"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
os.makedirs(OUTDIR, exist_ok=True)

p = pd.read_csv(f"{D}/expanded/filterpass_FINAL.csv")
p["key"] = p["sample"] + "|" + p["MGE"].astype(str)
m = p["methods"].fillna("")
sets = {
    "geNomad":    set(p.loc[m.str.contains("genomad"), "key"]),
    "VirSorter1": set(p.loc[m.str.contains("virsorter1"), "key"]),
    "VirSorter2": set(p.loc[m.str.contains("virsorter2"), "key"]),
    "VIBRANT":    set(p.loc[m.str.contains("vibrant"), "key"]),
    "PlasX":      set(p.loc[m.str.contains("plasx"), "key"]),
}
np_, nv = int((p.MGE_type == "plasmid").sum()), int((p.MGE_type == "virus").sum())
print(f"revised filter-passing: {len(p)} ({np_} plasmid, {nv} virus)")
for k, v in sets.items():
    print(f"  {k:12s} {len(v)}")

names = list(sets); rows = []
for r in range(1, len(names) + 1):
    for combo in combinations(names, r):
        inc = set.intersection(*[sets[n] for n in combo])
        exc = set().union(*[sets[n] for n in names if n not in combo]) if r < len(names) else set()
        only = inc - exc
        if only:
            rows.append({"region": " & ".join(combo), "n_callers": r, "count": len(only)})
pd.DataFrame(rows).sort_values(["n_callers", "count"], ascending=[True, False]).to_csv(
    f"{OUTDIR}/venn5_callers_filtered_sourcedata.csv", index=False)

cmap = ["#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"]
fig, ax = plt.subplots(figsize=(11, 10))
venn(sets, cmap=cmap, fontsize=9, legend_loc="upper left", ax=ax)
ax.set_title(f"Filter-passing ECEs across callers (64 metagenomes, revised criteria)\n"
             f"{len(p)} ECEs ({np_} plasmid, {nv} virus)", fontsize=12)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUTDIR}/venn5_callers_filtered.{ext}", bbox_inches="tight", dpi=200)
print(f"wrote {OUTDIR}/venn5_callers_filtered.pdf / .png / _sourcedata.csv")
