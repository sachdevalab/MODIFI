#!/usr/bin/env python3
"""5-set Venn of FILTER-PASSING ECEs across geNomad, VirSorter1/2, VIBRANT, PlasX.

Filter (loose metagenome criteria, "any strong caller" P2 auto-satisfied by combined-set
membership): pass = support_marker (P3) AND NOT flag_chromosomal (rRNA present) AND NOT
flag_artifact. Evidence is gathered from all available dirs (metagenome_all = geNomad;
expanded_new/delta/refresh = non-geNomad). ECEs without evidence yet (e.g. the 6 big soil
samples still running) are excluded -> this is the CURRENT filtered picture.

Outputs (next to the figure):
  tmp/rev_figs/ece_anno/venn5_callers_filtered.pdf / .png
  tmp/rev_figs/ece_anno/venn5_callers_filtered_sourcedata.csv  (per-region counts)
"""
import os, glob
from itertools import combinations
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from venn import venn

D = "/home/shuaiw/borg/revision/ece_anno"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
os.makedirs(OUTDIR, exist_ok=True)

comb = pd.read_csv(f"{D}/expanded/combined_metagenome_eces.csv")
comb["key"] = comb["sample"] + "|" + comb["MGE"].astype(str)

# --- gather evidence + loose pass rule ---
frames = []
for base in ("metagenome_all", "expanded_new", "expanded_delta", "expanded_refresh"):
    for tsv in glob.glob(f"{D}/{base}/per_sample/*/ece_evidence.tsv"):
        s = os.path.basename(os.path.dirname(tsv))
        try:
            d = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        if d.empty or "seq_name" not in d.columns:
            continue
        d["sample"] = s
        frames.append(d[["sample", "seq_name", "support_marker", "flag_artifact", "rrna_count"]])
ev = pd.concat(frames, ignore_index=True)
ev["key"] = ev["sample"] + "|" + ev["seq_name"].astype(str)
ev = ev.drop_duplicates("key")
p3 = ev["support_marker"].fillna(False).astype(bool)
flag_chr = pd.to_numeric(ev["rrna_count"], errors="coerce").fillna(0) > 0
flag_art = ev["flag_artifact"].fillna(False).astype(bool)
ev["pass"] = p3 & ~flag_chr & ~flag_art
pass_keys = set(ev.loc[ev["pass"], "key"])

# --- keep passing ECEs present in the combined set; report coverage ---
comb_ev = comb[comb["key"].isin(set(ev["key"]))]                 # have evidence
passing = comb[comb["key"].isin(pass_keys)].copy()
print(f"combined ECEs: {len(comb)}; with evidence: {len(comb_ev)}; "
      f"MISSING evidence (excluded): {len(comb) - len(comb_ev)}")
print(f"FILTER-PASSING (current): {len(passing)} "
      f"(plasmid {(passing.MGE_type=='plasmid').sum()}, virus {(passing.MGE_type=='virus').sum()})")

m = passing["methods"].fillna("")
sets = {
    "geNomad":    set(passing.loc[m.str.contains("genomad"), "key"]),
    "VirSorter1": set(passing.loc[m.str.contains("virsorter1"), "key"]),
    "VirSorter2": set(passing.loc[m.str.contains("virsorter2"), "key"]),
    "VIBRANT":    set(passing.loc[m.str.contains("vibrant"), "key"]),
    "PlasX":      set(passing.loc[m.str.contains("plasx"), "key"]),
}
print("\nfilter-passing ECEs per caller:")
for k, v in sets.items():
    print(f"  {k:12s} {len(v)}")

# --- source data: per-region counts ---
names = list(sets); rows = []
for r in range(1, len(names) + 1):
    for combo in combinations(names, r):
        inc = set.intersection(*[sets[n] for n in combo])
        exc = set().union(*[sets[n] for n in names if n not in combo]) if r < len(names) else set()
        only = inc - exc
        if only:
            rows.append({"region": " & ".join(combo), "n_callers": r, "count": len(only)})
sd = pd.DataFrame(rows).sort_values(["n_callers", "count"], ascending=[True, False])
sd.to_csv(f"{OUTDIR}/venn5_callers_filtered_sourcedata.csv", index=False)

# --- figure ---
cmap = ["#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756"]
fig, ax = plt.subplots(figsize=(11, 10))
venn(sets, cmap=cmap, fontsize=9, legend_loc="upper left", ax=ax)
ax.set_title("Filter-passing ECEs across callers (64 metagenomes)\n"
             f"{len(passing)} high-confidence ECEs "
             f"({(passing.MGE_type=='plasmid').sum()} plasmid, {(passing.MGE_type=='virus').sum()} virus)",
             fontsize=12)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUTDIR}/venn5_callers_filtered.{ext}", bbox_inches="tight", dpi=200)
print(f"\nwrote {OUTDIR}/venn5_callers_filtered.pdf / .png / _sourcedata.csv")
