#!/usr/bin/env python3
"""Aggregate the new-ECE host-linkage results and merge with the existing geNomad linkages.

For each sample with filter-passing new ECEs, read the fresh linkage output
ece_linkage/<sample>/meth/host_summary.csv. A new ECE is "linked" at final_score>0.5 &
specificity<0.01 (the host_summary.filter.csv rule); we also report the full score
distribution and a stricter >0.6 variant. Existing geNomad linkages are read (read-only)
from run2/<sample>/<mm>/host_summary.filter.csv for the merged network total.

Outputs: tables under ece_linkage/, figures + source data under tmp/rev_figs/ece_anno/."""
import os, glob
import pandas as pd
import numpy as np

L = "/home/shuaiw/borg/revision/ece_linkage"
RUN2 = "/home/shuaiw/borg/paper/run2"
D = "/home/shuaiw/borg/revision/ece_anno"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
os.makedirs(FIG, exist_ok=True)
SCORE_MIN, SPEC_MAX = 0.5, 0.01

# metadata for the new ECEs: caller methods + environment + type
comb = pd.read_csv(f"{D}/expanded/combined_metagenome_eces.csv")
comb["key"] = comb["sample"] + "|" + comb["MGE"].astype(str)
meta = comb.set_index("key")[["MGE_type", "environment", "methods"]].to_dict("index")


def read_summary(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    try:
        d = pd.read_csv(path)
    except Exception:
        return pd.DataFrame()
    return d if "MGE" in d.columns else pd.DataFrame()


# ---- new-ECE linkages: best host per new ECE ----
rows = []
samples = [s.strip() for s in open(f"{L}/link_samples_frozen.txt") if s.strip()]
extra = f"{L}/link_samples_extra.txt"
if os.path.exists(extra):
    samples += [s.strip() for s in open(extra) if s.strip()]
for s in sorted(set(samples)):
    d = read_summary(f"{L}/{s}/meth/host_summary.csv")
    if d.empty:
        continue
    d["final_score"] = pd.to_numeric(d["final_score"], errors="coerce")
    d["specificity"] = pd.to_numeric(d["specificity"], errors="coerce")
    best = d.sort_values("final_score", ascending=False).drop_duplicates("MGE")
    best["sample"] = s
    rows.append(best)
new = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
new["key"] = new["sample"] + "|" + new["MGE"].astype(str)
new["MGE_type"] = new["key"].map(lambda k: meta.get(k, {}).get("MGE_type", "NA"))
new["environment"] = new["key"].map(lambda k: meta.get(k, {}).get("environment", "NA"))
new["methods"] = new["key"].map(lambda k: meta.get(k, {}).get("methods", ""))
new["linked"] = (new["final_score"] > SCORE_MIN) & (new["specificity"] < SPEC_MAX)
new["linked_strict"] = (new["final_score"] > 0.6) & (new["specificity"] < SPEC_MAX)
new.to_csv(f"{L}/new_ece_linkage_all.csv", index=False)

print(f"=== new-ECE linkage (samples scored: {new['sample'].nunique()}) ===")
print(f"new ECEs scored: {len(new)}")
print(f"LINKED (score>{SCORE_MIN} & spec<{SPEC_MAX}): {int(new['linked'].sum())} "
      f"({int(((new.MGE_type=='plasmid')&new.linked).sum())} plasmid, "
      f"{int(((new.MGE_type=='virus')&new.linked).sum())} virus)")
print(f"LINKED strict (score>0.6): {int(new['linked_strict'].sum())}")
print("\nlinked new ECEs by environment:")
print(new[new.linked].groupby("environment").size().sort_values(ascending=False).to_string())
print("\nlinked new ECEs by caller (any of):")
for c in ["virsorter2", "vibrant", "virsorter1", "plasx"]:
    k = int((new.linked & new.methods.str.contains(c)).sum())
    print(f"  {c}: {k}")

# ---- existing geNomad linkages (read-only) for the merged total ----
gen_links = 0
for s in comb["sample"].unique():
    for mm in ("methylation4", "methylation3"):
        p = f"{RUN2}/{s}/{s}_{mm}/host_summary.filter.csv"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                gd = pd.read_csv(p)
            except Exception:
                break
            if "MGE" in gd.columns:
                gen_links += gd["MGE"].nunique()
            break
print(f"\n=== merged network ===")
print(f"existing geNomad-ECE linkages (filter.csv, distinct MGEs): {gen_links}")
print(f"new-ECE linkages added: {int(new['linked'].sum())}")
print(f"TOTAL expanded linkages: {gen_links + int(new['linked'].sum())}")

# ---- figure: new-ECE linkage score distribution ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fs = new["final_score"].dropna()
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].hist(fs, bins=25, range=(0, 1), color="#2c7fb8", edgecolor="white")
ax[0].axvline(SCORE_MIN, color="#d73027", ls="--", lw=1.5, label=f"link cutoff {SCORE_MIN}")
ax[0].set_xlabel("linkage score (final_score)"); ax[0].set_ylabel("new ECEs")
ax[0].set_title(f"New-ECE linkage scores (n={len(fs)})"); ax[0].legend(frameon=False)
lk = new[new.linked]
by = lk.groupby(["environment", "MGE_type"]).size().unstack(fill_value=0)
by.plot(kind="barh", stacked=True, ax=ax[1],
        color={"plasmid": "#1b9e77", "virus": "#7570b3"}, edgecolor="white")
ax[1].set_xlabel("linked new ECEs"); ax[1].set_ylabel(""); ax[1].set_title("Linked new ECEs by environment")
ax[1].legend(frameon=False)
fig.tight_layout()
out = f"{FIG}/new_ece_linkage"
fig.savefig(out + ".pdf", bbox_inches="tight"); fig.savefig(out + ".png", dpi=150, bbox_inches="tight")
new[["sample", "MGE", "MGE_type", "environment", "final_score", "specificity", "linked"]].to_csv(
    out + "_sourcedata.csv", index=False)
print(f"\nwrote {L}/new_ece_linkage_all.csv and {out}.pdf/.png + _sourcedata.csv")
