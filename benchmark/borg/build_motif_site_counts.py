#!/usr/bin/env python
"""
Build a per-contig, per-motif site-count table for the Borg/Methanoperedens profile.

The filtered profile (profile_profile_df_filtered.csv) carries only contig, motifString,
fraction. The number of modified sites per motif lives in the per-contig CROSS-PROFILE
files (one value for every motif scanned in every genome, not just the natively-called
motifs) at:
    gg_run3/<sample>/<sample>_methylation4/profiles/<contig>.motifs.profile.csv
with columns incl.: motifString, centerPos, motif_loci_num, motif_modified_num,
motif_modified_ratio, ...
The motif id used in the profile is "<motifString>_<centerPos>". We take the modified-site
count as `motif_modified_num` so that cross-detected motifs (present in a genome but not in
its native call set) still get a site count.

Output: profile4/motif_site_counts.csv  (columns: contig, motif, nDetected, fraction)
where nDetected = motif_modified_num and fraction = motif_modified_ratio.
"""
import os
import glob
import pandas as pd

GG_RUN = "/home/shuaiw/borg/paper/gg_run3"
SEQ_DIR = "/home/shuaiw/borg/paper/borg_data/profile4"
PROFILE = os.path.join(SEQ_DIR, "profile_profile_df_filtered.csv")
OUT = os.path.join(SEQ_DIR, "motif_site_counts.csv")

# Contigs in the analysis (all types; filtering by sample happens in the R script)
prof = pd.read_csv(PROFILE)
contigs = sorted(prof["contig"].unique())

# Map each contig name -> its per-contig cross-profile file (glob once)
paths = glob.glob(os.path.join(GG_RUN, "*", "*_methylation4", "profiles", "*.motifs.profile.csv"))
by_contig = {}
for p in paths:
    name = os.path.basename(p)[: -len(".motifs.profile.csv")]
    by_contig.setdefault(name, p)  # first match wins

rows = []
missing = []
for c in contigs:
    p = by_contig.get(c)
    if p is None:
        missing.append(c)
        continue
    df = pd.read_csv(p)
    if not {"motifString", "centerPos", "motif_modified_num"}.issubset(df.columns):
        missing.append(c)
        continue
    df = df.copy()
    df["motif"] = df["motifString"].astype(str) + "_" + df["centerPos"].astype(str)
    for _, r in df.iterrows():
        rows.append({
            "contig": c,
            "motif": r["motif"],
            "nDetected": r.get("motif_modified_num"),
            "fraction": r.get("motif_modified_ratio"),
        })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(f"Wrote {OUT}: {len(out)} rows, {out['contig'].nunique()} contigs, "
      f"{out['motif'].nunique()} motifs.")
if missing:
    print(f"[warn] no motif file found for {len(missing)} contigs (e.g. {missing[:3]})")
