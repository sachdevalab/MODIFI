#!/usr/bin/env python3
"""Build per-sample MGE list files for the FILTER-PASSING new (non-geNomad) ECEs, to feed
MODIFI's `--run_steps host` linkage step (custom --mge_file).

A new ECE passes when: it is non-geNomad (so P2 "any strong caller" is auto-satisfied) AND
support_marker (P3) AND not flag_chromosomal AND not flag_artifact, using the per-sample
evidence in expanded_new/per_sample/<sample>/ece_evidence.tsv.

Writes /home/shuaiw/borg/revision/ece_linkage/<sample>/new_eces.mge.tsv (columns: seq_name, length).
Read-only w.r.t. all inputs. Safe to re-run (overwrites only its own output under ece_linkage/)."""
import os, glob
import pandas as pd

D = "/home/shuaiw/borg/revision/ece_anno"
OUT = "/home/shuaiw/borg/revision/ece_linkage"
os.makedirs(OUT, exist_ok=True)

comb = pd.read_csv(f"{D}/expanded/combined_metagenome_eces.csv")
comb = comb[~comb.methods.str.contains("genomad")].copy()   # NEW non-geNomad only
comb["key"] = comb["sample"] + "|" + comb["MGE"].astype(str)
len_map = dict(zip(comb["key"], comb["mge_len"]))

# gather evidence (both the initial expanded_new run and the VS1-completion delta run)
frames = []
for base in ("expanded_new", "expanded_delta"):
    for tsv in glob.glob(f"{D}/{base}/per_sample/*/ece_evidence.tsv"):
        s = os.path.basename(os.path.dirname(tsv))
        try:
            d = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        if d.empty:
            continue
        d["sample"] = s
        frames.append(d)
ev = pd.concat(frames, ignore_index=True)
ev["key"] = ev["sample"] + "|" + ev["seq_name"].astype(str)
ev = ev.drop_duplicates("key")

# restrict evidence to the new (non-geNomad) set and apply the pass rule
ev = ev[ev["key"].isin(set(comb["key"]))].copy()
p3 = ev["support_marker"].fillna(False).astype(bool)
neg = ev["flag_chromosomal"].fillna(False).astype(bool) | ev["flag_artifact"].fillna(False).astype(bool)
ev["pass"] = p3 & ~neg
passing = ev[ev["pass"]].copy()
passing["length"] = passing["key"].map(len_map)

n_written = 0
for s, g in passing.groupby("sample"):
    sd = os.path.join(OUT, s)
    os.makedirs(sd, exist_ok=True)
    tsv = os.path.join(sd, "new_eces.mge.tsv")
    g[["seq_name", "length"]].to_csv(tsv, sep="\t", index=False)
    n_written += 1

print(f"evidence samples available: {ev['sample'].nunique()}")
print(f"filter-passing new ECEs: {len(passing)} "
      f"({(passing['type']=='plasmid').sum()} plasmid, {(passing['type']=='virus').sum()} virus)")
print(f"wrote {n_written} per-sample mge files under {OUT}/<sample>/new_eces.mge.tsv")
# a manifest of samples with new ECEs (for the SLURM array)
with open(os.path.join(OUT, "samples_with_new_eces.txt"), "w") as f:
    f.write("\n".join(sorted(passing["sample"].unique())) + "\n")
print(f"wrote {OUT}/samples_with_new_eces.txt")
