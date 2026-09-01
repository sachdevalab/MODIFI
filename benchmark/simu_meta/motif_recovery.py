#!/usr/bin/env python3
"""motif_recovery.py -- motif-detection recovery vs con-specific strain depth K.

Ground truth = each community genome's own solo MODIFI run motif set (all donor + background
genomes). Recovered = the motifs detected on that genome's contigs in the strain-mixture community.
The SAME filter (fraction >= MIN_FRAC and nDetected >= MIN_SITES) is applied to both sides; motif
strings are reverse-complement canonicalized. Recovery is pooled across genomes per rep, summarized
mean +/- 95% CI vs K. Pure re-analysis of files already on disk (no MODIFI runs).

Writes tmp/rev_figs/simu_meta/strain_het/strain_mix_motif_recovery_sourcedata.csv.
"""
import os
import glob
import numpy as np
import pandas as pd
from Bio.Seq import Seq

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
ISO_ROOTS = ["/home/shuaiw/borg/paper/isolation/batch2_results",
             "/home/shuaiw/borg/paper/isolation/bacteria"]
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/strain_het"
MIN_FRAC, MIN_SITES = 0.4, 100
KDEF = [("1", "strain_mix_k1"), ("2", "strain_mix_k2"), ("3", "strain_mix_k3"),
        ("4", "strain_mix_k4"), ("all", "strain_mix_kall")]
os.makedirs(OUT, exist_ok=True)


def canon(m):
    """RC-canonical motif string (IUPAC-aware), matching sample_object.get_unique_motifs."""
    try:
        rc = str(Seq(m).reverse_complement())
    except Exception:
        rc = m
    return min(m, rc)


def filt_set(df):
    """apply the shared filter to a genome-level motif table -> set of RC-canonical motif strings."""
    d = df[(df["fraction"] >= MIN_FRAC) & (df["nDetected"] >= MIN_SITES)]
    return {canon(m) for m in d["motifString"]}


def ground_truth(acc, cache):
    if acc in cache:
        return cache[acc]
    s = None
    for r in ISO_ROOTS:
        f = os.path.join(r, acc, f"{acc}_methylation4", "all.motifs.csv")
        if os.path.exists(f):
            try:
                s = filt_set(pd.read_csv(f))
            except Exception:
                s = set()
            break
    cache[acc] = s          # None => no ground-truth file found on disk
    return s


def rep_labels(base):
    return [l for l in [f"{base}_rep{r}" for r in range(1, 6)]
            if os.path.exists(f"{ROOT}/{l}/modifi/{l}/all.motifs.merged.csv")]


def rep_recovery(lab, gt_cache, missing):
    man = pd.read_csv(f"{ROOT}/{lab}/{lab}.manifest.csv")
    members = set(man["sample"].astype(str))            # ALL community genomes (donor + background)
    m = pd.read_csv(f"{ROOT}/{lab}/modifi/{lab}/all.motifs.merged.csv")
    m["iso"] = m["contig"].astype(str).str.split("_").str[0]
    m = m[m.iso.isin(members)]
    # aggregate community motifs to genome level per (isolate, motif)
    agg = (m.groupby(["iso", "motifString"], as_index=False)
             .agg(nDetected=("nDetected", "sum"), nGenome=("nGenome", "sum")))
    agg["fraction"] = agg.nDetected / agg.nGenome.replace(0, np.nan)
    inter = tot = n_don = 0
    per = []
    for acc, g in agg.groupby("iso"):
        G = ground_truth(acc, gt_cache)
        if G is None:
            missing.add(acc); continue
        if not G:
            continue
        D = filt_set(g)
        i = len(D & G)
        inter += i; tot += len(G); n_don += 1
        per.append(i / len(G))
    return (inter / tot if tot else np.nan), (np.mean(per) if per else np.nan), n_don


def agg(v):
    v = np.array([x for x in v if x == x], float)
    if not len(v):
        return np.nan, 0.0
    return v.mean(), (1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def main():
    gt_cache, missing = {}, set()
    rows = []
    for k, base in KDEF:
        labs = rep_labels(base)
        pooled, perdonor, ndon = [], [], []
        for lab in labs:
            p, pd_, nd = rep_recovery(lab, gt_cache, missing)
            pooled.append(p); perdonor.append(pd_); ndon.append(nd)
        pm, pc = agg(pooled)
        rows.append(dict(K=k, n_rep=len(labs), n_genome=int(np.mean(ndon)) if ndon else 0,
                         recovery_pooled_mean=pm, recovery_pooled_ci=pc,
                         recovery_perdonor_mean=agg(perdonor)[0]))
        print(f"  K={k:>3}: n_rep={len(labs)} n_genome~{rows[-1]['n_genome']} "
              f"pooled_recovery={pm:.3f} +/- {pc:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/strain_mix_motif_recovery_sourcedata.csv", index=False)
    print(f"wrote {OUT}/strain_mix_motif_recovery_sourcedata.csv")
    if missing:
        print(f"[warn] {len(missing)} donor accessions had no ground-truth all.motifs.csv: "
              f"{sorted(missing)[:8]}{' ...' if len(missing) > 8 else ''}")


if __name__ == "__main__":
    main()
