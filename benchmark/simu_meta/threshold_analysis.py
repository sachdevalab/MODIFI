#!/usr/bin/env python3
"""
threshold_analysis.py — data-driven calibration of MODIFI's linkage thresholds
(Part G / C5), prototyped on the completed rep-1 C1 communities.

Answers "how should linkage_score / specificity / motif-fraction / min-sites be set?"
using ground truth (an ECE-host pair is correct iff same SRA prefix). Instead of
asserting cutoffs, we show where correct and incorrect links separate and whether the
chosen values sit in a robust (flat) region.

Two threshold groups:
  - LINKAGE-DECISION (final_score, specificity): calibrated on the full candidate-host
    ranking (hosts/*.host_prediction.csv) — every ECE x host pair is a labeled example.
  - MOTIF-EVIDENCE (motif fraction, #modified sites): from motif_info of the motifs that
    actually drive correct links; show they sit far above 0.3 / 100 (robustness).
"""
import os
import re
import glob
import numpy as np
import pandas as pd

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
COMMUNITIES = ["ladder_25", "ladder_40", "ladder_58", "bg_80", "bg_150"]


def sra(name):
    return str(name).split("_")[0]


def load_pairs():
    """All (ECE, candidate host) pairs across communities, labeled correct/incorrect."""
    rows = []
    for c in COMMUNITIES:
        for f in glob.glob(f"{ROOT}/{c}/modifi/{c}/hosts/*.host_prediction.csv"):
            try:
                d = pd.read_csv(f)
            except Exception:
                continue
            if d.empty:
                continue
            d["community"] = c
            d["correct"] = d["MGE"].map(sra) == d["host"].map(sra)
            d["rank"] = range(1, len(d) + 1)
            rows.append(d)
    return pd.concat(rows, ignore_index=True)


def best_host(pairs):
    """Top-ranked candidate host per ECE (the actual decision unit)."""
    return pairs.sort_values("rank").groupby(["community", "MGE"], as_index=False).first()


def sweep_score(best, spec_cut=1.0):
    """Recall/precision/F1/FDR as final_score cutoff varies (best-host accepted if
    final_score>s AND specificity<spec_cut). Denominator = all ECEs (unaccepted=FN)."""
    n = len(best)
    out = []
    for s in np.round(np.arange(0.0, 1.001, 0.05), 3):
        acc = best[(best.final_score > s) & (best.specificity < spec_cut)]
        tp = int(acc.correct.sum())
        fp = len(acc) - tp
        recall = tp / n
        prec = tp / (tp + fp) if (tp + fp) else np.nan
        f1 = 2 * prec * recall / (prec + recall) if prec and recall else 0.0
        fdr = fp / (tp + fp) if (tp + fp) else 0.0
        out.append(dict(score_cut=s, tp=tp, fp=fp, recall=round(recall, 3),
                        precision=round(prec, 3) if prec == prec else np.nan,
                        f1=round(f1, 3), fdr=round(fdr, 3)))
    return pd.DataFrame(out)


def parse_motif_info(s):
    """motif_info entries: motif_cp:cp:host_total:host_meth:pl_total:pl_meth:occ:occ_len"""
    res = []
    if not isinstance(s, str):
        return res
    for e in s.split(";"):
        p = e.split(":")
        if len(p) >= 6:
            try:
                ht, hm, pt, pm = int(p[2]), int(p[3]), int(p[4]), int(p[5])
                res.append(dict(host_total=ht, host_meth=hm, pl_total=pt, pl_meth=pm,
                                host_frac=hm / ht if ht else 0,
                                pl_frac=pm / pt if pt else 0))
            except ValueError:
                pass
    return res


def main():
    pairs = load_pairs()
    best = best_host(pairs)
    print(f"[data] {len(pairs)} ECE x host pairs; {best.MGE.nunique()} ECEs; "
          f"{len(COMMUNITIES)} communities (rep 1)")

    # ---- separation of correct vs incorrect best-host calls ----
    corr = best[best.correct]
    inc = best[~best.correct]
    def q(x, col): return np.round(np.percentile(x[col], [10, 50, 90]), 4).tolist() if len(x) else []
    print("\n=== LINKAGE-DECISION thresholds: correct vs incorrect BEST-host calls ===")
    print(f"  correct best-hosts:   n={len(corr)}")
    print(f"    final_score [p10,med,p90] = {q(corr,'final_score')}")
    print(f"    specificity [p10,med,p90] = {q(corr,'specificity')}")
    print(f"  incorrect best-hosts: n={len(inc)}")
    print(f"    final_score [p10,med,p90] = {q(inc,'final_score')}")
    print(f"    specificity [p10,med,p90] = {q(inc,'specificity')}")

    # ---- final_score sweep (at specificity<0.01, the code default) ----
    sw = sweep_score(best, spec_cut=0.01)
    print("\n=== final_score sweep (specificity<0.01) ===")
    print(sw.to_string(index=False))
    # recommended: max-F1
    best_f1 = sw.loc[sw.f1.idxmax()]
    print(f"\n  max-F1 at final_score>{best_f1.score_cut}: "
          f"recall={best_f1.recall}, precision={best_f1.precision}, F1={best_f1.f1}")
    for cut in (0.5, 0.8):
        r = sw[sw.score_cut == cut]
        if len(r):
            r = r.iloc[0]
            print(f"  at final_score>{cut} (default/reviewer): recall={r.recall}, "
                  f"precision={r.precision}, FDR={r.fdr}")

    # ---- specificity sweep (at final_score>0.5) ----
    print("\n=== specificity sweep (final_score>0.5) ===")
    n = len(best)
    for p in [1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 1.0]:
        acc = best[(best.final_score > 0.5) & (best.specificity < p)]
        tp = int(acc.correct.sum()); fp = len(acc) - tp
        print(f"  specificity<{p:>7}: recall={tp/n:.3f}  precision="
              f"{tp/(tp+fp) if tp+fp else float('nan'):.3f}  (tp={tp}, fp={fp})")

    # ---- MOTIF-EVIDENCE: fraction & #modified sites of link-driving motifs ----
    print("\n=== MOTIF-EVIDENCE thresholds (motifs driving CORRECT links) ===")
    host_fracs, host_meths, pl_fracs = [], [], []
    for _, r in corr.iterrows():
        for m in parse_motif_info(r.motif_info):
            # only motifs the ECE is actually methylated at (the ones that count)
            if m["pl_frac"] >= 0.4:
                host_fracs.append(m["host_frac"]); host_meths.append(m["host_meth"])
                pl_fracs.append(m["pl_frac"])
    host_fracs, host_meths = np.array(host_fracs), np.array(host_meths)
    if len(host_fracs):
        print(f"  n matched motifs in correct links: {len(host_fracs)}")
        print(f"  host methylation FRACTION [p10,med,p90] = "
              f"{np.round(np.percentile(host_fracs,[10,50,90]),3).tolist()}  (threshold cited: >0.3)")
        print(f"    fraction below 0.3: {100*np.mean(host_fracs<0.3):.1f}% of link-driving motifs")
        print(f"  host #MODIFIED SITES [p10,med,p90] = "
              f"{np.round(np.percentile(host_meths,[10,50,90])).tolist()}  (threshold cited: >=100)")
        print(f"    sites below 100: {100*np.mean(host_meths<100):.1f}% of link-driving motifs")

    outdir = f"{ROOT}/eval"
    os.makedirs(outdir, exist_ok=True)
    sw.to_csv(f"{outdir}/threshold_score_sweep.csv", index=False)
    best.to_csv(f"{outdir}/besthost_calls.csv", index=False)
    print(f"\n[out] wrote {outdir}/threshold_score_sweep.csv, besthost_calls.csv")


if __name__ == "__main__":
    main()
