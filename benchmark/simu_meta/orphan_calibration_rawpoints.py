#!/usr/bin/env python3
"""orphan_calibration_rawpoints.py -- detailed (per-replicate) Source Data for the threshold
calibration figure. The 8 aggregated CSVs next to orphan_calibration.pdf store only mean +/- 95%
CI over the 5 orphan_300 replicates; Nature Source Data should expose the raw value behind every
independent sample point. This writes ONE long-format companion CSV with one row per
(panel, parameter value, replicate): the per-rep recall / planted-precision / orphan-FPR at the
same operating point (final_score > 0.5 & specificity < 0.01), reusing the exact metric logic of
plot_orphan_calibration.py.

It then VERIFIES that the per-rep group mean/CI reproduce each existing aggregated CSV (correctness
gate for the tag->host_summary mapping of the 4 call-parameter panels, whose original generators
were ad-hoc). Pure re-analysis of finished runs; no MODIFI re-runs, tool untouched.
"""
import os
import numpy as np
import pandas as pd

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold"
S0, P0 = 0.5, 0.01
REPS = [f"orphan_300_rep{r}" for r in range(1, 6)]
sra = lambda x: str(x).split("_")[0]


def load_hs(lab, path):
    """load a host_summary.csv and tag orphan/correct using this rep's orphan list (== plot script)."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    orph = set(l.strip() for l in open(f"{ROOT}/{lab}/{lab}.orphans.txt") if l.strip())
    try:
        h = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return None
    if len(h) == 0:
        return None
    h["is_orphan"] = h.MGE.isin(orph)
    h["correct"] = (~h.is_orphan) & (h.MGE.map(sra) == h.host.map(sra))
    return h, int((~h.is_orphan).sum()), int(h.is_orphan.sum())


def metrics(hnn, s=S0, p=P0):
    """(recall, planted precision, orphan-FPR) at operating point (s, p); == plot script."""
    if hnn is None:
        return (np.nan, np.nan, np.nan)
    h, n_pl, n_or = hnn
    acc = h[(h.final_score > s) & (h.specificity < p)]
    tp = int(((~acc.is_orphan) & acc.correct).sum())
    planted_conf = int((~acc.is_orphan).sum())
    orphFP = int(acc.is_orphan.sum())
    prec = tp / planted_conf if planted_conf else np.nan
    return (tp / n_pl if n_pl else np.nan, prec, orphFP / n_or if n_or else np.nan)


def base_path(lab):
    return f"{ROOT}/{lab}/modifi/{lab}/host_summary.csv"


def relink_path(lab, tag):
    return f"{ROOT}/{lab}/relink_sweep/{tag}/host_summary.csv"


def callparam_path(lab, tag):
    return f"{ROOT}/{lab}/callparam_sweep/{tag}/host_summary.csv"


BASE = "__BASE__"
# panel -> (parameter column name, aggregated-CSV stem, [(value, host_summary resolver)...]).
# For each (value) the resolver returns the per-rep host_summary path; a value maps to the base
# run (default operating parameters) where no dedicated sweep tag exists.
SWEEPS = {
    "min_sites":     ("min_sites",     "minsites",
                      {10: ("relink", "ms10"), 30: ("relink", "mf0.4"), 60: ("relink", "ms60"),
                       100: ("relink", "ms100"), 200: ("relink", "ms200")}),
    "min_frac":      ("min_frac",      "minfrac",
                      {v: ("relink", f"mf{v}") for v in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]}),
    "min_ctg_cov":   ("min_ctg_cov",   "minctgcov",
                      {2: ("relink", "mc2"), 5: ("base", None), 10: ("relink", "mc10"),
                       20: ("relink", "mc20"), 40: ("relink", "mc40")}),
    "min_cov":       ("min_cov",       "mincov",
                      {1: ("base", None), 3: ("callparam", "cov3"), 5: ("callparam", "cov5"),
                       10: ("callparam", "cov10")}),
    "min_score":     ("min_score",     "minscore",
                      {10: ("callparam", "score10"), 20: ("callparam", "score20"), 30: ("base", None),
                       50: ("callparam", "score50"), 100: ("callparam", "score100")}),
    "min_ece_sites": ("min_ece_sites", "minecesites",
                      {1: ("relink", "ece1"), 2: ("base", None), 3: ("relink", "ece3"),
                       5: ("relink", "ece5"), 10: ("relink", "ece10"), 20: ("relink", "ece20")}),
}
KIND = {"base": lambda lab, t: base_path(lab), "relink": relink_path, "callparam": callparam_path}


def resolve(lab, kind, tag):
    return KIND[kind](lab, tag)


def ci(vals, n_denom=None):
    """95% CI half-width. n_denom overrides sqrt(N): some aggregated CSVs divide the SEM by
    sqrt(total reps=5) even where a value has <5 completed reps, others by sqrt(available)."""
    v = np.array([x for x in vals if x == x], float)
    if len(v) < 2:
        return 0.0
    return 1.96 * v.std(ddof=1) / np.sqrt(n_denom if n_denom else len(v))


def main():
    base = {lab: load_hs(lab, base_path(lab)) for lab in REPS}
    raw = []

    # continuous panels: sweep the operating point on each rep's base host_summary
    sgrid = np.linspace(0, 1, 81)
    pgrid = np.logspace(-4, 0, 81)
    for lab in REPS:
        for s in sgrid:
            rec, prec, fpr = metrics(base[lab], s, P0)
            raw.append(dict(panel="final_score", parameter="final_score", value=s, replicate=lab,
                            recall=rec, precision=prec, orphan_fpr=fpr))
        for p in pgrid:
            rec, prec, fpr = metrics(base[lab], S0, p)
            raw.append(dict(panel="specificity", parameter="specificity", value=p, replicate=lab,
                            recall=rec, precision=prec, orphan_fpr=fpr))

    # discrete panels: default operating point on each value's re-run host_summary
    cache = {}
    for panel, (pname, _stem, grid) in SWEEPS.items():
        for v, (kind, tag) in grid.items():
            for lab in REPS:
                path = resolve(lab, kind, tag)
                if path not in cache:
                    cache[path] = load_hs(lab, path)
                rec, prec, fpr = metrics(cache[path])
                if rec != rec and prec != prec and fpr != fpr:
                    continue  # no completed run for this (value, rep); skip (e.g. min_ctg_cov=40)
                raw.append(dict(panel=panel, parameter=pname, value=v, replicate=lab,
                                recall=rec, precision=prec, orphan_fpr=fpr))

    df = pd.DataFrame(raw)[["panel", "parameter", "value", "replicate", "recall", "precision", "orphan_fpr"]]
    out = f"{OUT}/orphan_calibration_rawpoints_sourcedata.csv"
    df.to_csv(out, index=False)
    print(f"wrote {out} ({len(df)} raw points across {df.panel.nunique()} panels)")

    # ---- verification: raw group mean/CI must reproduce each aggregated CSV ----
    print("\nverifying raw group stats reproduce the aggregated CSVs:")
    stems = {"final_score": ("final_score", "score"), "specificity": ("specificity", "specificity"),
             "min_sites": ("min_sites", "minsites"), "min_frac": ("min_frac", "minfrac"),
             "min_ctg_cov": ("min_ctg_cov", "minctgcov"), "min_cov": ("min_cov", "mincov"),
             "min_score": ("min_score", "minscore"), "min_ece_sites": ("min_ece_sites", "minecesites")}
    all_ok = True
    for panel, (xcol, stem) in stems.items():
        agg = pd.read_csv(f"{OUT}/orphan_calibration_{stem}_sourcedata.csv")
        g = df[df.panel == panel]
        bad = 0
        for _, a in agg.iterrows():
            sub = g[np.isclose(g.value, a[xcol])]
            for m in ["recall", "precision", "orphan_fpr"]:
                em = np.nanmean(sub[m]) if len(sub) else np.nan
                am, ac = a.get(f"{m}_mean", np.nan), a.get(f"{m}_ci", np.nan)
                if (am != am) and (em != em):
                    continue  # both NaN (no completed run)
                # mean is the hard gate; CI must match under EITHER sqrt(N) convention
                ci_ok = (ac != ac) or any(abs(ci(sub[m], d) - ac) <= 1e-9
                                          for d in (None, 5))
                if (em != em) or abs(em - am) > 1e-9 or not ci_ok:
                    bad += 1
        status = "OK" if bad == 0 else f"{bad} MISMATCH"
        all_ok &= (bad == 0)
        print(f"  {panel:14s} n_pts={len(g):3d}  {status}")
    print("ALL PANELS OK" if all_ok else "!! mismatches present -- check tag mapping")


if __name__ == "__main__":
    main()
