#!/usr/bin/env python3
"""Stage 3: modification-fraction trajectory of Type I (bipartite) motifs across
time points, for each phase-variation candidate strain.

For every candidate cluster, walk its member contigs at each time point, read the
per-contig motif table, and record the modification fraction of each bipartite
(gapped, Type I-like) motif. A motif that crosses the ON<->OFF boundary across
time points is the epigenetic readout of the invertible specificity locus.

Login-node-safe: reads small per-contig motifs.csv files.

Outputs to /home/shuaiw/borg/revision/phase_variation/:
  motif_trajectory_long.tsv  : one row per (strain, motif, timepoint)
  motif_switch_summary.tsv   : one row per (strain, motif) with min/max/flag
"""

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict

RUN2 = Path("/groups/banfield/projects/multienv/methylation/data/borg/paper/run2")
BASE = Path("/home/shuaiw/borg/revision/phase_variation")
STRAINS = BASE / "longitudinal_strains.tsv"
REGIONS = BASE / "candidate_inversion_regions.tsv"

MIN_COV = 10          # meanCoverage floor for a trustworthy fraction
BIPARTITE_N = 3       # >=3 consecutive N in motifString => Type I-like bipartite
STRONG_HI, STRONG_LO = 0.7, 0.3
SOFT_DELTA = 0.3

# Scope: "all" (default) = every longitudinal strain; "candidates" = only the
# multi-SP Type I clusters in candidate_inversion_regions.tsv.
SCOPE = sys.argv[1] if len(sys.argv) > 1 else "all"


def meth_dir(sample):
    for m in ("methylation4", "methylation3"):
        d = RUN2 / sample / f"{sample}_{m}"
        if d.exists():
            return d
    return None


def is_bipartite(motif_string):
    return bool(re.search("N{%d,}" % BIPARTITE_N, motif_string))


def read_motifs(sample, contig):
    d = meth_dir(sample)
    if d is None:
        return {}
    f = d / "motifs" / f"{contig}.motifs.csv"
    if not f.exists():
        return {}
    out = {}
    for r in csv.DictReader(open(f)):
        ms = r["motifString"]
        mid = f"{ms}_{r['centerPos']}"
        try:
            frac = float(r["fraction"])
            cov = float(r.get("meanCoverage", "0") or 0)
            ndet = int(r.get("nDetected", "0") or 0)
            ngen = int(r.get("nGenome", "0") or 0)
        except ValueError:
            continue
        out[mid] = dict(fraction=frac, coverage=cov, nDetected=ndet, nGenome=ngen,
                        motif_type="bipartite" if is_bipartite(ms) else "simple")
    return out


def main():
    strains = list(csv.DictReader(open(STRAINS), delimiter="\t"))
    cand_clusters = sorted({r["cluster"] for r in csv.DictReader(open(REGIONS), delimiter="\t")})

    long_rows = []
    # traj[(cluster, subject, motif)] = {timepoint: dict}
    traj = defaultdict(dict)
    meta = {}
    mtype = {}
    for st in strains:
        if SCOPE == "candidates" and st["cluster"] not in cand_clusters:
            continue
        for member in st["members"].split(";"):
            tp, contig = member.split(":", 1)
            mm = re.match(r"^(?P<sample>.+)_\d+_[CL]$", contig)
            if not mm:
                continue
            sample = mm.group("sample")
            motifs = read_motifs(sample, contig)
            for mid, d in motifs.items():
                long_rows.append(dict(cluster=st["cluster"], subject=st["subject"],
                                      species=st["species"], timepoint=tp,
                                      sample=sample, contig=contig, motif=mid,
                                      **d))
                # keep the highest-coverage observation per timepoint for the trajectory
                k = (st["cluster"], st["subject"], mid)
                meta[k] = st["species"]
                mtype[k] = d["motif_type"]
                prev = traj[k].get(tp)
                if prev is None or d["coverage"] > prev["coverage"]:
                    traj[k][tp] = d

    # write long
    lp = BASE / "motif_trajectory_long.tsv"
    with open(lp, "w", newline="") as f:
        cols = ["cluster", "subject", "species", "timepoint", "sample", "contig",
                "motif", "motif_type", "fraction", "coverage", "nDetected", "nGenome"]
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in long_rows:
            w.writerow(r)

    # summary per (strain, motif) using only adequately covered timepoints
    sp = BASE / "motif_switch_summary.tsv"
    n_strong = 0
    with open(sp, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["cluster", "subject", "species", "motif", "motif_type",
                    "n_tp_covered", "timepoints", "fractions", "min_frac",
                    "max_frac", "delta", "flag"])
        n_soft = 0
        for k, tpmap in sorted(traj.items()):
            cluster, subject, mid = k
            good = {tp: d for tp, d in tpmap.items() if d["coverage"] >= MIN_COV}
            if len(good) < 2:
                continue
            tps = sorted(good, key=lambda x: float(x))
            fr = [good[tp]["fraction"] for tp in tps]
            mn, mx = min(fr), max(fr)
            delta = mx - mn
            if mx >= STRONG_HI and mn <= STRONG_LO:
                flag = "strong_switch"
                n_strong += 1
            elif delta >= SOFT_DELTA:
                flag = "soft_shift"
                n_soft += 1
            else:
                flag = "stable"
            w.writerow([cluster, subject, meta[k], mid, mtype[k], len(good),
                        ",".join(tps), ",".join(f"{x:.3f}" for x in fr),
                        f"{mn:.3f}", f"{mx:.3f}", f"{delta:.3f}", flag])

    print(f"Wrote {lp} ({len(long_rows)} rows)  scope={SCOPE}")
    print(f"Wrote {sp}")
    print(f"  strong ON<->OFF switches: {n_strong};  soft shifts: {n_soft}")


if __name__ == "__main__":
    main()
