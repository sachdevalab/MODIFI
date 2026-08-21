#!/usr/bin/env python3
"""Combine Stage 2b (read-level orientation flip) with Stage 3 (modification
shift) into a per-candidate phase-variation verdict.

Verdict:
  CONFIRMED  = orientation flips across timepoints AND >=1 bipartite motif shifts
  CANDIDATE  = exactly one of {orientation flip, motif shift}
  INCONCLUSIVE = orientation not assessable (too few spanning reads/timepoints)
                 and no motif shift
  NEGATIVE   = neither flip nor shift, orientation assessable
"""

import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("/home/shuaiw/borg/revision/phase_variation")
ORI = BASE / "orientation" / "orientation_ratios.csv"
MOT = BASE / "motif_switch_summary.tsv"
MIN_READS = 5   # per-timepoint (raw+inv) floor to trust an orientation call


def main():
    # orientation per cluster
    ori = defaultdict(list)
    meta = {}
    for r in csv.DictReader(open(ORI)):
        cl = r["cluster"]
        meta[cl] = (r["species"], r["subject"])
        try:
            raw, inv = int(r["raw_crossing"]), int(r["inv_crossing"])
        except ValueError:
            raw = inv = 0
        tot = raw + inv
        frac = (inv / tot) if tot else None
        ori[cl].append((r["timepoint"], raw, inv, tot, frac))

    # motif shifts per cluster
    mot = defaultdict(list)
    for r in csv.DictReader(open(MOT), delimiter="\t"):
        if r["flag"] in ("strong_switch", "soft_shift"):
            mot[r["cluster"]].append((r["motif"], r["flag"], r["fractions"], r["timepoints"]))

    out = BASE / "phase_variation_events.tsv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["cluster", "species", "subject", "orientation_by_tp",
                    "n_tp_assessable", "orientation_flip", "shifting_motifs",
                    "verdict"])
        for cl in sorted(ori, key=lambda c: (meta[c][0])):
            rows = sorted(ori[cl], key=lambda x: float(x[0]))
            assessable = [x for x in rows if x[3] >= MIN_READS]
            fracs = [x[4] for x in assessable]
            if len(assessable) >= 2 and any(fr is not None for fr in fracs):
                flip = "yes" if (min(fracs) < 0.5 and max(fracs) > 0.5) else "no"
            else:
                flip = "inconclusive"
            shifts = mot.get(cl, [])
            shift_str = ";".join(f"{m}({flag})" for m, flag, _, _ in shifts) or "none"

            if flip == "yes" and shifts:
                verdict = "CONFIRMED"
            elif flip == "yes" or shifts:
                verdict = "CANDIDATE"
            elif flip == "inconclusive":
                verdict = "INCONCLUSIVE"
            else:
                verdict = "NEGATIVE"

            ori_str = ",".join(
                f"tp{tp}:{'NA' if fr is None else format(fr, '.2f')}(n={tot})"
                for tp, raw, inv, tot, fr in rows)
            w.writerow([cl, meta[cl][0], meta[cl][1], ori_str,
                        len(assessable), flip, shift_str, verdict])
    # echo
    print(open(out).read())


if __name__ == "__main__":
    main()
