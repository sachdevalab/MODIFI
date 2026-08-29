#!/usr/bin/env python3
"""Standalone figure data/stats: paired host-vs-ECE modification-motif recognition-site density in
the strict isolate ECE set (known ECE-host linkage). For each ECE, x = its host chromosome's
host-motif density, y = the ECE's host-motif density. Writes host_motif_density_sourcedata.csv
(read by plot_host_motif_density.R); no figure rendered."""
import os
import numpy as np
import pandas as pd
from scipy import stats
import ece_plot_common as C

STEM = "fig_host_motif_density"


def main():
    df = C.load_master()
    ev = df[(df.dataset == "isolate") & df.passes_full_gate
            & df.motif_density_per_kb.notna() & df.host_own_density_per_kb.notna()].copy()

    out = ev[["type", "host_own_density_per_kb", "motif_density_per_kb"]].rename(
        columns={"host_own_density_per_kb": "host_density", "motif_density_per_kb": "ece_density"})
    out = out[out.type.isin(["plasmid", "virus"])]
    os.makedirs(C.OUT, exist_ok=True)
    out.to_csv(os.path.join(C.OUT, f"{STEM}_sourcedata.csv"), index=False)

    for t in ["plasmid", "virus"]:
        sub = out[out.type == t]
        h, e = sub.host_density.values, sub.ece_density.values
        # paired Wilcoxon signed-rank: is the ECE depleted relative to its own host chromosome?
        p = stats.wilcoxon(h, e).pvalue
        below = float(np.mean(e < h))
        print(f"{t:8s} n={len(sub)}  host median {np.median(h):.2f}  ECE median {np.median(e):.2f}  "
              f"below y=x {below*100:.0f}%  Wilcoxon paired p={p:.2e}")
    print(f"wrote {STEM}_sourcedata.csv")


if __name__ == "__main__":
    main()
