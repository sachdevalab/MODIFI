#!/usr/bin/env python3
"""agg_hosts.py — build a host_summary from per-ECE hosts/*.host_prediction.csv.

MODIFI's bin-mode linkage writes correct per-ECE ranked host lists (bin-level), but its
final host_summary-write step crashes in this fragmented/bin setup. This bypasses it:
take the top-scoring host (bin) per ECE. Output has the columns score_mag.py needs.

Usage: agg_hosts.py <hosts_dir> <out.host_summary.csv>
"""
import sys, os, glob
import pandas as pd


def main():
    hosts_dir, out = sys.argv[1], sys.argv[2]
    rows = []
    for f in glob.glob(os.path.join(hosts_dir, "*.host_prediction.csv")):
        d = pd.read_csv(f)
        if len(d):
            rows.append(d.sort_values("final_score", ascending=False).iloc[0])
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"[agg_hosts] {len(df)} ECEs -> {out}")


if __name__ == "__main__":
    main()
