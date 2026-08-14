#!/usr/bin/env python3
"""Break down the metagenome ECE-host CRISPR-spacer validation rate.

The overall figure (17/379 = 4.49%) pools two very different populations, so
here we report it split by MGE type (plasmid vs virus) and by environment.

A linkage (MGE, host contig) is counted as CRISPR-validated when the linked host
contig carries a CRISPR spacer that targets that MGE. This uses the same spacer
hits (My_sample.read_spacer, 0 mismatch) as compare_spacer.py; the linkage set is
the 379-linkage table in mge_host_gc_cov.csv.
"""
import os
import sys
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "isolation"))
from sample_object import My_sample

ALL_DIR = "/home/shuaiw/borg/paper/run2/"
CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"


def validate_linkages(csv=CSV, all_dir=ALL_DIR, mismatch_allowed=0):
    """Return a per-linkage DataFrame with a `validated` (0/1) flag."""
    df = pd.read_csv(csv)
    rows = []
    for prefix, grp in df.groupby("sample"):
        obj = My_sample(prefix, all_dir)
        # MGE -> set(host contigs carrying a spacer that hits this MGE)
        spacer = obj.read_spacer(mismatch_allowed=mismatch_allowed)
        for _, r in grp.iterrows():
            mge, host = r["MGE"], r["host"]
            validated = int(mge in spacer and host in spacer[mge])
            rows.append((prefix, r["environment"], r["MGE_type"], mge, host, validated))
    return pd.DataFrame(
        rows, columns=["sample", "environment", "type", "MGE", "host", "validated"]
    )


def _fmt(label, sub):
    n = len(sub)
    k = int(sub["validated"].sum())
    rate = f"{100 * k / n:5.2f}%" if n else "  n/a"
    return f"{label:<28} {k:>3}/{n:<4} = {rate}"


def summarize(v):
    lines = ["=" * 60, "OVERALL", _fmt("all linkages", v), "", "BY MGE TYPE"]
    for t in ["plasmid", "virus"]:
        lines.append(_fmt(t, v[v.type == t]))
    lines.append("")
    lines.append("BY ENVIRONMENT")
    for env in sorted(v.environment.unique()):
        lines.append(_fmt(env, v[v.environment == env]))
    lines.append("")
    lines.append("BY ENVIRONMENT x TYPE")
    for env in sorted(v.environment.unique()):
        for t in ["plasmid", "virus"]:
            sub = v[(v.environment == env) & (v.type == t)]
            if len(sub):
                lines.append(_fmt(f"  {env} / {t}", sub))
    lines.append("=" * 60)
    return "\n".join(lines)


def main():
    v = validate_linkages()
    print(summarize(v))
    out = os.path.join(os.path.dirname(__file__), "crispr_validation_breakdown.tsv")
    v.to_csv(out, sep="\t", index=False)
    print(f"\n[✔] Per-linkage table written to {out}")


if __name__ == "__main__":
    main()
