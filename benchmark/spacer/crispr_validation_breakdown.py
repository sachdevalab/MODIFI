#!/usr/bin/env python3
"""Break down the metagenome ECE-host CRISPR-spacer validation rate.

The overall figure pools populations CRISPR interrogates very differently, so we
report it split by MGE type (plasmid vs virus) and by environment, and on two
denominators:
  - ALL linkages
  - only linkages whose HOST encodes a CRISPR array (the only hosts that could
    possibly yield spacer-based validation)

A linkage (MGE, host contig) is counted as CRISPR-validated when the linked host
contig carries a CRISPR spacer that targets that MGE. This uses the same spacer
hits (My_sample.read_spacer, 0 mismatch) as compare_spacer.py. The linkage set is
the FINAL strict 317-linkage table (revised ECE set).
"""
import os
import sys
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "isolation"))
from sample_object import My_sample

ALL_DIR = "/home/shuaiw/borg/paper/run2/"
# FINAL strict linkage set (revised ECEs) -> 317 linkages
CSV = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"


def hosts_with_array(prefix, all_dir=ALL_DIR):
    """Set of host contigs that encode >=1 CRISPR array (minced spacer sources)."""
    fa = f"{all_dir}/{prefix}/spacer/{prefix}_spacers.fa"
    s = set()
    if not os.path.exists(fa):
        return s
    for line in open(fa):
        if line.startswith(">"):
            s.add(line[1:].split("_CRISPR_")[0])
    return s


def validate_linkages(csv=CSV, all_dir=ALL_DIR, mismatch_allowed=0):
    """Return a per-linkage DataFrame with `validated` and `host_has_array` flags."""
    df = pd.read_csv(csv)
    rows = []
    for prefix, grp in df.groupby("sample"):
        obj = My_sample(prefix, all_dir)
        # MGE -> set(host contigs carrying a spacer that hits this MGE)
        spacer = obj.read_spacer(mismatch_allowed=mismatch_allowed)
        arrays = hosts_with_array(prefix, all_dir)
        for _, r in grp.iterrows():
            mge, host = r["MGE"], r["host"]
            validated = int(mge in spacer and host in spacer[mge])
            rows.append((prefix, r["environment"], r["type"], mge, host,
                         int(host in arrays), validated))
    return pd.DataFrame(
        rows,
        columns=["sample", "environment", "type", "MGE", "host",
                 "host_has_array", "validated"],
    )


def _fmt(label, sub):
    n = len(sub)
    k = int(sub["validated"].sum())
    rate = f"{100 * k / n:5.2f}%" if n else "  n/a"
    return f"{label:<28} {k:>3}/{n:<4} = {rate}"


def _block(title, v):
    out = [title, _fmt("all linkages", v)]
    for t in ["plasmid", "virus"]:
        out.append(_fmt("  " + t, v[v.type == t]))
    return out


def summarize(v):
    arr = v[v.host_has_array == 1]
    lines = ["=" * 60]
    lines += _block(f"DENOMINATOR = ALL linkages ({len(v)})", v)
    lines.append("")
    lines += _block(
        f"DENOMINATOR = host encodes a CRISPR array ({len(arr)}/{len(v)})", arr)
    lines.append("")
    lines.append("BY ENVIRONMENT (all-linkage denominator)")
    for env in sorted(v.environment.unique()):
        lines.append(_fmt(env, v[v.environment == env]))
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
