#!/usr/bin/env python3
"""
Step 2 (v2 whole-metagenome): select contigs with length >= 100 kb AND HiFi depth
>= 10x from the MODIFI HiFi full-metagenome run, and emit the per-contig bam paths
for the downstream per-contig SOTA runs.

Reads MODIFI's own mean_depth.csv (cols: contig,depth,length) from the HiFi run;
cross-references the subread run's depth and both runs' per-contig bams.

Writes: compare_all_meta/selected_contigs.tsv
"""
from __future__ import annotations
import csv, os, sys

OUT = "/home/shuaiw/borg/paper/ipdsummary/compare_all_meta"
MODIFI_HIFI = f"{OUT}/modifi_hifi_full"
MODIFI_SUB = f"{OUT}/modifi_subread_full"
MIN_LEN = 100_000
MIN_HIFI_DEPTH = 10.0


def load_depth(path):
    """contig -> (depth, length) from a MODIFI mean_depth.csv."""
    d = {}
    if not os.path.isfile(path):
        return d
    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                d[row["contig"]] = (float(row["depth"]), int(float(row["length"])))
            except (KeyError, ValueError):
                continue
    return d


def main():
    hifi = load_depth(f"{MODIFI_HIFI}/mean_depth.csv")
    if not hifi:
        sys.exit(f"ERROR: no {MODIFI_HIFI}/mean_depth.csv yet (MODIFI HiFi run not finished?)")
    sub = load_depth(f"{MODIFI_SUB}/mean_depth.csv")

    rows = []
    for ctg, (hdep, hlen) in hifi.items():
        if hlen < MIN_LEN or hdep < MIN_HIFI_DEPTH:
            continue
        hbam = f"{MODIFI_HIFI}/bams/{ctg}.bam"
        sbam = f"{MODIFI_SUB}/bams/{ctg}.bam"
        sdep = sub.get(ctg, (None, None))[0]
        rows.append({
            "contig": ctg,
            "length": hlen,
            "hifi_depth": round(hdep, 2),
            "subread_depth": round(sdep, 2) if sdep is not None else "NA",
            "is_circular": ctg.endswith("_C"),
            "hifi_bam": hbam if os.path.isfile(hbam) else "MISSING",
            "subread_bam": sbam if os.path.isfile(sbam) else "MISSING",
        })
    rows.sort(key=lambda r: -r["hifi_depth"])

    out = f"{OUT}/selected_contigs.tsv"
    fields = ["contig", "length", "hifi_depth", "subread_depth", "is_circular",
              "hifi_bam", "subread_bam"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    nc = sum(1 for r in rows if r["is_circular"])
    miss_h = sum(1 for r in rows if r["hifi_bam"] == "MISSING")
    miss_s = sum(1 for r in rows if r["subread_bam"] == "MISSING")
    print(f"Wrote {out}: {len(rows)} contigs (len>={MIN_LEN}, HiFi depth>={MIN_HIFI_DEPTH}x); "
          f"{nc} circular, {len(rows)-nc} linear.")
    if miss_h or miss_s:
        print(f"  WARNING: {miss_h} missing HiFi bams, {miss_s} missing subread bams "
              f"(MODIFI --no-clean kept bams only for contigs with depth>=min_ctg_cov).")


if __name__ == "__main__":
    main()
