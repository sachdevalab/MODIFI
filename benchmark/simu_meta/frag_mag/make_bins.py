#!/usr/bin/env python3
"""make_bins.py — build a MODIFI --bin_file for one (completeness, contamination) cell.

Each host isolate's MAG bin = a random subset of its fragments summing to `completeness`%
of the genome (incompleteness), PLUS foreign fragments from other isolates summing to
`contamination`% of the retained length (contamination). ECE contigs are NOT binned (they
are the query). A foreign fragment is simply also listed under the contaminated host's bin
- no sequence duplication; MODIFI's merge_bin_profile aggregates it into that bin's profile.

Note: load_bin() reads with header=0, so the file MUST start with a header row.
Usage: make_bins.py <fragment_map.tsv> <out.bin_file> <completeness%> <contamination%> [seed]
"""
import sys
import pandas as pd


def main():
    fmap, out, comp, contam = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 42
    m = pd.read_csv(fmap, sep="\t")
    lines = ["contig\tbin_id"]                         # header (load_bin uses header=0)
    for iso, g in m.groupby("isolate"):
        frags = g.sample(frac=1, random_state=seed).reset_index(drop=True)
        target, cum, kept = comp / 100.0 * g["length"].sum(), 0, []
        for _, r in frags.iterrows():
            if cum >= target and kept:
                break
            kept.append(r["frag"]); cum += r["length"]
        for f in kept:
            lines.append(f"{f}\t{iso}")
        if contam > 0:                                  # foreign fragments -> this host's bin
            others = m[m["isolate"] != iso].sample(frac=1, random_state=seed + 7)
            ctarget, cc = contam / 100.0 * cum, 0
            for _, r in others.iterrows():
                if cc >= ctarget:
                    break
                lines.append(f"{r['frag']}\t{iso}"); cc += r["length"]
    with open(out, "w") as o:
        o.write("\n".join(lines) + "\n")
    print(f"[make_bins] comp={comp:.0f}% contam={contam:.0f}% -> {out} "
          f"({len(lines)-1} contig-bin assignments, {m['isolate'].nunique()} MAGs)")


if __name__ == "__main__":
    main()
