"""
Compare MODIFI per-contig GFFs vs a single ipdSummary GFF on modification sites.

Per-contig metrics (intersection / union, recoveries, Jaccard) with score >= cutoff.
Only the first MAX_CONTIGS contigs (sorted by MODIFI GFF filename) are evaluated.
Results go to OUTPUT_CSV. Run:

  python benchmark/compare/gff_modifi_ipd_consistency.py

All ipd GFF feature types pass the same score filter. Sites match on position + strand
within a contig.
"""

from __future__ import annotations

import csv
import glob
import os
from typing import Dict, Optional, Set

# --- configuration ---
MODIFI_GFF_DIR = "/home/shuaiw/borg/paper/ipdsummary/soil_1/modifi.out/test_100/gffs"
IPD_GFF = "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out/test_100.gff"
SCORE_CUTOFF = 10  # keep rows with score >= SCORE_CUTOFF
MAX_CONTIGS = 30

OUTPUT_CSV = (
    "/home/shuaiw/MODIFI/tmp/figures/base_benchmark/modifi_ipd_consistency_per_contig.csv"
)


def _local_site_key(fields: list[str]) -> Optional[str]:
    """pos:strand within one contig."""
    if len(fields) < 7:
        return None
    try:
        pos = int(fields[3])
    except ValueError:
        return None
    return f"{pos}:{fields[6]}"


def _parse_score(fields: list[str]) -> Optional[float]:
    if len(fields) < 6:
        return None
    try:
        return float(fields[5])
    except ValueError:
        return None


def load_modifi_one_gff(path: str) -> Dict[str, float]:
    """pos:strand -> max score for one contig GFF."""
    best: Dict[str, float] = {}
    with open(path) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            score = _parse_score(parts)
            if score is None or score < SCORE_CUTOFF:
                continue
            key = _local_site_key(parts)
            if key is None:
                continue
            prev = best.get(key)
            if prev is None or score > prev:
                best[key] = score
    return best


def stream_ipd_for_contigs(path: str, contigs: Set[str]) -> Dict[str, Dict[str, float]]:
    """contig -> (pos:strand -> max score), only for seqnames in contigs."""
    out: Dict[str, Dict[str, float]] = {c: {} for c in contigs}
    with open(path) as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seq = parts[0]
            if seq not in contigs:
                continue
            score = _parse_score(parts)
            if score is None or score < SCORE_CUTOFF:
                continue
            key = _local_site_key(parts)
            if key is None:
                continue
            bucket = out[seq]
            prev = bucket.get(key)
            if prev is None or score > prev:
                bucket[key] = score
    return out


def _ratio(num: int, den: int) -> str:
    if den == 0:
        return ""
    return f"{num / den:.6f}"


def main() -> None:
    pattern = os.path.join(MODIFI_GFF_DIR, "*.gff")
    paths = sorted(glob.glob(pattern))[:MAX_CONTIGS]
    if not paths:
        raise SystemExit(f"No *.gff under {MODIFI_GFF_DIR!r}")

    contigs = [os.path.splitext(os.path.basename(p))[0] for p in paths]
    contig_set = set(contigs)

    modifi_by_ctg: Dict[str, Dict[str, float]] = {}
    for p, ctg in zip(paths, contigs):
        modifi_by_ctg[ctg] = load_modifi_one_gff(p)

    ipd_by_ctg = stream_ipd_for_contigs(IPD_GFF, contig_set)

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    fieldnames = [
        "contig",
        "n_modifi",
        "n_ipd",
        "n_intersection",
        "n_union",
        "modifi_recovery",
        "ipd_recovery",
        "jaccard",
    ]
    rows = []
    for ctg in contigs:
        m_keys = set(modifi_by_ctg[ctg].keys())
        i_keys = set(ipd_by_ctg[ctg].keys())
        inter = m_keys & i_keys
        union = m_keys | i_keys
        n_m, n_i, n_int, n_uni = len(m_keys), len(i_keys), len(inter), len(union)
        rows.append(
            {
                "contig": ctg,
                "n_modifi": n_m,
                "n_ipd": n_i,
                "n_intersection": n_int,
                "n_union": n_uni,
                "modifi_recovery": _ratio(n_int, n_m),
                "ipd_recovery": _ratio(n_int, n_i),
                "jaccard": _ratio(n_int, n_uni),
            }
        )

    with open(OUTPUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"SCORE_CUTOFF (>={SCORE_CUTOFF}), MAX_CONTIGS={MAX_CONTIGS}")
    print(f"Wrote {OUTPUT_CSV}")
    for r in rows:
        print(
            f"  {r['contig']}: MODIFI={r['n_modifi']} ipd={r['n_ipd']} "
            f"Jaccard={r['jaccard'] or 'NA'}"
        )


if __name__ == "__main__":
    main()
