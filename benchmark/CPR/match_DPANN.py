#!/usr/bin/env python3
"""
DPANN archaea analysis: same workflow as CPR (match_CPR.py) but for DPANN.
DPANN are Archaea; genomes are taken from GTDB ar53 (archaeal) summary, not bac120.
Reports which GTDB phylum names are treated as DPANN (align with your GTDB version).
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


work_dir = "/home/shuaiw/borg/paper/gg_run3/"

# GTDB phylum names (p__*) that belong to the DPANN superphylum.
# Align this list with your GTDB release (see GTDB tree for d__Archaea).
DPANN_PHYLA_GTDB = [
    "p__Aenigmarchaeota",
    "p__Altarchaeota",
    "p__Diapherotrites",
    "p__Hadarchaeota",
    "p__Huberarchaeota",
    "p__Mamarchaeota",
    "p__Micrarchaeota",
    "p__Nanoarchaeota",
    "p__Nanobdellota",
    "p__Nanohaloarchaeota",
    "p__Pacearchaeota",
    "p__Parvarchaeota",
    "p__Undinarchaeota",
    "p__Woesearchaeota",
]


@dataclass(frozen=True)
class SamplePaths:
    sample: str
    sample_dir: Path
    gtdb_ar53_summary_tsv: Path
    gtdb_bac120_summary_tsv: Path | None
    motif_profile_csv: Path


def iter_sample_dirs(work_dir_path: Path) -> Iterable[Path]:
    for p in sorted(work_dir_path.iterdir()):
        if p.is_dir():
            yield p


def find_gtdb_ar53_summary(sample_dir: Path) -> Path | None:
    direct = sample_dir / "GTDB" / "gtdbtk.ar53.summary.tsv"
    if direct.exists():
        return direct
    for pattern in ("**/GTDB/**/gtdbtk.ar53.summary.tsv", "**/gtdbtk.ar53.summary.tsv"):
        candidates = list(sample_dir.glob(pattern))
        if candidates:
            return candidates[0]
    return None


def find_gtdb_bac120_summary(sample_dir: Path) -> Path | None:
    direct = sample_dir / "GTDB" / "gtdbtk.bac120.summary.tsv"
    if direct.exists():
        return direct
    for pattern in ("**/GTDB/**/gtdbtk.bac120.summary.tsv", "**/gtdbtk.bac120.summary.tsv"):
        candidates = list(sample_dir.glob(pattern))
        if candidates:
            return candidates[0]
    return None


def find_motif_profile(sample_dir: Path, motif_profile_name: str) -> Path | None:
    preferred = list(sample_dir.glob(f"**/*methylation4*/{motif_profile_name}"))
    if preferred:
        return preferred[0]
    any_hits = list(sample_dir.glob(f"**/{motif_profile_name}"))
    if any_hits:
        return any_hits[0]
    return None


def discover_samples(work_dir_path: Path, motif_profile_name: str) -> list[SamplePaths]:
    samples: list[SamplePaths] = []
    for sample_dir in iter_sample_dirs(work_dir_path):
        ar53 = find_gtdb_ar53_summary(sample_dir)
        motif = find_motif_profile(sample_dir, motif_profile_name)
        if ar53 is None or motif is None:
            continue
        bac120 = find_gtdb_bac120_summary(sample_dir)
        samples.append(
            SamplePaths(
                sample=sample_dir.name,
                sample_dir=sample_dir,
                gtdb_ar53_summary_tsv=ar53,
                gtdb_bac120_summary_tsv=bac120,
                motif_profile_csv=motif,
            )
        )
    return samples


def read_dpann_genomes(
    gtdb_ar53_summary_tsv: Path,
    dpann_phyla: list[str],
) -> dict[str, str]:
    """Return {genome_id: taxonomy_string} for genomes in any of the given DPANN phyla."""
    hits: dict[str, str] = {}
    with gtdb_ar53_summary_tsv.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genome = (row.get("user_genome") or "").strip()
            if not genome:
                continue
            taxonomy = (row.get("classification") or row.get("pplacer_taxonomy") or "").strip()
            if any(ph in taxonomy for ph in dpann_phyla):
                hits[genome] = taxonomy
    return hits


def read_all_taxonomy(summary_tsv: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with summary_tsv.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genome = (row.get("user_genome") or "").strip()
            if not genome:
                continue
            taxonomy = (row.get("classification") or row.get("pplacer_taxonomy") or "").strip()
            if taxonomy:
                mapping[genome] = taxonomy
    return mapping


def extract_phylum(taxonomy: str) -> str | None:
    if not taxonomy or not taxonomy.strip():
        return None
    for part in taxonomy.strip().split(";"):
        part = part.strip()
        if part.startswith("p__") and len(part) > 3:
            return part
    return None


def is_unclassified(taxonomy: str) -> bool:
    if not taxonomy or not taxonomy.strip():
        return True
    t = taxonomy.strip()
    if t.upper() in ("N/A", "NA", ""):
        return True
    if t.startswith("Unclassified"):
        return True
    return False


def load_presence_matrix(
    motif_profile_csv: Path,
    presence_threshold: float,
) -> tuple[list[str], np.ndarray]:
    df = pd.read_csv(motif_profile_csv)
    if df.shape[1] < 2:
        raise ValueError(f"Unexpected motif_profile.csv format: {motif_profile_csv}")
    genome_ids = [c for c in df.columns if c != "motif_identifier"]
    values = df[genome_ids].to_numpy(dtype=np.float32, copy=False)
    presence = values > presence_threshold
    return genome_ids, presence


def compute_jaccard_hits_for_targets(
    genome_ids: list[str],
    presence: np.ndarray,
    target_genomes: Iterable[str],
    jaccard_threshold: float,
) -> list[tuple[str, str, float, int, int]]:
    if presence.dtype != np.bool_:
        presence = presence.astype(bool, copy=False)
    id_to_idx = {g: i for i, g in enumerate(genome_ids)}
    counts = presence.sum(axis=0).astype(np.int32, copy=False)
    results: list[tuple[str, str, float, int, int]] = []
    for tgt in target_genomes:
        tgt_idx = id_to_idx.get(tgt)
        if tgt_idx is None:
            continue
        tgt_vec = presence[:, tgt_idx]
        tgt_count = int(counts[tgt_idx])
        inter = (tgt_vec[:, None] & presence).sum(axis=0).astype(np.int32, copy=False)
        union = tgt_count + counts - inter
        with np.errstate(divide="ignore", invalid="ignore"):
            jac = np.where(union > 0, inter / union, 0.0).astype(np.float32, copy=False)
        hit_idx = np.where(jac >= jaccard_threshold)[0]
        for j in hit_idx:
            if j == tgt_idx:
                continue
            results.append((tgt, genome_ids[j], float(jac[j]), int(inter[j]), int(union[j])))
    results.sort(key=lambda r: (r[0], -r[2], r[1]))
    return results


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Collect GTDB archaeal genomes in DPANN phyla (from ar53 summary) and, "
            "for each, find genomes with Jaccard similarity over a threshold. "
            "Reports GTDB phylum names treated as DPANN."
        )
    )
    ap.add_argument("--work-dir", default=work_dir, help="Directory containing per-sample folders")
    ap.add_argument(
        "--motif-profile-name",
        default="motif_profile.csv",
        help="Motif profile filename to search for under each sample",
    )
    ap.add_argument(
        "--dpann-phyla",
        nargs="*",
        default=None,
        help="GTDB phylum names (p__X) to treat as DPANN; default: built-in list",
    )
    ap.add_argument(
        "--presence-threshold",
        type=float,
        default=0.3,
        help="motif_profile value must be > this to count as present",
    )
    ap.add_argument(
        "--jaccard-threshold",
        type=float,
        default=0.8,
        help="Report genome pairs with Jaccard similarity >= this",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output TSV path (default: <work-dir>/dpann_jaccard_ge{thr}.tsv)",
    )
    ap.add_argument(
        "--report-dpann-phyla",
        default=None,
        help="Write list of GTDB DPANN phylum names to this file (default: next to main TSV)",
    )

    args = ap.parse_args()
    work_dir_path = Path(args.work_dir)
    if not work_dir_path.exists():
        raise FileNotFoundError(work_dir_path)

    dpann_phyla = args.dpann_phyla if args.dpann_phyla else DPANN_PHYLA_GTDB

    out_path = (
        Path(args.out)
        if args.out
        else work_dir_path / f"dpann_jaccard_ge{args.jaccard_threshold}.tsv"
    )

    samples = discover_samples(work_dir_path, args.motif_profile_name)
    if not samples:
        raise RuntimeError(
            f"No samples found under {work_dir_path} with both GTDB ar53 and {args.motif_profile_name}"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    report_path = args.report_dpann_phyla or out_path.with_name("dpann_GTDB_phyla_used.txt")
    with open(report_path, "w") as f:
        f.write("GTDB phylum names treated as DPANN (archaeal ar53):\n")
        for p in sorted(dpann_phyla):
            f.write(f"  {p}\n")
        f.write("\nAlign this list with your GTDB release if needed.\n")
    print(f"DPANN phyla (GTDB): {', '.join(sorted(dpann_phyla))}", file=sys.stderr)
    print(f"Wrote: {report_path}")

    linked_per_sample: list[tuple[str, str, str]] = []

    with out_path.open("w", newline="") as out_f:
        w = csv.writer(out_f, delimiter="\t")
        w.writerow(
            [
                "sample",
                "target_genome",
                "target_taxonomy",
                "other_genome",
                "other_taxonomy",
                "jaccard",
                "intersection_motifs",
                "union_motifs",
                "motif_profile_csv",
                "gtdb_ar53_summary_tsv",
            ]
        )

        for sp in samples:
            dpann = read_dpann_genomes(sp.gtdb_ar53_summary_tsv, dpann_phyla)
            if not dpann:
                continue

            all_taxonomy: dict[str, str] = {}
            all_taxonomy.update(read_all_taxonomy(sp.gtdb_ar53_summary_tsv))
            if sp.gtdb_bac120_summary_tsv and sp.gtdb_bac120_summary_tsv.exists():
                all_taxonomy.update(read_all_taxonomy(sp.gtdb_bac120_summary_tsv))

            genome_ids, presence = load_presence_matrix(sp.motif_profile_csv, args.presence_threshold)
            hits = compute_jaccard_hits_for_targets(
                genome_ids=genome_ids,
                presence=presence,
                target_genomes=dpann.keys(),
                jaccard_threshold=args.jaccard_threshold,
            )

            for target, other, jac, inter, union in hits:
                other_tax = all_taxonomy.get(other, "")
                w.writerow(
                    [
                        sp.sample,
                        target,
                        all_taxonomy.get(target, dpann.get(target, "")),
                        other,
                        other_tax,
                        f"{jac:.6f}",
                        inter,
                        union,
                        str(sp.motif_profile_csv),
                        str(sp.gtdb_ar53_summary_tsv),
                    ]
                )
                linked_per_sample.append((sp.sample, other, other_tax))

    summary_path = out_path.with_name(out_path.stem + "_linked_phyla_summary.tsv")
    TAX_SEP = " ;; "
    phylum_entries: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for sample, other_genome, other_tax in linked_per_sample:
        if any(ph in other_tax for ph in dpann_phyla):
            continue
        if is_unclassified(other_tax):
            continue
        phylum = extract_phylum(other_tax)
        if phylum is None:
            continue
        key = (sample, phylum)
        phylum_entries[key].append((other_genome, other_tax))

    with summary_path.open("w", newline="") as out_f:
        w = csv.writer(out_f, delimiter="\t")
        w.writerow(["sample", "phylum", "count", "genome_ids", "taxonomy_strings"])
        for (sample, phylum), entries in sorted(phylum_entries.items()):
            entries_sorted = sorted(entries, key=lambda e: e[0])
            genome_ids = ";".join(e[0] for e in entries_sorted)
            taxonomy_strings = TAX_SEP.join(e[1] for e in entries_sorted)
            w.writerow([sample, phylum, len(entries_sorted), genome_ids, taxonomy_strings])

    print(f"Wrote: {out_path}")
    print(f"Wrote: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
