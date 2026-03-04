#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


work_dir = "/home/shuaiw/borg/paper/gg_run3/"


@dataclass(frozen=True)
class SamplePaths:
    sample: str
    sample_dir: Path
    gtdb_bac120_summary_tsv: Path
    motif_profile_csv: Path


def iter_sample_dirs(work_dir_path: Path) -> Iterable[Path]:
    for p in sorted(work_dir_path.iterdir()):
        if p.is_dir():
            yield p


def find_gtdb_bac120_summary(sample_dir: Path) -> Path | None:
    # Most runs: {sample}/GTDB/gtdbtk.bac120.summary.tsv (may be symlink)
    direct = sample_dir / "GTDB" / "gtdbtk.bac120.summary.tsv"
    if direct.exists():
        return direct

    # Fallback patterns in case the symlink/layout differs
    candidates = list(sample_dir.glob("**/GTDB/**/gtdbtk.bac120.summary.tsv"))
    if candidates:
        return candidates[0]
    candidates = list(sample_dir.glob("**/gtdbtk.bac120.summary.tsv"))
    if candidates:
        return candidates[0]
    return None


def find_motif_profile(sample_dir: Path, motif_profile_name: str) -> Path | None:
    # Prefer methylation4-style outputs
    preferred = list(sample_dir.glob(f"**/*methylation4*/{motif_profile_name}"))
    if preferred:
        return preferred[0]

    # Otherwise, accept any motif_profile.csv under the sample
    any_hits = list(sample_dir.glob(f"**/{motif_profile_name}"))
    if any_hits:
        return any_hits[0]
    return None


def discover_samples(work_dir_path: Path, motif_profile_name: str) -> list[SamplePaths]:
    samples: list[SamplePaths] = []
    for sample_dir in iter_sample_dirs(work_dir_path):
        gtdb = find_gtdb_bac120_summary(sample_dir)
        motif = find_motif_profile(sample_dir, motif_profile_name)
        if gtdb is None or motif is None:
            continue
        samples.append(
            SamplePaths(
                sample=sample_dir.name,
                sample_dir=sample_dir,
                gtdb_bac120_summary_tsv=gtdb,
                motif_profile_csv=motif,
            )
        )
    return samples


def read_patescibacteriota_genomes(
    gtdb_bac120_summary_tsv: Path,
    phylum_token: str = "p__Patescibacteriota",
) -> dict[str, str]:
    """
    Returns {genome_id: taxonomy_string} for genomes matching phylum_token.
    """
    hits: dict[str, str] = {}
    with gtdb_bac120_summary_tsv.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genome = (row.get("user_genome") or "").strip()
            if not genome:
                continue
            # classification is the final GTDB taxonomy assignment; fall back to pplacer_taxonomy
            taxonomy = (row.get("classification") or row.get("pplacer_taxonomy") or "").strip()
            if phylum_token in taxonomy:
                hits[genome] = taxonomy
    return hits


def read_all_taxonomy(summary_tsv: Path) -> dict[str, str]:
    """
    Returns {genome_id: taxonomy_string} for all genomes in a GTDB summary TSV.
    """
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
    """
    Extract phylum from GTDB taxonomy string (e.g. d__Bacteria;p__Proteobacteria;... -> p__Proteobacteria).
    Returns None if no phylum-level rank is present.
    """
    if not taxonomy or not taxonomy.strip():
        return None
    for part in taxonomy.strip().split(";"):
        part = part.strip()
        if part.startswith("p__") and len(part) > 3:
            return part
    return None


def is_unclassified(taxonomy: str) -> bool:
    """True if taxonomy is empty, N/A, or indicates Unclassified (Bacteria/Archaea)."""
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
    """
    Returns (genome_ids, presence_matrix) where presence_matrix is shape (motifs, genomes) bool.

    Assumes motif_profile.csv has:
      - first column: motif_identifier
      - remaining columns: genomes, numeric values
    """
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
    """
    For each target genome, compute Jaccard similarity vs all genomes.

    Returns rows of (target, other, jaccard, intersection, union).
    """
    if presence.dtype != np.bool_:
        presence = presence.astype(bool, copy=False)

    id_to_idx = {g: i for i, g in enumerate(genome_ids)}
    # motifs x genomes
    counts = presence.sum(axis=0).astype(np.int32, copy=False)

    results: list[tuple[str, str, float, int, int]] = []
    for tgt in target_genomes:
        tgt_idx = id_to_idx.get(tgt)
        if tgt_idx is None:
            continue

        tgt_vec = presence[:, tgt_idx]
        tgt_count = int(counts[tgt_idx])

        # intersection with all genomes: for boolean, dot is intersection count
        inter = (tgt_vec[:, None] & presence).sum(axis=0).astype(np.int32, copy=False)
        union = tgt_count + counts - inter

        # avoid divide-by-zero (e.g., genomes with no present motifs)
        with np.errstate(divide="ignore", invalid="ignore"):
            jac = np.where(union > 0, inter / union, 0.0).astype(np.float32, copy=False)

        hit_idx = np.where(jac >= jaccard_threshold)[0]
        for j in hit_idx:
            if j == tgt_idx:
                continue
            results.append((tgt, genome_ids[j], float(jac[j]), int(inter[j]), int(union[j])))

    # sort: target then decreasing similarity
    results.sort(key=lambda r: (r[0], -r[2], r[1]))
    return results


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Collect GTDB genomes annotated as p__Patescibacteriota and, "
            "for each, find genomes with Jaccard similarity over a threshold "
            "based on motif_profile.csv presence/absence."
        )
    )
    ap.add_argument("--work-dir", default=work_dir, help="Directory containing per-sample folders")
    ap.add_argument(
        "--motif-profile-name",
        default="motif_profile.csv",
        help="Motif profile filename to search for under each sample",
    )
    ap.add_argument(
        "--phylum-token",
        default="p__Patescibacteriota",
        help="Substring token to match within GTDB taxonomy",
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
        help="Output TSV path (default: <work-dir>/patescibacteriota_jaccard_ge{thr}.tsv)",
    )

    args = ap.parse_args()
    work_dir_path = Path(args.work_dir)
    if not work_dir_path.exists():
        raise FileNotFoundError(work_dir_path)

    out_path = (
        Path(args.out)
        if args.out
        else work_dir_path / f"patescibacteriota_jaccard_ge{args.jaccard_threshold}.tsv"
    )

    samples = discover_samples(work_dir_path, args.motif_profile_name)
    if not samples:
        raise RuntimeError(
            f"No samples found under {work_dir_path} with both GTDB and {args.motif_profile_name}"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Collect (sample, other_genome, other_taxonomy) for linked-genome summary
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
                "gtdb_bac120_summary_tsv",
            ]
        )

        for sp in samples:
            patesci = read_patescibacteriota_genomes(sp.gtdb_bac120_summary_tsv, args.phylum_token)
            if not patesci:
                continue

            # Collect taxonomy for all genomes in both bacterial and archaeal summaries (if present)
            all_taxonomy: dict[str, str] = {}
            all_taxonomy.update(read_all_taxonomy(sp.gtdb_bac120_summary_tsv))
            ar53_summary = sp.gtdb_bac120_summary_tsv.with_name("gtdbtk.ar53.summary.tsv")
            if ar53_summary.exists():
                all_taxonomy.update(read_all_taxonomy(ar53_summary))

            genome_ids, presence = load_presence_matrix(sp.motif_profile_csv, args.presence_threshold)
            hits = compute_jaccard_hits_for_targets(
                genome_ids=genome_ids,
                presence=presence,
                target_genomes=patesci.keys(),
                jaccard_threshold=args.jaccard_threshold,
            )

            for target, other, jac, inter, union in hits:
                other_tax = all_taxonomy.get(other, "")
                w.writerow(
                    [
                        sp.sample,
                        target,
                        all_taxonomy.get(target, patesci.get(target, "")),
                        other,
                        other_tax,
                        f"{jac:.6f}",
                        inter,
                        union,
                        str(sp.motif_profile_csv),
                        str(sp.gtdb_bac120_summary_tsv),
                    ]
                )
                linked_per_sample.append((sp.sample, other, other_tax))

    # Summary of linked genomes: exclude Patescibacteriota and unclassified; only phylum-annotated
    summary_path = out_path.with_name(
        out_path.stem + "_linked_phyla_summary.tsv"
    )
    # (sample, phylum) -> list of (genome_id, full_taxonomy)
    TAX_SEP = " ;; "
    phylum_entries: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    for sample, other_genome, other_tax in linked_per_sample:
        if args.phylum_token in other_tax:
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