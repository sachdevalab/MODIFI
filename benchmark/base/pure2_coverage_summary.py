#!/usr/bin/env python3
"""
Count motif number and depth per contig in 96plex subsample folders under pure2.
- Base dir: /home/shuaiw/borg/paper/linkage/pure2
- Process each subfolder except m64004_210929_143746.p10_test
- Depth from mean_depth.csv, motif count via get_unique_motifs (unique, dereplicated)
- Only keep contigs that have at least 1 motif in p100 (discard contigs with no motif even in p100)
"""

import sys
from pathlib import Path

import pandas as pd

_benchmark_dir = Path(__file__).resolve().parent.parent
if str(_benchmark_dir) not in sys.path:
    sys.path.insert(0, str(_benchmark_dir))
from isolation.sample_object import get_unique_motifs

PURE2_BASE = "/home/shuaiw/borg/paper/linkage/pure2"
EXCLUDE_FOLDER = "m64004_210929_143746.p10_test"
MIN_FRAC = 0.6   # only count motifs with fraction > 0.6
MIN_SITES = 500  # and modified sites (nDetected) > 500


def get_folders(base: Path, exclude: str):
    """Return sorted list of subdirs that look like run folders (have motifs + mean_depth), excluding exclude."""
    out = []
    for p in base.iterdir():
        if not p.is_dir() or p.name == exclude:
            continue
        if (p / "mean_depth.csv").exists() and (p / "motifs").is_dir():
            out.append(p)
    return sorted(out, key=lambda x: x.name)


def count_motifs_unique(motif_file: Path, min_frac=MIN_FRAC, min_sites=MIN_SITES) -> int:
    """Count unique motifs (dereplicated) using sample_object.get_unique_motifs."""
    if not motif_file.exists():
        return 0
    df = pd.read_csv(motif_file)
    if df.empty:
        return 0
    n, _, _ = get_unique_motifs(df, min_frac=min_frac, min_sites=min_sites)
    return n


def get_depth_for_contig(mean_depth_path: Path, contig: str):
    """Return depth for contig from mean_depth.csv or None."""
    if not mean_depth_path.exists():
        return None
    df = pd.read_csv(mean_depth_path)
    if "contig" not in df.columns or "depth" not in df.columns:
        return None
    match = df[df["contig"] == contig]
    if match.empty:
        return None
    return float(match["depth"].iloc[0])


def main():
    base = Path(PURE2_BASE)
    if not base.exists():
        print(f"Base dir not found: {base}", file=sys.stderr)
        sys.exit(1)

    folders = get_folders(base, EXCLUDE_FOLDER)
    p100_name = "m64004_210929_143746.p100"
    p100_path = base / p100_name
    if not p100_path.exists() or not (p100_path / "motifs").exists():
        print(f"p100 folder or motifs not found: {p100_path}", file=sys.stderr)
        sys.exit(1)

    # Contigs that have at least 1 motif in p100
    motif_dir_p100 = p100_path / "motifs"
    contigs_in_p100 = set()
    contigs_with_motif_in_p100 = []
    for f in motif_dir_p100.glob("*.motifs.csv"):
        contig = f.name.replace(".motifs.csv", "")
        contigs_in_p100.add(contig)
        n = count_motifs_unique(f, min_frac=MIN_FRAC, min_sites=MIN_SITES)
        if n >= 1:
            contigs_with_motif_in_p100.append(contig)
    contigs_keep = sorted(set(contigs_with_motif_in_p100))
    # Keep only one contig per species: the one with _1 (e.g. B_cepacia_UCB-717_1)
    # contigs_keep = [c for c in contigs_keep if c.endswith("_1")]
    print(f"Contigs with ≥1 motif in p100 (one per species, _1 only): {len(contigs_keep)}", file=sys.stderr)

    rows = []
    for folder in folders:
        name = folder.name
        mean_depth_file = folder / "mean_depth.csv"
        motifs_dir = folder / "motifs"
        for contig in contigs_keep:
            depth = get_depth_for_contig(mean_depth_file, contig)
            motif_file = motifs_dir / f"{contig}.motifs.csv"
            n_motifs = count_motifs_unique(motif_file, min_frac=MIN_FRAC, min_sites=MIN_SITES)
            rows.append({
                "folder": name,
                "contig": contig,
                "depth": depth,
                "motif_count": n_motifs,
            })

    df = pd.DataFrame(rows)
    out_file = base / "coverage_motif_summary.csv"
    df.to_csv(out_file, index=False)
    print(f"Wrote {len(df)} rows to {out_file}")

    # Pivot tables for readability (optional print)
    depth_pivot = df.pivot(index="contig", columns="folder", values="depth")
    motif_pivot = df.pivot(index="contig", columns="folder", values="motif_count")
    depth_file = base / "coverage_depth_pivot.csv"
    motif_file = base / "coverage_motif_count_pivot.csv"
    depth_pivot.to_csv(depth_file)
    motif_pivot.to_csv(motif_file)
    print(f"Depth pivot: {depth_file}")
    print(f"Motif count pivot: {motif_file}")


if __name__ == "__main__":
    main()
