#!/usr/bin/env python3
"""
Coverage evaluation: subsample aligned BAMs to different depths, run methylation
detection, and record motif counts vs depth for selected contigs.

Usage:
  1. Select 10 contigs from motif_num_all_samples.csv with depth > min_depth (default 3) from run2 mean_depth.csv.
  2. For each contig, subsample its BAM to target depths (e.g. 5, 10, 20, 50, full).
  3. Run MODIFI main.py (methylation detection) for each (contig, depth) with sample-specific control.
  4. Count detected motifs per run and output a summary table. The depth column in the summary
     is read from each run's modification result folder (work_dir/mean_depth.csv).
"""

import os
import sys
import random
import subprocess
import argparse
from pathlib import Path

import pandas as pd

# Import motif loading and dereplication from isolation/sample_object
_benchmark_dir = Path(__file__).resolve().parent.parent
if str(_benchmark_dir) not in sys.path:
    sys.path.insert(0, str(_benchmark_dir))
from isolation.sample_object import get_unique_motifs

# Default paths
MOTIF_CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/motif_num_all_samples.csv"
ALL_DIR = "/home/shuaiw/borg/paper/run2"
WORK_BASE = "/home/shuaiw/borg/paper/coverage"
MAIN_PY = "/home/shuaiw/MODIFI/main.py"

# Target depths to evaluate (only those <= original depth will be run)
DEFAULT_TARGET_DEPTHS = [5, 10, 20, 50]  # and "full" added per contig

# MODIFI parameters (match call_methy.sh)
MIN_LEN = 1000
MIN_COV = 0
MIN_FRAC = 0.3
MIN_SCORE = 30
MIN_SITES = 100
READ_TYPE = "hifi"
THREADS = 10


def load_contigs_with_depth(motif_csv: str, all_dir: str, min_depth: float, seed: int):
    """
    Load motif CSV, join with mean_depth per sample, keep rows with depth > min_depth.
    Return unique (sample, contig, depth) and randomly select 10.
    """
    df = pd.read_csv(motif_csv)
    # Columns: sample, motif_num, environment, contig, ...
    if "sample" not in df.columns or "contig" not in df.columns:
        raise ValueError(f"Expected columns 'sample' and 'contig' in {motif_csv}")

    # Unique (sample, contig) from motif CSV
    pairs = df[["sample", "contig"]].drop_duplicates()
    # Load mean_depth once per sample
    all_depths = []
    for sample in pairs["sample"].unique():
        sample_dir = Path(all_dir) / sample / f"{sample}_methylation4"
        depth_file = sample_dir / "mean_depth.csv"
        if not depth_file.exists():
            continue
        depth_df = pd.read_csv(depth_file)
        if "contig" not in depth_df.columns or "depth" not in depth_df.columns:
            continue
        depth_df = depth_df[["contig", "depth"]].copy()
        depth_df["sample"] = sample
        all_depths.append(depth_df)
    if not all_depths:
        raise SystemExit("No mean_depth.csv files found under all_dir.")
    depth_merged = pd.concat(all_depths, ignore_index=True)
    # Join: keep (sample, contig) that exist in both and depth > min_depth
    merged = pairs.merge(depth_merged, on=["sample", "contig"], how="inner")
    candidates_df = merged[merged["depth"] > min_depth].drop_duplicates()
    if candidates_df.empty:
        raise SystemExit("No contigs found with depth > min_depth.")

    random.seed(seed)
    n_select = min(10, len(candidates_df))
    selected = candidates_df.sample(n=n_select, random_state=seed)
    return selected[["sample", "contig", "depth"]]


def get_sample_paths(all_dir: str, sample: str):
    """Return bams_dir, contigs_dir, control_dir, and paths to kmer mean/num dat for sample."""
    base = Path(all_dir) / sample / f"{sample}_methylation4"
    bams_dir = base / "bams"
    contigs_dir = base / "contigs"
    control_dir = base / "control"
    kmer_mean = control_dir / "control_db.up7.down3.mean.dat"
    kmer_num = control_dir / "control_db.up7.down3.num.dat"
    return {
        "bams_dir": bams_dir,
        "contigs_dir": contigs_dir,
        "control_dir": control_dir,
        "kmer_mean_db": kmer_mean,
        "kmer_num_db": kmer_num,
    }


def subsample_bam(input_bam: str, output_bam: str, original_depth: float, target_depth: float, seed: int = 42):
    """
    Subsample BAM so expected depth is target_depth.
    Uses samtools view --subsample FRAC --subsample-seed INT.
    """
    if target_depth >= original_depth:
        frac = 1.0
    else:
        frac = target_depth / original_depth
    if frac < 0.01:
        frac = 0.01
    cmd = ["samtools", "view", "-b", "--subsample", str(frac), "--subsample-seed", str(seed), "-o", output_bam, input_bam]
    subprocess.run(cmd, check=True)
    # Index for pbmm2/pbindex usage if needed; MODIFI split_bam may use pbindex
    subprocess.run(["samtools", "index", output_bam], check=True)
    pbi = Path(output_bam + ".pbi")
    if not pbi.exists():
        # Try pbindex if available (required by MODIFI for some steps)
        try:
            subprocess.run(["pbindex", output_bam], check=True, capture_output=True)
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass


def run_methylation(
    work_dir: str,
    whole_bam: str,
    whole_ref: str,
    kmer_mean_db: str,
    kmer_num_db: str,
    main_py: str = MAIN_PY,
    threads: int = THREADS,
    min_ctg_cov: int = 1,
):
    """Run MODIFI main.py with --whole_bam and given ref + control DBs."""
    os.makedirs(work_dir, exist_ok=True)
    cmd = [
        sys.executable,
        main_py,
        "--work_dir", work_dir,
        "--whole_bam", whole_bam,
        "--whole_ref", whole_ref,
        "--read_type", READ_TYPE,
        "--min_len", str(MIN_LEN),
        "--min_cov", str(MIN_COV),
        "--min_frac", str(MIN_FRAC),
        "--min_score", str(MIN_SCORE),
        "--min_sites", str(MIN_SITES),
        "--min_ctg_cov", str(min_ctg_cov),
        "--kmer_mean_db", kmer_mean_db,
        "--kmer_num_db", kmer_num_db,
        "--threads", str(threads),
    ]
    subprocess.run(cmd, check=True)


def count_motifs(work_dir: str, contig: str, min_frac: float = MIN_FRAC, min_sites: int = MIN_SITES) -> int:
    """Load motifs from work_dir/motifs/{contig}.motifs.csv, dereplicate (unique motifs only, reverse-complement aware), and return count. Uses get_unique_motifs from isolation/sample_object."""
    motif_file = Path(work_dir) / "motifs" / f"{contig}.motifs.csv"
    if not motif_file.exists():
        return 0
    df = pd.read_csv(motif_file)
    if df.empty:
        return 0
    n_unique, _, _ = get_unique_motifs(df, min_frac=min_frac, min_sites=min_sites)
    return n_unique


def get_depth_from_result(work_dir: str, contig: str):
    """Read actual mean depth from modification result folder (work_dir/mean_depth.csv). Returns float or None."""
    depth_file = Path(work_dir) / "mean_depth.csv"
    if not depth_file.exists():
        return None
    df = pd.read_csv(depth_file)
    if "contig" not in df.columns or "depth" not in df.columns:
        return None
    match = df[df["contig"] == contig]
    if match.empty:
        return None
    return float(match["depth"].iloc[0])


def _q(s: str) -> str:
    """Quote and escape for double-quoted shell string."""
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def build_sbatch_wrap(
    *,
    input_bam: str,
    contig_fa: str,
    work_dir: str,
    kmer_mean_db: str,
    kmer_num_db: str,
    original_depth: float,
    target_d: int,
    contig: str,
    work_base: str,
    run_name: str,
    seed: int,
    main_py: str,
    threads: int,
) -> str:
    """Build the shell commands for one sbatch job: subsample (if needed) then run main.py."""
    parts = []
    if target_d < original_depth:
        frac = target_d / original_depth
        if frac < 0.01:
            frac = 0.01
        subsampled_dir = os.path.join(work_base, "subsampled", run_name)
        use_bam = os.path.join(subsampled_dir, f"{contig}.bam")
        parts.append(f"mkdir -p {_q(subsampled_dir)}")
        parts.append(f"samtools view -b --subsample {frac} --subsample-seed {seed} -o {_q(use_bam)} {_q(input_bam)}")
        parts.append(f"samtools index {_q(use_bam)}")
    else:
        use_bam = input_bam
    py_cmd = (
        f"python {_q(main_py)}"
        f" --work_dir {_q(work_dir)}"
        f" --whole_bam {_q(use_bam)}"
        f" --whole_ref {_q(contig_fa)}"
        f" --read_type {READ_TYPE}"
        f" --min_len {MIN_LEN}"
        f" --min_cov 3"
        f" --min_frac {MIN_FRAC}"
        f" --min_score {MIN_SCORE}"
        f" --min_sites {MIN_SITES}"
        f" --min_ctg_cov 3"
        f" --kmer_mean_db {_q(kmer_mean_db)}"
        f" --kmer_num_db {_q(kmer_num_db)}"
        f" --threads {threads}"
    )
    parts.append(py_cmd)
    return " && ".join(parts)


def run_test_depth3(args):
    """Run a single task: soil_1_743_C at target depth 3 (subsample + methylation + count)."""
    sample = "soil_1"
    contig = "soil_1_743_C"
    target_d = 3
    work_base = Path(args.work_base)
    work_base.mkdir(parents=True, exist_ok=True)

    paths = get_sample_paths(args.all_dir, sample)
    for key in ["kmer_mean_db", "kmer_num_db"]:
        if not paths[key].exists():
            print(f"[test_depth3] Missing: {paths[key]}", file=sys.stderr)
            sys.exit(1)
    input_bam = paths["bams_dir"] / f"{contig}.bam"
    contig_fa = paths["contigs_dir"] / f"{contig}.fa"
    if not input_bam.exists():
        print(f"[test_depth3] BAM not found: {input_bam}", file=sys.stderr)
        sys.exit(1)
    if not contig_fa.exists():
        print(f"[test_depth3] FASTA not found: {contig_fa}", file=sys.stderr)
        sys.exit(1)

    depth_file = paths["control_dir"].parent / "mean_depth.csv"
    if not depth_file.exists():
        print(f"[test_depth3] mean_depth.csv not found: {depth_file}", file=sys.stderr)
        sys.exit(1)
    depth_df = pd.read_csv(depth_file)
    match = depth_df[depth_df["contig"] == contig]
    if match.empty:
        print(f"[test_depth3] Contig {contig} not in mean_depth.csv", file=sys.stderr)
        sys.exit(1)
    original_depth = float(match["depth"].iloc[0])

    run_name = f"{contig}_depth{target_d}"
    work_dir = work_base / run_name
    subsampled_dir = work_base / "subsampled" / run_name
    subsampled_dir.mkdir(parents=True, exist_ok=True)
    use_bam = str(subsampled_dir / f"{contig}.bam")

    print(f"[test_depth3] Subsampling to depth {target_d} (original {original_depth:.1f})...")
    subsample_bam(str(input_bam), use_bam, original_depth, target_d, args.seed)
    print(f"[test_depth3] Running methylation detection...")
    run_methylation(
        str(work_dir),
        use_bam,
        str(contig_fa),
        str(paths["kmer_mean_db"]),
        str(paths["kmer_num_db"]),
        main_py=args.main_py,
        threads=args.threads,
        min_ctg_cov=1,
    )
    depth_from_result = get_depth_from_result(str(work_dir), contig)
    motif_count = count_motifs(str(work_dir), contig)
    print(f"[test_depth3] Done. depth (from result): {depth_from_result}, unique motif count: {motif_count}")


def main():
    ap = argparse.ArgumentParser(description="Coverage evaluation: subsample BAMs, run methylation, count motifs vs depth.")
    ap.add_argument("--motif_csv", default=MOTIF_CSV, help="Path to motif_num_all_samples.csv")
    ap.add_argument("--all_dir", default=ALL_DIR, help="Base dir for run2 sample dirs (e.g. run2/96plex/96plex_methylation4)")
    ap.add_argument("--work_base", default=WORK_BASE, help="Base work dir for coverage runs")
    ap.add_argument("--min_depth", type=float, default=100, help="Only consider contigs with depth > this")
    ap.add_argument("--target_depths", type=str, default="3,5,10,20,50,100",
                    help="Comma-separated target depths (e.g. 5,10,20,50). 'full' is always added per contig.")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for selecting contigs and subsampling")
    ap.add_argument("--threads", type=int, default=THREADS)
    ap.add_argument("--skip_run", action="store_true", help="Only collect motif counts from existing runs (no subsample/run)")
    ap.add_argument("--generate_sbatch", action="store_true", help="Write a shell script with one sbatch job per (contig, depth) instead of running in Python")
    ap.add_argument("-o", "--output_script", type=str, default=None, help="Output path for generated sbatch script (default: work_base/run_coverage_jobs.sh)")
    ap.add_argument("--partition", type=str, default="standard", help="Slurm partition for sbatch")
    ap.add_argument("--list_contigs", action="store_true", help="Only print the 10 selected contigs and exit")
    ap.add_argument("--test_depth3", action="store_true", help="Small test: run only soil_1_743_C at depth 3 (subsample + methylation + count)")
    ap.add_argument("--main_py", default=MAIN_PY)
    args = ap.parse_args()

    # Small test: single contig soil_1_743_C at depth 3
    if args.test_depth3:
        run_test_depth3(args)
        return 0

    target_depths_list = DEFAULT_TARGET_DEPTHS
    if args.target_depths:
        target_depths_list = [int(x.strip()) for x in args.target_depths.split(",")]

    selected = load_contigs_with_depth(args.motif_csv, args.all_dir, args.min_depth, args.seed)
    if args.list_contigs:
        print(selected.to_string(index=False))
        return 0

    work_base = Path(args.work_base)
    work_base.mkdir(parents=True, exist_ok=True)

    # Build list of (contig, depth) tasks for either sbatch generation or run
    tasks = []
    for _, row in selected.iterrows():
        sample = row["sample"]
        contig = row["contig"]
        original_depth = row["depth"]
        paths = get_sample_paths(args.all_dir, sample)
        if not paths["kmer_mean_db"].exists() or not paths["kmer_num_db"].exists():
            print(f"Skip {contig}: control DBs not found for sample {sample}", file=sys.stderr)
            continue
        input_bam = paths["bams_dir"] / f"{contig}.bam"
        contig_fa = paths["contigs_dir"] / f"{contig}.fa"
        if not input_bam.exists():
            print(f"Skip {contig}: BAM not found {input_bam}", file=sys.stderr)
            continue
        if not contig_fa.exists():
            print(f"Skip {contig}: FASTA not found {contig_fa}", file=sys.stderr)
            continue
        depths_to_run = [d for d in target_depths_list if d <= original_depth]
        if original_depth not in depths_to_run:
            depths_to_run.append(int(round(original_depth)))
        depths_to_run = sorted(set(depths_to_run))
        for target_d in depths_to_run:
            tasks.append({
                "sample": sample,
                "contig": contig,
                "original_depth": original_depth,
                "target_d": target_d,
                "input_bam": str(input_bam),
                "contig_fa": str(contig_fa),
                "kmer_mean_db": str(paths["kmer_mean_db"]),
                "kmer_num_db": str(paths["kmer_num_db"]),
                "run_name": f"{contig}_depth{target_d}",
                "work_dir": str(work_base / f"{contig}_depth{target_d}"),
            })

    if args.generate_sbatch:
        out_script = Path(args.output_script or work_base / "run_coverage_jobs.sh")
        with open(out_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("# Generated by coverage_evaluation.py --generate_sbatch\n")
            f.write(f"# Submit with: bash {out_script.name}\n\n")
            for t in tasks:
                wrap = build_sbatch_wrap(
                    input_bam=t["input_bam"],
                    contig_fa=t["contig_fa"],
                    work_dir=t["work_dir"],
                    kmer_mean_db=t["kmer_mean_db"],
                    kmer_num_db=t["kmer_num_db"],
                    original_depth=t["original_depth"],
                    target_d=t["target_d"],
                    contig=t["contig"],
                    work_base=str(work_base),
                    run_name=t["run_name"],
                    seed=args.seed,
                    main_py=args.main_py,
                    threads=args.threads,
                )
                job_name = f"cov_{t['contig']}_d{t['target_d']}".replace(".", "_")
                wrap_escaped = wrap.replace("\\", "\\\\").replace('"', '\\"')
                f.write(f'sbatch --partition {args.partition} --wrap "{wrap_escaped}" --job-name={job_name}\n')
        print(f"Wrote {len(tasks)} sbatch commands to {out_script}")
        print(f"Submit with: bash {out_script}")
        return 0

    results = []
    for t in tasks:
        contig = t["contig"]
        sample = t["sample"]
        original_depth = t["original_depth"]
        target_d = t["target_d"]
        work_dir = t["work_dir"]
        if target_d < original_depth:
            use_bam = str(work_base / "subsampled" / t["run_name"] / f"{contig}.bam")
            if not args.skip_run:
                subsample_bam(t["input_bam"], use_bam, original_depth, target_d, args.seed)
            elif not Path(use_bam).exists():
                print(f"Skip {t['run_name']}: subsampled BAM missing (run without --skip_run)", file=sys.stderr)
                continue
        else:
            use_bam = t["input_bam"]
        if not args.skip_run:
            run_methylation(
                work_dir,
                use_bam,
                t["contig_fa"],
                t["kmer_mean_db"],
                t["kmer_num_db"],
                main_py=args.main_py,
                threads=args.threads,
                min_ctg_cov=1,
            )
        motif_count = count_motifs(work_dir, contig)
        depth_from_result = get_depth_from_result(work_dir, contig)
        results.append({
            "contig": contig,
            "sample": sample,
            "original_depth": original_depth,
            "target_depth": target_d,
            "depth": depth_from_result,
            "motif_count": motif_count,
        })
        depth_str = f"{depth_from_result:.2f}" if depth_from_result is not None else "NA"
        print(f"{contig}\tsample={sample}\ttarget_depth={target_d}\tdepth={depth_str}\tmotifs={motif_count}")

    out_df = pd.DataFrame(results)
    out_file = work_base / "motif_count_vs_depth.csv"
    out_df.to_csv(out_file, index=False)
    print(f"Results written to {out_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
