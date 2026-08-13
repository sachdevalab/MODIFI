#!/usr/bin/env python3
"""
ECE validation evidence — batch driver.

Reads the linkage CSV (jaccard_same_sample.csv), groups ECEs by sample (`prefix`),
runs the single-sample engine for each, concatenates the per-sample tables, and
writes a cross-sample summary + rebuttal figure.

  Data tables -> /home/shuaiw/borg/revision/ece_anno/
  Figures     -> /home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/

Usage:
  python run_ece_evidence.py [--ece_csv ...] [--jobs 8] [--threads 8] [--resume]
                             [--limit N] [--summary-only]
"""
import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, "ece_evidence.py")
PY = sys.executable

DEFAULT_CSV = "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv"
BATCH_ROOT = "/home/shuaiw/borg/paper/isolation/batch2_results"
DATA_DIR = "/home/shuaiw/borg/revision/ece_anno"
FIG_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
DB_DIR = os.path.join(DATA_DIR, "db")

SUPPORT_COLS = ["support_circular", "support_genomad", "support_marker", "support_coverage"]


def run_one(sample, ece_csv, threads, resume):
    work_dir = os.path.join(BATCH_ROOT, sample)
    out_dir = os.path.join(DATA_DIR, "per_sample", sample)
    out_tsv = os.path.join(out_dir, "ece_evidence.tsv")
    if resume and os.path.exists(out_tsv):
        return sample, "skipped (resume)", out_tsv
    if not os.path.isdir(work_dir):
        return sample, "MISSING work_dir", None
    cmd = [PY, ENGINE, "--sample", sample, "--work_dir", work_dir,
           "--ece_csv", ece_csv, "--out_dir", out_dir, "--db_dir", DB_DIR,
           "--threads", str(threads)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return sample, "ok", out_tsv
    except subprocess.CalledProcessError as e:
        return sample, f"FAILED: {e.stderr.decode()[-300:]}", None


def aggregate(per_sample_dir):
    frames = []
    for tsv in sorted(glob.glob(os.path.join(per_sample_dir, "*", "ece_evidence.tsv"))):
        sample = os.path.basename(os.path.dirname(tsv))
        try:
            d = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        if d.empty:
            continue
        d.insert(0, "sample", sample)
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def summarize(all_df, out_summary):
    lines = []
    n = len(all_df)
    lines.append(f"total_ECEs\t{n}")
    lines.append(f"samples\t{all_df['sample'].nunique()}")
    lines.append("")
    lines.append("# confidence x type")
    ct = all_df.groupby(["type", "confidence"]).size().unstack(fill_value=0)
    lines.append(ct.to_csv(sep="\t").rstrip())
    lines.append("")
    lines.append("# ECEs supported by each evidence line")
    for c in SUPPORT_COLS:
        if c in all_df.columns:
            k = int(all_df[c].astype(bool).sum())
            lines.append(f"{c}\t{k}\t{k/n:.3f}")
    lines.append("")
    lines.append("# negative flags")
    for c in ["flag_chromosomal", "flag_artifact"]:
        k = int(all_df[c].astype(bool).sum())
        lines.append(f"{c}\t{k}\t{k/n:.3f}")
    lines.append("")
    ge2 = int((all_df["n_positive_lines"] >= 2).sum())
    ge1 = int((all_df["n_positive_lines"] >= 1).sum())
    lines.append(f"# HEADLINE: ECEs with >=1 independent support\t{ge1}\t{ge1/n:.3f}")
    lines.append(f"# HEADLINE: ECEs with >=2 independent supports\t{ge2}\t{ge2/n:.3f}")
    lines.append("")
    lines.append("# median cov_ratio by type")
    med = all_df.groupby("type")["cov_ratio"].median()
    lines.append(med.to_csv(sep="\t").rstrip())
    with open(out_summary, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def make_figure(all_df, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conf_order = ["High", "Medium", "Low"]
    conf_colors = {"High": "#2c7fb8", "Medium": "#7fcdbb", "Low": "#edf8b1"}
    types = sorted(all_df["type"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # (a) stacked confidence per type
    ax = axes[0]
    bottom = [0] * len(types)
    for conf in conf_order:
        vals = [int(((all_df["type"] == t) & (all_df["confidence"] == conf)).sum())
                for t in types]
        ax.bar(types, vals, bottom=bottom, label=conf, color=conf_colors[conf],
               edgecolor="white")
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_title("ECE confidence by type")
    ax.set_ylabel("number of ECEs")
    ax.legend(title="confidence", frameon=False)

    # (b) fraction of ECEs supported by each evidence line (+ negative flags)
    ax = axes[1]
    labels = ["circular", "geNomad\nrobust", "plasmid/viral\nmarker", "coverage\ndistinct",
              "chromosomal\nflag (neg)"]
    n = len(all_df)
    fracs = [all_df["support_circular"].mean(), all_df["support_genomad"].mean(),
             all_df["support_marker"].mean(), all_df["support_coverage"].mean(),
             all_df["flag_chromosomal"].mean()]
    colors = ["#2c7fb8"] * 4 + ["#d73027"]
    ax.barh(labels, fracs, color=colors)
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.set_xlabel("fraction of ECEs")
    ax.set_title(f"Evidence support (n={n})")
    for i, v in enumerate(fracs):
        ax.text(min(v + 0.02, 0.9), i, f"{v:.0%}", va="center", fontsize=9)

    # (c) distribution of independent supports
    ax = axes[2]
    counts = all_df["n_positive_lines"].value_counts().sort_index()
    ax.bar(counts.index.astype(str), counts.values, color="#41b6c4", edgecolor="white")
    ax.set_xlabel("# independent supporting lines")
    ax.set_ylabel("number of ECEs")
    ax.set_title("Independent support per ECE")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out_pdf}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ece_csv", default=DEFAULT_CSV)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="only first N samples (debug)")
    ap.add_argument("--summary-only", action="store_true",
                    help="skip running; just aggregate + summarize + plot")
    args = ap.parse_args()

    os.makedirs(os.path.join(DATA_DIR, "per_sample"), exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    if not args.summary_only:
        csv = pd.read_csv(args.ece_csv)
        samples = sorted(csv["prefix"].astype(str).unique())
        if args.limit:
            samples = samples[: args.limit]
        print(f"running {len(samples)} samples, jobs={args.jobs}, threads={args.threads}")
        done = failed = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_one, s, args.ece_csv, args.threads, args.resume): s
                    for s in samples}
            for fut in as_completed(futs):
                sample, status, _ = fut.result()
                if status.startswith("FAILED") or status.startswith("MISSING"):
                    failed += 1
                    print(f"  [{sample}] {status}", file=sys.stderr)
                else:
                    done += 1
                if (done + failed) % 25 == 0:
                    print(f"  progress {done+failed}/{len(samples)} (ok={done}, fail={failed})")
        print(f"finished: ok={done}, fail={failed}")

    all_df = aggregate(os.path.join(DATA_DIR, "per_sample"))
    if all_df.empty:
        print("no per-sample results to aggregate.", file=sys.stderr)
        return
    all_tsv = os.path.join(DATA_DIR, "ece_evidence_all.tsv")
    all_df.to_csv(all_tsv, sep="\t", index=False)
    print(f"wrote {all_tsv} ({len(all_df)} ECEs)")
    summarize(all_df, os.path.join(DATA_DIR, "ece_evidence_summary.tsv"))
    make_figure(all_df, os.path.join(FIG_DIR, "ece_evidence.pdf"))


if __name__ == "__main__":
    main()
