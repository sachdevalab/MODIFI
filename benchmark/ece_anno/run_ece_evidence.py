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

DATA_DIR = "/home/shuaiw/borg/revision/ece_anno"
FIG_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
DB_DIR = os.path.join(DATA_DIR, "db")

# Per-dataset config: CSV, sample-dir root, and CSV column mapping.
DATASETS = {
    "isolate": {
        "csv": "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv",
        "batch_root": "/home/shuaiw/borg/paper/isolation/batch2_results",
        "criteria": "strict",
        "cols": {"sample": "prefix", "seqname": "mge_contig", "type": "mge_type",
                 "length": "mge_length", "mgedepth": "mge_depth",
                 "hostdepth": "host_depth", "lineage": "host_lineage"},
    },
    "metagenome": {
        "csv": "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv",
        "batch_root": "/home/shuaiw/borg/paper/run2",
        "criteria": "loose",
        "cols": {"sample": "sample", "seqname": "MGE", "type": "MGE_type",
                 "length": "mge_len", "mgedepth": "MGE_cov",
                 "hostdepth": "host_cov", "lineage": "host_taxa"},
    },
    # ALL depth>=5 & length>=1kb plasmid/virus ECEs in the 59 metagenomes (linked + unlinked)
    "metagenome_all": {
        "csv": "/home/shuaiw/borg/revision/ece_anno/all_metagenome_eces.csv",
        "batch_root": "/home/shuaiw/borg/paper/run2",
        "criteria": "loose",
        "cols": {"sample": "sample", "seqname": "MGE", "type": "MGE_type",
                 "length": "mge_len", "mgedepth": "MGE_cov",
                 "hostdepth": "host_cov", "lineage": "host_taxa"},
    },
    # Multi-caller expansion: the combined-set ECEs (>=5x) that are NOT already in the
    # metagenome_all evidence (NEW non-geNomad calls + geNomad calls from the 5 samples
    # not covered before). Evidence computed here; the "any-strong-caller" P2 is applied
    # in post (needs the `methods` column), so raw criteria=loose just fills the columns.
    "expanded_new": {
        "csv": "/home/shuaiw/borg/revision/ece_anno/expanded/needs_evidence.csv",
        "batch_root": "/home/shuaiw/borg/paper/run2",
        "criteria": "loose",
        "cols": {"sample": "sample", "seqname": "MGE", "type": "MGE_type",
                 "length": "mge_len", "mgedepth": "MGE_cov",
                 "hostdepth": "host_cov", "lineage": "host_taxa"},
    },
}

SUPPORT_COLS = ["support_circular", "support_genomad", "support_marker", "support_coverage"]
DROP_LEGACY = ["confidence", "confidence_reason"]


def reclassify(d, criteria="strict"):
    """(Re)derive very_high_confidence from the stored support/flag columns.
    strict = all four support lines; loose = geNomad + marker only. In loose mode the
    chromosomal flag drops the SCMG component (rRNA-only), recomputed from stored counts."""
    for c in SUPPORT_COLS + ["flag_chromosomal", "flag_artifact"]:
        d[c] = d[c].astype(bool)
    if criteria == "loose":
        flag_chr = pd.to_numeric(d.get("rrna_count"), errors="coerce").fillna(0) > 0
        vh = (d["support_genomad"] & d["support_marker"]
              & ~flag_chr & ~d["flag_artifact"])
        d["flag_chromosomal"] = flag_chr
        d["very_high_confidence"] = vh
    else:
        d["very_high_confidence"] = (d[SUPPORT_COLS].all(axis=1)
                                     & ~d["flag_chromosomal"] & ~d["flag_artifact"])
    return d.drop(columns=[c for c in DROP_LEGACY if c in d.columns])


def run_one(sample, ece_csv, batch_root, per_sample_root, cols, threads, resume,
            criteria="strict"):
    work_dir = os.path.join(batch_root, sample)
    out_dir = os.path.join(per_sample_root, sample)
    out_tsv = os.path.join(out_dir, "ece_evidence.tsv")
    if resume and os.path.exists(out_tsv):
        return sample, "skipped (resume)", out_tsv
    if not os.path.isdir(work_dir):
        return sample, "MISSING work_dir", None
    cmd = [PY, ENGINE, "--sample", sample, "--work_dir", work_dir,
           "--ece_csv", ece_csv, "--out_dir", out_dir, "--db_dir", DB_DIR,
           "--threads", str(threads), "--criteria", criteria]
    for k, v in cols.items():
        cmd += [f"--col_{k}", v]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        return sample, "ok", out_tsv
    except subprocess.CalledProcessError as e:
        return sample, f"FAILED: {e.stderr.decode()[-300:]}", None


def aggregate(per_sample_dir, rewrite=True, criteria="strict"):
    frames = []
    for tsv in sorted(glob.glob(os.path.join(per_sample_dir, "*", "ece_evidence.tsv"))):
        sample = os.path.basename(os.path.dirname(tsv))
        try:
            d = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        if d.empty:
            continue
        d = reclassify(d, criteria)
        if rewrite:  # normalize the per-sample table to the new schema on disk
            d.to_csv(tsv, sep="\t", index=False)
        d.insert(0, "sample", sample)
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def summarize(all_df, out_summary):
    lines = []
    n = len(all_df)
    vh = all_df["very_high_confidence"].astype(bool)
    lines.append(f"total_ECEs\t{n}")
    lines.append(f"samples\t{all_df['sample'].nunique()}")
    lines.append("")
    lines.append("# very_high_confidence (all 4 lines pass, no negative flag)")
    lines.append(f"very_high_confidence\t{int(vh.sum())}\t{vh.mean():.3f}")
    lines.append("")
    lines.append("# very_high_confidence by type")
    ct = all_df.assign(very_high_confidence=vh).groupby(
        ["type", "very_high_confidence"]).size().unstack(fill_value=0)
    lines.append(ct.to_csv(sep="\t").rstrip())
    lines.append("")
    lines.append("# ECEs supported by each evidence line (fraction)")
    for c in SUPPORT_COLS:
        k = int(all_df[c].astype(bool).sum())
        lines.append(f"{c}\t{k}\t{k/n:.3f}")
    lines.append("")
    lines.append("# negative flags")
    for c in ["flag_chromosomal", "flag_artifact"]:
        k = int(all_df[c].astype(bool).sum())
        lines.append(f"{c}\t{k}\t{k/n:.3f}")
    lines.append("")
    lines.append("# distribution of n_positive_lines (0-4)")
    lines.append(all_df["n_positive_lines"].value_counts().sort_index().to_csv(sep="\t").rstrip())
    lines.append("")
    lines.append("# median cov_ratio by type")
    lines.append(all_df.groupby("type")["cov_ratio"].median().to_csv(sep="\t").rstrip())
    with open(out_summary, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def make_figure(all_df, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    all_df = all_df.copy()
    all_df["very_high_confidence"] = all_df["very_high_confidence"].astype(bool)
    types = sorted(all_df["type"].unique())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # (a) very-high-confidence vs rest, per type
    ax = axes[0]
    vh = [int(((all_df["type"] == t) & all_df["very_high_confidence"]).sum()) for t in types]
    rest = [int(((all_df["type"] == t) & ~all_df["very_high_confidence"]).sum()) for t in types]
    ax.bar(types, vh, label="very-high (all 4 lines)", color="#2c7fb8", edgecolor="white")
    ax.bar(types, rest, bottom=vh, label="not very-high", color="#d9d9d9", edgecolor="white")
    ax.set_title("Very-high-confidence ECEs by type")
    ax.set_ylabel("number of ECEs")
    ax.legend(frameon=False)

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
    ap.add_argument("--dataset", choices=list(DATASETS), default="isolate",
                    help="which ECE set: isolate (default) or metagenome")
    ap.add_argument("--ece_csv", default=None, help="override the dataset's CSV")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="only first N samples (debug)")
    ap.add_argument("--summary-only", action="store_true",
                    help="skip running; just aggregate + summarize + plot")
    args = ap.parse_args()

    cfg = DATASETS[args.dataset]
    ece_csv = args.ece_csv or cfg["csv"]
    batch_root = cfg["batch_root"]
    cols = cfg["cols"]
    criteria = cfg.get("criteria", "strict")
    # isolate keeps the original DATA_DIR layout; other datasets get a subdir.
    out_root = DATA_DIR if args.dataset == "isolate" else os.path.join(DATA_DIR, args.dataset)
    per_sample_root = os.path.join(out_root, "per_sample")
    fig_name = "ece_evidence.pdf" if args.dataset == "isolate" \
        else f"ece_evidence_{args.dataset}.pdf"
    os.makedirs(per_sample_root, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)

    if not args.summary_only:
        csv = pd.read_csv(ece_csv)
        samples = sorted(csv[cols["sample"]].astype(str).unique())
        if args.limit:
            samples = samples[: args.limit]
        print(f"[{args.dataset}] running {len(samples)} samples, "
              f"jobs={args.jobs}, threads={args.threads}")
        done = failed = 0
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_one, s, ece_csv, batch_root, per_sample_root, cols,
                              args.threads, args.resume, criteria): s for s in samples}
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

    all_df = aggregate(per_sample_root, criteria=criteria)
    if all_df.empty:
        print("no per-sample results to aggregate.", file=sys.stderr)
        return
    all_tsv = os.path.join(out_root, "ece_evidence_all.tsv")
    all_df.to_csv(all_tsv, sep="\t", index=False)
    print(f"wrote {all_tsv} ({len(all_df)} ECEs)")
    summarize(all_df, os.path.join(out_root, "ece_evidence_summary.tsv"))
    make_figure(all_df, os.path.join(FIG_DIR, fig_name))


if __name__ == "__main__":
    main()
