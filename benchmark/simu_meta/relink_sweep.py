#!/usr/bin/env python3
"""relink_sweep.py — linkage-only re-run of MODIFI at a chosen (min_frac, min_sites),
reusing a finished run's cached methylation profiles. Used to calibrate the two
motif-evidence thresholds on the orphan_300 control WITHOUT recalling modifications.

Writes ONLY into an isolated --out_dir (fresh hosts/ + a copied host_summary); it
symlinks the shared inputs (motif_profile.csv, mean_depth.csv) read-only and never
touches the source run's files. If --out_dir already has host_summary.csv it exits
without doing anything (never overwrite).

Inputs reused from --work_dir (a finished modifi/<label> dir): profiles/, motifs/,
motif_profile.csv, mean_depth.csv, and the community ref (.fa + .fai) + mge_list.
"""
import argparse
import os
import shutil
import sys

sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
import estimate_linkage as EL  # noqa: E402
from estimate_linkage import batch_MGE_invade  # noqa: E402

# The finished runs auto-cleaned work_dir/contigs, which report_gc needs only to add
# cosmetic MGE_gc/host_gc/cos_sim columns. Those are unused for threshold calibration,
# and final_score/specificity/total_sites/host_motif_num are computed upstream and are
# unaffected. No-op report_gc so summary_host completes without the per-contig fastas.
def _no_gc(data, *a, **k):
    for c in ("MGE_gc", "host_gc", "cos_sim"):
        if c not in data.columns:
            data[c] = 0.0
    return data
EL.report_gc = _no_gc

# sort_top_by_cos_sim() breaks ties among 2-3 equal-top-score hosts using per-contig FASTA
# k-mer similarity; with contigs/ cleaned, kmer_freq_sim_bin_worker hits bin_1_len==0 and
# calls sys.exit(1), silently killing the whole run for any plasmid that has a top-score tie
# (get_kmer_freq.py:178-180). We don't stage contigs, so neutralize the tie-break: keep the
# tied hosts in their existing (final_score) order. Best-host / operating-point metrics are
# unaffected (tied hosts share the same final_score).
EL.sort_top_by_cos_sim = lambda data, *a, **k: data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work_dir", required=True, help="finished modifi/<label> dir (cached profiles)")
    ap.add_argument("--mge_file", required=True)
    ap.add_argument("--ref", required=True, help="community ref .fa (needs .fai next to it)")
    ap.add_argument("--out_dir", required=True, help="isolated output dir for this (min_frac,min_sites)")
    ap.add_argument("--min_frac", type=float, default=0.4)
    ap.add_argument("--min_sites", type=int, default=30)
    ap.add_argument("--min_ctg_cov", type=int, default=5)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    out_summary = os.path.join(a.out_dir, "host_summary.csv")
    if os.path.exists(out_summary):
        print(f"[skip] {out_summary} already exists, not overwriting.")
        return

    hosts = os.path.join(a.out_dir, "hosts")
    os.makedirs(hosts, exist_ok=True)
    # symlink the shared inputs the linkage step reads from out_dir (host_dir/..)
    for f in ("motif_profile.csv", "mean_depth.csv"):
        src, dst = os.path.join(a.work_dir, f), os.path.join(a.out_dir, f)
        if not os.path.exists(dst):
            os.symlink(src, dst)

    print(f"[relink] {a.out_dir} min_frac={a.min_frac} min_sites={a.min_sites} threads={a.threads}")
    batch_MGE_invade(a.mge_file, os.path.join(a.work_dir, "profiles"), hosts, a.ref,
                     bin_file=None, min_frac=a.min_frac, threads=a.threads,
                     min_ctg_cov=a.min_ctg_cov, min_detect=a.min_sites)
    if os.path.exists(out_summary):
        print(f"[done] wrote {out_summary}")
    else:
        print(f"[warn] {out_summary} not produced")


if __name__ == "__main__":
    main()
