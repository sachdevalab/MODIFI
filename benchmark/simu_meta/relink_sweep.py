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

# ECE-side minimum motif-site cutoff sweep, WITHOUT editing the tool. estimate_linkage.py:60
# hardcodes `if p_total < 2` (the min occurrences of a motif on the ECE for it to count). We
# override linkage_score_from_counts2 at runtime with a verbatim copy whose ONLY change is that
# the cutoff is read from env var MODIFI_MIN_ECE_SITES (default 2 => identical to the shipped tool;
# ece=2 must reproduce the base run). All other logic is copied unchanged from estimate_linkage.py.
from math import log  # noqa: E402
_spec_weight = EL.specificity_weight


def _linkage_score_ece_cutoff(motif_data, min_frac, max_sites=5000):
    min_ece = int(os.environ.get("MODIFI_MIN_ECE_SITES", 2))
    scores = []; weights = []; total_sites = 0
    mge_methy_sites = 0; total_confidence = 0; total_p_meth = 0
    probability = 1; match_p = 1; motif_count = 0; miss_penalty = 0
    for m in motif_data:
        h_total = m['host_total']; h_meth = m['host_meth']
        p_total = m['plasmid_total']; p_meth = m['plasmid_meth']
        if h_total == 0:
            continue
        if p_total == 0:
            continue
        if p_total < min_ece:            # <-- only change: env-var cutoff (default 2)
            continue
        f_host = h_meth / h_total; f_plasmid = p_meth / p_total
        weight = _spec_weight(float(m["occurrence_ratio"]), s_max=6)
        if f_plasmid >= min_frac:
            motif_score = 1
            probability *= float(m["occurrence_ratio"]) * 0.01
            match_p *= float(m["occurrence_ratio"]) * 0.01
            total_sites += h_total + p_total
            motif_count += 1
            mge_methy_sites += p_meth
        else:
            motif_score = 0
            probability *= (1 - float(m["occurrence_ratio"]) * 0.01)
            miss_penalty += 2 ** (f_host - f_plasmid) - 1
        scores.append(motif_score * weight); weights.append(weight)
        total_confidence += 10000000 / m['occurrence_len']; total_p_meth += p_meth
    motif_confidence = min(log(1 + motif_count) / log(1 + 3), 1)
    confidence = min(log(1 + total_sites) / log(1 + max_sites), 1)
    if motif_count > 0:
        normalized_probability = probability ** (1.0 / motif_count)
    else:
        normalized_probability = 1.0
    linkage_score = round(match_p, 4)
    final_score = (1 - normalized_probability) * confidence * motif_confidence - miss_penalty
    final_score = max(final_score, 0)
    return {'specificity': round(linkage_score, 4), 'confidence': round(confidence, 4),
            'final_score': round(final_score, 4), 'total_sites': mge_methy_sites,
            'motif_confidence': round(motif_confidence, 4), 'host_motif_num': len(motif_data)}
EL.linkage_score_from_counts2 = _linkage_score_ece_cutoff


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work_dir", required=True, help="finished modifi/<label> dir (cached profiles)")
    ap.add_argument("--mge_file", required=True)
    ap.add_argument("--ref", required=True, help="community ref .fa (needs .fai next to it)")
    ap.add_argument("--out_dir", required=True, help="isolated output dir for this (min_frac,min_sites)")
    ap.add_argument("--min_frac", type=float, default=0.4)
    ap.add_argument("--min_sites", type=int, default=30)
    ap.add_argument("--min_ctg_cov", type=int, default=5)
    ap.add_argument("--min_ece_sites", type=int, default=2,
                    help="ECE-side per-motif min occurrences (estimate_linkage default 2)")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.environ["MODIFI_MIN_ECE_SITES"] = str(a.min_ece_sites)  # read by the overridden linkage fn

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

    print(f"[relink] {a.out_dir} min_frac={a.min_frac} min_sites={a.min_sites} "
          f"min_ctg_cov={a.min_ctg_cov} min_ece_sites={a.min_ece_sites} threads={a.threads}")
    batch_MGE_invade(a.mge_file, os.path.join(a.work_dir, "profiles"), hosts, a.ref,
                     bin_file=None, min_frac=a.min_frac, threads=a.threads,
                     min_ctg_cov=a.min_ctg_cov, min_detect=a.min_sites)
    if os.path.exists(out_summary):
        print(f"[done] wrote {out_summary}")
    else:
        print(f"[warn] {out_summary} not produced")


if __name__ == "__main__":
    main()
