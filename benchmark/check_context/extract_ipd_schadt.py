"""
Extract per-position IPD from a PacBio subreads BAM using the *Schadt* normalization
(for the k-mer context variance-decomposition analysis), NOT MODIFI's standard
normalization.

Reuses only the BAM-reading mechanics from scripts/standard_load7.py (pbcore
AlignmentSet, 1/FrameRate scaling, readsInRange chunking, the aligned-position
`matched` mask). The normalization is replaced with:
    1. rawIpd = aln.IPD() * (1 / FrameRate)          # seconds
    2. logIpd = log(rawIpd + 0.01)
    3. subtract the per-subread mean of logIpd        # remove per-molecule speed
    then, per genomic position:
    4. Tukey (boxplot) outlier removal
    5. keep positions with >= --min-reads reads
    6. y_i = mean of the remaining centered log-IPD

Run in the pbcore env:
    conda run -n MODIFI_subreads python benchmark/check_context/extract_ipd_schadt.py

Output: a compact CSV cache (strand, tpl, y, n), one row per kept incorporation site.
"""

import argparse
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# ---- defaults (match the file the analysis was opened from) ------------------
DEFAULT_REF = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
DEFAULT_BAM = "/home/shuaiw/borg/paper/base/pure/control/bams/CP011331.1.bam"
DEFAULT_OUT = "/home/shuaiw/borg/revision/context/c224/c227_ipd_schadt.csv"

# ---- hit filters (same thresholds as standard_load7._loadRawIpds) ------------
MAP_QV_THRESHOLD = 0
MIN_IDENTITY = 0.0
MIN_READLENGTH = 50
MAX_REGION = 100000  # bp per chunk, as in standard_load7


def schadt_position_means(s_dict, min_reads):
    """Aggregate {pos: [centered log-IPD, ...]} -> list of (pos, y, n).

    Tukey boxplot outlier removal, then require n >= min_reads.
    """
    out = []
    for pos, vals in s_dict.items():
        d = np.asarray(vals, dtype=np.float64)
        if d.size < min_reads:
            continue
        q1, q3 = np.percentile(d, [25, 75])
        iqr = q3 - q1
        lo, hi = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        d = d[(d >= lo) & (d <= hi)]
        n = d.size
        if n < min_reads:
            continue
        out.append((pos, float(d.mean()), int(n)))
    return out


# ---- per-worker global state (set in _init_worker) ---------------------------
_W = {}


def _init_worker(bam, ref, contig, min_reads):
    from pbcore.io import AlignmentSet
    alignments = AlignmentSet(bam, referenceFastaFname=ref)
    each_ref = None
    for my_ref in alignments.referenceInfoTable:
        if contig is None or my_ref.Name == contig:
            each_ref = my_ref
            break
    _W["aln"] = alignments
    _W["factor"] = 1.0 / alignments.readGroupTable[0].FrameRate
    _W["rgid"] = alignments.referenceInfo(each_ref.Name).Name
    _W["min_reads"] = min_reads


def _worker_chunk(bounds):
    start, end = bounds
    return process_chunk(_W["aln"], _W["rgid"], start, end, _W["factor"],
                         _W["min_reads"])


def process_chunk(alignments, ref_group_id, start, end, factor, min_reads):
    """Read all subreads overlapping [start, end), Schadt-normalize, aggregate."""
    hits = [
        h for h in alignments.readsInRange(ref_group_id, max(start, 0), end)
        if (h.mapQV >= MAP_QV_THRESHOLD
            and h.identity >= MIN_IDENTITY
            and h.readLength >= MIN_READLENGTH)
    ]

    fwd, rev = {}, {}  # pos -> list of centered log-IPD
    for aln in hits:
        rawIpd = aln.IPD() * factor
        read = np.array([c != '-' for c in aln.read()])
        ref = np.array([c != '-' for c in aln.reference()])
        matched = read & ref & ~np.isnan(rawIpd)          # aligned, non-gap, valid
        if matched.sum() < 2:
            continue

        # Schadt per-subread normalization: log then center by the subread mean.
        # Center over the whole subread (all matched positions), like standard_load7
        # computes its normalization before restricting to the interval.
        log_ipd = np.log(rawIpd + 0.01)
        center = log_ipd[matched].mean()

        ref_pos = aln.referencePositions()
        in_win = matched & (ref_pos >= start) & (ref_pos < end)
        if not in_win.any():
            continue

        vals = log_ipd[in_win] - center
        tpls = ref_pos[in_win]
        target = fwd if aln.isForwardStrand else rev
        for p, v in zip(tpls, vals):
            target.setdefault(int(p), []).append(float(v))

    rows = []
    for p, y, n in schadt_position_means(fwd, min_reads):
        rows.append(("+", p, y, n))
    for p, y, n in schadt_position_means(rev, min_reads):
        rows.append(("-", p, y, n))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", default=DEFAULT_BAM)
    ap.add_argument("--ref", default=DEFAULT_REF)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--contig", default=None,
                    help="contig name; default = first contig in the BAM")
    ap.add_argument("--min-reads", type=int, default=25,
                    help="min reads per position after outlier removal (default 25)")
    ap.add_argument("--max-pos", type=int, default=None,
                    help="stop after this many bp (smoke test)")
    ap.add_argument("--threads", type=int, default=8,
                    help="parallel worker processes (default 8)")
    args = ap.parse_args()

    from pbcore.io import AlignmentSet
    t0 = time.time()
    alignments = AlignmentSet(args.bam, referenceFastaFname=args.ref)

    # resolve contig
    each_ref = None
    for my_ref in alignments.referenceInfoTable:
        if args.contig is None or my_ref.Name == args.contig:
            each_ref = my_ref
            break
    if each_ref is None:
        raise SystemExit(f"contig {args.contig} not found in {args.bam}")
    contig = each_ref.Name
    length = each_ref.Length
    del alignments  # workers re-open their own handle
    if args.max_pos is not None:
        length = min(length, args.max_pos)
    print(f"contig={contig} length={length} min_reads={args.min_reads} "
          f"threads={args.threads}", flush=True)

    n_chunks = max(1, int(length / MAX_REGION))
    bounds = []
    for i in range(n_chunks):
        start = i * MAX_REGION
        end = length if i == n_chunks - 1 else min((i + 1) * MAX_REGION, length)
        bounds.append((start, end))

    all_rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.threads, initializer=_init_worker,
                             initargs=(args.bam, args.ref, contig, args.min_reads)) as ex:
        for rows in ex.map(_worker_chunk, bounds):
            all_rows.extend(rows)
            done += 1
            print(f"  chunk {done}/{n_chunks} kept={len(rows)} total={len(all_rows)} "
                  f"t={round(time.time()-t0)}s", flush=True)

    df = pd.DataFrame(all_rows, columns=["strand", "tpl", "y", "n"])
    df.to_csv(args.out, index=False)

    # sanity summary
    n_fwd = int((df["strand"] == "+").sum())
    n_rev = int((df["strand"] == "-").sum())
    print(f"\nwrote {len(df)} sites -> {args.out}", flush=True)
    print(f"  forward={n_fwd}  reverse={n_rev}", flush=True)
    print(f"  y mean={df['y'].mean():.4f} (should be ~0)  y std={df['y'].std():.4f}",
          flush=True)
    print(f"  coverage n: median={df['n'].median()} p10={df['n'].quantile(.1):.0f} "
          f"p90={df['n'].quantile(.9):.0f}", flush=True)
    print(f"done in {round(time.time()-t0)}s", flush=True)


if __name__ == "__main__":
    main()
