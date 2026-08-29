#!/usr/bin/env python3
"""Per-ECE zero-coverage fraction for the isolate ECEs, from batch2_results/<s>/<s>.align.bam.
zero_cov_frac = fraction of contig positions with depth 0. Parallel across samples (<=64 procs).
Writes isolate_all/zerocov.csv (sample,contig,zero_cov_frac)."""
import csv, os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pysam

B = "/home/shuaiw/borg/paper/isolation/batch2_results"
OUT = "/home/shuaiw/borg/revision/ece_anno/isolate_all"
CSV = f"{OUT}/isolate_all_eces.csv"

by_sample = defaultdict(list)
for r in csv.DictReader(open(CSV)):
    by_sample[r["prefix"]].append(r["mge_contig"])

def one_sample(sample):
    bam = f"{B}/{sample}/{sample}.align.bam"
    out = []
    if not os.path.exists(bam):
        return [(sample, c, "") for c in by_sample[sample]]
    try:
        af = pysam.AlignmentFile(bam, "rb")
        refs = set(af.references)
        for c in by_sample[sample]:
            if c not in refs:
                out.append((sample, c, "")); continue
            L = af.get_reference_length(c)
            cov = np.zeros(L, dtype=np.int32)
            for col in af.pileup(c, truncate=True, min_base_quality=0, stepper="nofilter"):
                cov[col.reference_pos] = col.nsegments
            zf = float((cov == 0).sum()) / L if L else ""
            out.append((sample, c, f"{zf:.4f}" if zf != "" else ""))
        af.close()
    except Exception:
        out = [(sample, c, "") for c in by_sample[sample]]
    return out

rows = []
with ProcessPoolExecutor(max_workers=64) as ex:
    futs = {ex.submit(one_sample, s): s for s in by_sample}
    done = 0
    for f in as_completed(futs):
        rows.extend(f.result()); done += 1
        if done % 100 == 0:
            print(f"  {done}/{len(by_sample)} samples")

with open(f"{OUT}/zerocov.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample", "contig", "zero_cov_frac"])
    w.writerows(rows)
n_hi = sum(1 for _, _, z in rows if z and float(z) > 0.10)
print(f"wrote {OUT}/zerocov.csv ({len(rows)} ECEs); zero_cov>0.10: {n_hi}")
