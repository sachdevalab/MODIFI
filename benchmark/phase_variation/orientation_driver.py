#!/usr/bin/env python3
"""Stage 2b: read-level orientation of the invertible TRD across time points.

For each multi-SP Type I candidate:
  - choose one reference contig (prefer a circular '_C' member),
  - build a raw reference and an in-place-inverted reference (the TRD region
    reverse-complemented at the same coordinates),
  - realign every time point's reads (from that time point's per-contig filtered
    BAM) to BOTH references with minimap2 map-hifi,
  - count reads that cleanly cross the region boundaries in each orientation,
  - orientation ratio = inv / (inv + raw) per time point.

A shifting ratio across time points is the read-level evidence that the element
actually inverts in vivo (phase variation), generalizing the E. faecalis workup.

HEAVY: run under SLURM (full node). No subprocess timeouts (SLURM manages walltime).

Outputs to /home/shuaiw/borg/revision/phase_variation/orientation/:
  refs/<cluster>/{raw.fa,inv.fa}, aln/<cluster>/<tp>.{raw,inv}.bam,
  orientation_ratios.csv (combined)
"""

import csv
import os
import re
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
import pysam

RUN2 = Path("/groups/banfield/projects/multienv/methylation/data/borg/paper/run2")
BASE = Path("/home/shuaiw/borg/revision/phase_variation")
REGIONS = BASE / "candidate_inversion_regions.tsv"
OUT = BASE / "orientation"
THREADS = int(os.environ.get("SLURM_CPUS_ON_NODE", "8"))


def sh(cmd, **kw):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True, **kw)


def meth_dir(sample):
    for m in ("methylation4", "methylation3"):
        d = RUN2 / sample / f"{sample}_{m}"
        if d.exists():
            return d
    return None


def filtered_bam(sample, contig):
    d = meth_dir(sample)
    for name in (f"{contig}.filtered.bam", f"{contig}.bam"):
        b = d / "bams" / name
        if b.exists():
            return b
    return None


def contig_fa(sample, contig):
    return meth_dir(sample) / "contigs" / f"{contig}.fa"


def build_inverted(ref_fa, out_raw, out_inv, contig, start, end):
    """out_raw = single-record raw ref; out_inv = same but region [start,end] RC'd."""
    rec = next(SeqIO.parse(str(ref_fa), "fasta"))
    seq = str(rec.seq)
    rec.id = contig
    rec.description = ""
    SeqIO.write(rec, str(out_raw), "fasta")
    s0, e0 = start - 1, end  # 0-based half-open
    inv = seq[:s0] + str(Seq(seq[s0:e0]).reverse_complement()) + seq[e0:]
    rec.seq = Seq(inv)
    SeqIO.write(rec, str(out_inv), "fasta")


def align(reads_bam, ref_fa, out_bam):
    fq = out_bam.with_suffix(".fq")
    sh(["samtools", "fastq", "-@", str(THREADS), str(reads_bam)],
       stdout=open(fq, "w"))
    p1 = subprocess.Popen(
        ["minimap2", "-t", str(THREADS), "-ax", "map-hifi", str(ref_fa), str(fq)],
        stdout=subprocess.PIPE)
    sh(["samtools", "sort", "-@", str(THREADS), "-o", str(out_bam)], stdin=p1.stdout)
    p1.wait()
    sh(["samtools", "index", str(out_bam)])
    fq.unlink(missing_ok=True)


def count_crossing(bam, contig, start, end, pad=20):
    """Count reads whose single primary alignment spans the WHOLE invertible
    segment colinearly (reference_start < start-pad AND reference_end > end+pad).

    Full-span is required (not per-junction): the segment is flanked by inverted
    repeats, so a read crossing only one junction aligns colinearly on BOTH the
    raw and inverted references (the IR arms match either way) and cannot
    discriminate orientation. Only a read that traverses the entire segment maps
    colinearly on the reference matching its internal orientation and clips/splits
    on the other. The segment (~SP span, ~3-6 kb) is below the median HiFi read
    length, so full-span reads are plentiful.
    """
    n = 0
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for r in bf.fetch(contig, max(0, start - 200), end + 200):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.reference_start < start - pad and r.reference_end > end + pad:
                n += 1
    return n


def choose_ref_row(rows):
    circ = [r for r in rows if r["contig"].endswith("_C")]
    pool = circ or rows
    # widest SP span (the invertible specificity segment), then latest timepoint
    return max(pool, key=lambda r: (int(r["sp_span_end"]) - int(r["sp_span_start"]),
                                    int(r["timepoint"])))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(REGIONS), delimiter="\t"))
    by_cluster = {}
    for r in rows:
        by_cluster.setdefault(r["cluster"], []).append(r)

    only = sys.argv[1] if len(sys.argv) > 1 else None  # optional single-cluster run
    results = []
    for cluster, crows in by_cluster.items():
        if only and cluster != only:
            continue
        ref = choose_ref_row(crows)
        contig = ref["contig"]
        # invertible segment = the SP (specificity) span; count reads crossing its junctions
        rs, re_ = int(ref["sp_span_start"]), int(ref["sp_span_end"])
        species = ref["species"]
        print(f"\n=== {cluster} {species}: ref {contig}:{rs}-{re_} "
              f"(tp {ref['timepoint']}, {ref['sample']}) ===", flush=True)

        rdir = OUT / "refs" / cluster
        rdir.mkdir(parents=True, exist_ok=True)
        raw_fa, inv_fa = rdir / "raw.fa", rdir / "inv.fa"
        build_inverted(contig_fa(ref["sample"], contig), raw_fa, inv_fa,
                       contig, rs, re_)
        for fa in (raw_fa, inv_fa):
            if not Path(str(fa) + ".fai").exists():
                sh(["samtools", "faidx", str(fa)])

        adir = OUT / "aln" / cluster
        adir.mkdir(parents=True, exist_ok=True)
        # one alignment per time point (dedup timepoints; use its own contig's reads)
        seen = {}
        for r in sorted(crows, key=lambda x: int(x["timepoint"])):
            tp = r["timepoint"]
            if tp in seen:
                continue
            seen[tp] = r
            rb = filtered_bam(r["sample"], r["contig"])
            if rb is None:
                print(f"  [skip tp{tp}] no read bam for {r['contig']}", flush=True)
                continue
            raw_bam = adir / f"tp{tp}.raw.bam"
            inv_bam = adir / f"tp{tp}.inv.bam"
            align(rb, raw_fa, raw_bam)
            align(rb, inv_fa, inv_bam)
            raw_c = count_crossing(raw_bam, contig, rs, re_)
            inv_c = count_crossing(inv_bam, contig, rs, re_)
            ratio = inv_c / (inv_c + raw_c) if (inv_c + raw_c) else float("nan")
            print(f"  tp{tp} {r['sample']}: raw={raw_c} inv={inv_c} inv_frac={ratio:.4f}",
                  flush=True)
            results.append(dict(cluster=cluster, species=species, subject=ref["subject"],
                                timepoint=tp, sample=r["sample"], ref_contig=contig,
                                region_start=rs, region_end=re_,
                                raw_crossing=raw_c, inv_crossing=inv_c,
                                inv_fraction=ratio))

    out_csv = OUT / ("orientation_ratios.csv" if not only
                     else f"orientation_ratios_{only}.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()) if results else
                           ["cluster"])
        w.writeheader()
        for row in results:
            w.writerow(row)
    print(f"\nWrote {out_csv} ({len(results)} rows)")


if __name__ == "__main__":
    main()
