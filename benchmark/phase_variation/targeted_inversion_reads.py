#!/usr/bin/env python3
"""Read-level verification of the assembly-flagged R-M inversions.

For each locus classified 'inversion' by rm_sv_assembly.py:
  - find the exact inverted segment (minus-strand blast block of the operon
    region against the query-timepoint contig), in ref-contig coordinates,
  - build a raw reference and an in-place-inverted reference,
  - align EVERY time point's reads to both and count reads that span the
    segment colinearly in each orientation.

Interpretation:
  - orientation flips across time points  -> population-level phase variation
  - a single time point carries BOTH orientations -> active inversion in vivo
  - each time point locked to its own assembly's orientation only -> consistent
    with a fixed strain difference / strain turnover (cannot call phase variation)

HEAVY (alignments) - run detached / SLURM. No subprocess timeouts.
Output: /home/shuaiw/borg/revision/phase_variation/targeted_inversion_reads.csv
"""

import csv
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

import orientation_driver as od  # reuse align/count_crossing/build_inverted/paths

RUN2 = od.RUN2
BASE = Path("/home/shuaiw/borg/revision/phase_variation")
SV = BASE / "rm_sv_assembly.tsv"
STRAINS = BASE / "longitudinal_strains.tsv"
OUT = BASE / "targeted_inversion"
PAD = 20


def sample_of(contig):
    m = re.match(r"^(.+)_\d+_[CL]$", contig)
    return m.group(1) if m else None


def find_segment(ref_fa, contig, reg_start, reg_end, query_fa):
    """Return (seg_start, seg_end) in contig coords = largest minus-strand block."""
    with tempfile.TemporaryDirectory() as td:
        region = Path(td) / "reg.fa"
        with open(region, "w") as o:
            subprocess.run(["seqkit", "subseq", "-r", f"{reg_start}:{reg_end}", str(ref_fa)],
                           check=True, stdout=o, stderr=subprocess.DEVNULL)
        cmd = ["blastn", "-query", str(region), "-subject", str(query_fa),
               "-outfmt", "6 qstart qend sstrand length", "-max_hsps", "50",
               "-evalue", "1e-10"]
        r = subprocess.run(cmd, capture_output=True, text=True)
    minus = []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 4 or f[2] != "minus":
            continue
        qs, qe = int(f[0]), int(f[1])
        minus.append((min(qs, qe), max(qs, qe)))
    if not minus:
        return None
    minus.sort()
    merged = [list(minus[0])]
    for s, e in minus[1:]:
        if s <= merged[-1][1] + 200:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    seg = max(merged, key=lambda x: x[1] - x[0])
    # region pos 1 == contig reg_start
    return (reg_start + seg[0] - 1, reg_start + seg[1] - 1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sv_rows = [r for r in csv.DictReader(open(SV), delimiter="\t")
               if r["sv_class"] == "inversion"]
    # unique locus = (cluster, operon, ref_contig, ref_start, ref_end); keep one query
    loci = {}
    for r in sv_rows:
        k = (r["cluster"], r["operon"], r["ref_contig"], r["ref_start"], r["ref_end"])
        loci.setdefault(k, r)

    # all timepoints/members per strain
    members = defaultdict(list)
    for st in csv.DictReader(open(STRAINS), delimiter="\t"):
        for m in st["members"].split(";"):
            tp, contig = m.split(":", 1)
            members[st["cluster"]].append((int(tp), contig))

    results = []
    for (cluster, operon, ref_contig, rs, re_), r in loci.items():
        rs, re_ = int(rs), int(re_)
        species = r["species"]
        ref_sample = sample_of(ref_contig)
        ref_fa = od.contig_fa(ref_sample, ref_contig)
        q_contig = r["query_contig"].split(";")[0]
        q_fa = od.contig_fa(sample_of(q_contig), q_contig)
        if not (ref_fa and ref_fa.exists() and q_fa and q_fa.exists()):
            print(f"[skip] {cluster} {operon}: missing fasta")
            continue
        seg = find_segment(ref_fa, ref_contig, rs, re_, q_fa)
        if seg is None:
            print(f"[skip] {cluster} {operon}: no minus block")
            continue
        seg_s, seg_e = seg
        tag = f"{cluster}_{operon.replace(' ', '').replace('#','')}"
        print(f"\n=== {cluster} {species} {operon}: ref {ref_contig} "
              f"inverted segment {seg_s}-{seg_e} ({seg_e-seg_s}bp) ===", flush=True)
        rdir = OUT / tag
        rdir.mkdir(parents=True, exist_ok=True)
        raw_fa, inv_fa = rdir / "raw.fa", rdir / "inv.fa"
        od.build_inverted(ref_fa, raw_fa, inv_fa, ref_contig, seg_s, seg_e)
        for fa in (raw_fa, inv_fa):
            if not Path(str(fa) + ".fai").exists():
                subprocess.run(["samtools", "faidx", str(fa)], check=True)
        for tp, contig in sorted(set(members[cluster])):
            samp = sample_of(contig)
            rb = od.filtered_bam(samp, contig)
            if rb is None:
                continue
            raw_bam = rdir / f"tp{tp}.raw.bam"
            inv_bam = rdir / f"tp{tp}.inv.bam"
            od.align(rb, raw_fa, raw_bam)
            od.align(rb, inv_fa, inv_bam)
            rc = od.count_crossing(raw_bam, ref_contig, seg_s, seg_e, PAD)
            ic = od.count_crossing(inv_bam, ref_contig, seg_s, seg_e, PAD)
            tot = rc + ic
            frac = ic / tot if tot else float("nan")
            mixed = "yes" if (min(rc, ic) >= 3 and min(rc, ic) / tot >= 0.1) else "no"
            print(f"  tp{tp} {samp}: raw={rc} inv={ic} inv_frac={frac:.3f} mixed={mixed}",
                  flush=True)
            results.append(dict(cluster=cluster, species=species, operon=operon,
                                ref_contig=ref_contig, seg_start=seg_s, seg_end=seg_e,
                                timepoint=tp, sample=samp, raw_crossing=rc,
                                inv_crossing=ic, inv_fraction=frac, mixed_in_tp=mixed))

    outcsv = BASE / "targeted_inversion_reads.csv"
    with open(outcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(results[0].keys()) if results else ["cluster"])
        w.writeheader()
        for row in results:
            w.writerow(row)
    print(f"\nWrote {outcsv} ({len(results)} rows)")


if __name__ == "__main__":
    main()
