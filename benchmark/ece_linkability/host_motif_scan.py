#!/usr/bin/env python3
"""Sequence-scan builder for the paired host-vs-ECE motif-density figure, over the FULL strict
isolate ECE set (all 1,254 pass ECEs, no depth/profile filter).

Rationale: the ECE value is NOT a methylation detection; it is a pure sequence scan for the host's
recognition motifs. So it needs no read depth and no per-ECE modification profile - only the ECE
DNA sequence. We therefore count host-motif occurrences (IUPAC-aware, both strands) directly in
each ECE (and host chromosome) sequence.

Per sample: host motif set = the host chromosome's default-filtered detected motifs
(modified_ratio > 0.4 & modified_num > 30, from the high-depth host contig profiles). x = host-set
occurrences/kb on the host chromosome; y = host-set occurrences/kb on each ECE.

Writes fig_host_motif_density_sourcedata.csv (type, host_density, ece_density) for plot_host_motif_density.R.
"""
import os
import csv
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import nt_search
from Bio.Seq import Seq

import ece_plot_common as C
from ece_linkability import (read_profile, filtered_set, read_all_mge,
                             ISO_BASE, ISO_VARIANT, load_iso_strict)

STRICT = "/home/shuaiw/borg/revision/ece_anno/isolate_all/filterpass_isolate_FINAL.csv"
STEM = "fig_host_motif_density"
STRICT_MAP = None   # {sample: {ece: (type, len, depth)}}, set in main (inherited by fork)


def count_occ(seq, motif):
    """Occurrences of a IUPAC motif on both strands (matches profile motif_loci_num semantics:
    for_loci + rev_loci, no palindrome dedup)."""
    n = len(nt_search(seq, motif)) - 1
    rc = str(Seq(motif).reverse_complement())
    n += len(nt_search(seq, rc)) - 1
    return n


def string_loci(profile):
    """{motifString: motif_loci_num} from a profile dict (occurrences are per-string, so collapse
    over centerPos by taking the max)."""
    d = {}
    for (ms, cp), (loci, mn, mr) in profile.items():
        d[ms] = max(d.get(ms, 0), loci)
    return d


def worker(sample):
    try:
        base = os.path.join(ISO_BASE, sample)
        fa = os.path.join(base, f"{sample}.hifiasm.p_ctg.rename.fa")
        wd = os.path.join(base, sample + ISO_VARIANT)
        if not os.path.exists(fa):
            return []
        seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(fa, "fasta")}
        mge_names = set(read_all_mge(os.path.join(base, "all_mge.tsv")).keys())
        strict = STRICT_MAP.get(sample, {})
        mge_names |= set(strict.keys())

        # host motif set = union of default-filtered motifs over host (non-MGE) contigs.
        # Host occurrences come from the pipeline profiles (motif_loci_num), no re-scan needed.
        host_ctgs = [c for c in seqs if c not in mge_names]
        hset = set()
        host_loci_by_ctg = {}
        for c in host_ctgs:
            prof = read_profile(os.path.join(wd, "profiles", f"{c}.motifs.profile.csv"))
            hset |= filtered_set(prof)
            host_loci_by_ctg[c] = string_loci(prof)
        motifs = sorted({ms for ms, cp in hset})   # unique recognition strings
        if not motifs:
            return []

        # host chromosome density of the whole host motif set (from profile motif_loci_num)
        host_len = sum(len(seqs[c]) for c in host_ctgs)
        if not host_len:
            return []
        host_occ = sum(host_loci_by_ctg[c].get(m, 0) for c in host_ctgs for m in motifs)
        host_density = host_occ / (host_len / 1000.0)

        # one point per host-ECE link: x = host-set density on host, y = host-set density on ECE
        rows = []
        for ece, (typ, ln, dep) in strict.items():
            if ece not in seqs:
                continue
            s = seqs[ece]
            L = len(s)
            if not L:
                continue
            ece_occ = sum(count_occ(s, m) for m in motifs)
            rows.append(dict(sample=sample, ece=ece, type=typ,
                             host_density=host_density, ece_density=ece_occ / (L / 1000.0)))
        return rows
    except Exception as e:
        return [dict(sample=sample, ece="__ERROR__", type=str(e),
                     host_density=np.nan, ece_density=np.nan)]


def main():
    global STRICT_MAP
    STRICT_MAP = load_iso_strict(STRICT)
    samples = sorted(STRICT_MAP)
    print(f"scanning {sum(len(v) for v in STRICT_MAP.values())} strict ECEs across {len(samples)} isolates")
    recs = []
    with Pool(64) as pool:
        for r in pool.imap_unordered(worker, samples, chunksize=4):
            recs.extend(r)
    df = pd.DataFrame(recs)
    err = df[df.ece == "__ERROR__"] if "ece" in df.columns else pd.DataFrame()
    if len(err):
        print(f"WARNING {len(err)} sample errors, e.g. {err.type.iloc[0]}")
        df = df[df.ece != "__ERROR__"]
    df = df[df.type.isin(["plasmid", "virus"])].dropna(subset=["host_density", "ece_density"])

    os.makedirs(C.OUT, exist_ok=True)
    df[["type", "host_density", "ece_density"]].to_csv(
        os.path.join(C.OUT, f"{STEM}_sourcedata.csv"), index=False)

    from scipy import stats
    print(f"ECEs plotted: {len(df)}  {df.type.value_counts().to_dict()}")
    for t in ["plasmid", "virus"]:
        sub = df[df.type == t]
        p = stats.wilcoxon(sub.host_density, sub.ece_density).pvalue
        below = float(np.mean(sub.ece_density < sub.host_density)) * 100
        print(f"  {t}: n={len(sub)} host {sub.host_density.median():.2f} ECE {sub.ece_density.median():.2f} "
              f"below={below:.0f}% Wilcoxon p={p:.2e}")


if __name__ == "__main__":
    main()
