#!/usr/bin/env python3
"""Data builder for the metagenome combined 3-panel figure, over the strict metagenome ECE set
(filterpass_FINAL.csv, 3,880 pass ECEs across 64 run2 metagenomes). Writes three Source-Data CSVs
read by plot_fig_metagenome.R. Read-only on all pipeline outputs.

Panels:
  a  length distribution of linked vs unlinked strict ECEs (all 3,880)
  b  paired scatter, linked ECEs only: host motif-set density on the linked host (x) vs on the ECE (y)
  c  density of the whole metagenome motif set (all.motifs.csv, "any motif pass filter") on
     linked host contigs vs linked ECEs vs unlinked ECEs

Definitions: confident link = final_score > 0.5 & specificity < 0.01 (from host_summary.csv).
Occurrence density = sum of profile motif_loci_num over a motif set / (contig_len/1000).
"""
import os
import csv
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

import ece_plot_common as C
from ece_linkability import (read_profile, filtered_set, set_density,
                             read_host_summary, gff_mod_count)

STRICT = "/home/shuaiw/borg/revision/ece_anno/expanded/filterpass_FINAL.csv"
# authoritative high-confidence linkage set (317 strict ECEs); linked = MGE present here
LINKED_TABLE = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/mge_host_gc_cov_final.csv"
RUN = "/home/shuaiw/borg/paper/run2"
VARIANTS = ["_methylation3", "_methylation4", "_methylation2", "_methylation"]
STRICT_MAP = None   # {sample: {MGE: (type, mge_len)}}, set in main (inherited by fork)
LINKED_MGE = None   # set of high-confidence linked MGE names (curated 317), set in main
LOOSE = False       # if True, linked = host_summary final_score>0.5 & specificity<0.01 (no curation)
S0, P0 = 0.5, 0.01  # loose confident-link thresholds


def load_strict(path):
    df = pd.read_csv(path)
    df = df[df["pass"] == 1]
    m = {}
    for _, r in df.iterrows():
        try:
            ln = int(float(r["mge_len"]))
        except (ValueError, TypeError):
            ln = None
        m.setdefault(str(r["sample"]), {})[str(r["MGE"])] = (str(r["MGE_type"]), ln)
    return m


def read_fai(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                try:
                    out[f[0]] = int(f[1])
                except ValueError:
                    pass
    return out


def read_host_ctgs(path):
    out = set()
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            c = line.split("\t")[0].strip()
            if c:
                out.add(c)
    return out


def read_meta_motifset(path):
    """all.motifs.csv -> set of (motifString, centerPos) = the sample's detected motif repertoire."""
    s = set()
    if not os.path.exists(path):
        return s
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                s.add((row["motifString"], int(row["centerPos"])))
            except (KeyError, ValueError, TypeError):
                continue
    return s


def resolve_wd(sample):
    for v in VARIANTS:
        p = os.path.join(RUN, sample, sample + v)
        if os.path.exists(os.path.join(p, "host_summary.csv")):
            return p
    return None


def worker(sample):
    try:
        wd = resolve_wd(sample)
        if wd is None:
            return []
        base = os.path.join(RUN, sample)
        fai = read_fai(os.path.join(base, f"{sample}.hifiasm.p_ctg.rename.fa.fai"))
        hs = read_host_summary(os.path.join(wd, "host_summary.csv"))   # {MGE:(host,fs,sp,len)}
        meta_set = read_meta_motifset(os.path.join(wd, "all.motifs.csv"))
        strict = STRICT_MAP.get(sample, {})

        pcache = {}

        def prof(ctg):
            if ctg not in pcache:
                pcache[ctg] = read_profile(os.path.join(wd, "profiles", f"{ctg}.motifs.profile.csv"))
            return pcache[ctg]

        def clen(ctg, fallback=None):
            return fai.get(ctg, fallback)

        a_rows, b_rows, c_rows = [], [], []
        linked_hosts = set()

        for ece, (typ, mlen) in strict.items():
            row = hs.get(ece)
            host = row[0] if row else ""
            if LOOSE:   # loose: any host_summary link passing final_score>0.5 & specificity<0.01
                fs = row[1] if row else np.nan
                sp = row[2] if row else np.nan
                linked = bool(row) and fs == fs and sp == sp and fs > S0 and sp < P0
            else:       # curated high-confidence linkage set (317)
                linked = ece in LINKED_MGE
            ece_len = clen(ece, mlen)

            # modification sites/frequency from the ECE's GFF (kinModCall, score >= 30)
            nmod = gff_mod_count(os.path.join(wd, "gffs", f"{ece}.gff"))
            mod_density = (nmod / (ece_len / 1000.0)) if (nmod is not None and ece_len) else np.nan
            a_rows.append(dict(sample=sample, ece=ece, type=typ, mge_len=mlen, linked=linked,
                               n_mod_sites=nmod if nmod is not None else np.nan,
                               mod_density_per_kb=mod_density))

            # panel c: metagenome-set occurrences + density on every ECE
            if meta_set and ece_len:
                n, d = set_density(prof(ece), meta_set, ece_len)
                c_rows.append(dict(sample=sample, ece=ece, ece_len=ece_len,
                                   group="linked_ECE" if linked else "unlinked_ECE",
                                   density_per_kb=d, n_occ=n))

            if linked and host:
                linked_hosts.add(host)
                # panel b: host motif set (filtered detected motifs on the linked host contig)
                hset = filtered_set(prof(host))
                host_len = clen(host)
                if hset and host_len and ece_len:
                    _, hd = set_density(prof(host), hset, host_len)
                    _, ed = set_density(prof(ece), hset, ece_len)
                    b_rows.append(dict(sample=sample, ece=ece, type=typ,
                                       host_density=hd, ece_density=ed))

        # panel c: metagenome-set density on each (deduped) linked host contig
        for hctg in linked_hosts:
            hl = clen(hctg)
            if meta_set and hl:
                n, d = set_density(prof(hctg), meta_set, hl)
                c_rows.append(dict(sample=sample, ece=hctg, ece_len=hl,
                                   group="linked_host", density_per_kb=d, n_occ=n))

        return [("a", r) for r in a_rows] + [("b", r) for r in b_rows] + [("c", r) for r in c_rows]
    except Exception as e:
        return [("err", dict(sample=sample, error=str(e)))]


def main():
    global STRICT_MAP, LINKED_MGE, LOOSE
    ap = argparse.ArgumentParser()
    ap.add_argument("--loose", action="store_true",
                    help="linked = host_summary final_score>0.5 & specificity<0.01 (no curation)")
    args = ap.parse_args()
    LOOSE = args.loose
    prefix = "loose_" if LOOSE else ""

    STRICT_MAP = load_strict(STRICT)
    strict_mge = {e for v in STRICT_MAP.values() for e in v}
    lk = pd.read_csv(LINKED_TABLE)
    LINKED_MGE = set(lk.MGE.astype(str)) & strict_mge
    samples = sorted(STRICT_MAP)
    n = sum(len(v) for v in STRICT_MAP.values())
    mode = "LOOSE (host_summary score>0.5 & spec<0.01)" if LOOSE else \
           f"curated high-confidence ({len(LINKED_MGE)})"
    print(f"strict metagenome ECEs: {n} across {len(samples)} samples; linked mode: {mode}")

    A, B, Cc, errs = [], [], [], []
    with Pool(64) as pool:
        for res in pool.imap_unordered(worker, samples, chunksize=1):
            for tag, r in res:
                (A if tag == "a" else B if tag == "b" else Cc if tag == "c" else errs).append(r)
    if errs:
        print(f"WARNING {len(errs)} sample errors, e.g. {errs[0]}")

    dfa, dfb, dfc = pd.DataFrame(A), pd.DataFrame(B), pd.DataFrame(Cc)
    os.makedirs(C.OUT, exist_ok=True)
    dfa.to_csv(os.path.join(C.OUT, f"fig_metagenome_{prefix}a_length_sourcedata.csv"), index=False)
    dfb.to_csv(os.path.join(C.OUT, f"fig_metagenome_{prefix}b_scatter_sourcedata.csv"), index=False)
    dfc.to_csv(os.path.join(C.OUT, f"fig_metagenome_{prefix}c_density_sourcedata.csv"), index=False)

    print(f"a: {len(dfa)} ECEs  linked={int(dfa.linked.sum())} unlinked={int((~dfa.linked).sum())} "
          f"types={dfa.type.value_counts().to_dict()}")
    print(f"b: {len(dfb)} linked ECEs with host motif set")
    print(f"c groups: {dfc.group.value_counts().to_dict()}")
    if len(dfb):
        from scipy import stats
        p = stats.wilcoxon(dfb.host_density, dfb.ece_density).pvalue
        below = float(np.mean(dfb.ece_density < dfb.host_density)) * 100
        print(f"   panel b: {below:.0f}% below y=x, paired Wilcoxon p={p:.2e}")
    if len(dfc):
        med = dfc.groupby("group").density_per_kb.median()
        print(f"   panel c medians/kb: {med.to_dict()}")


if __name__ == "__main__":
    main()
