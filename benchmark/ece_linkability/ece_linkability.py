#!/usr/bin/env python3
"""
ece_linkability.py -- assemble a per-ECE master table for the linkage-ability
reviewer analysis, across three datasets:

  1. isolate   : the paper's 1,420 pure high-depth isolates (batch2_results, _methylation4)
  2. simulated : the bg300 complexity community, 5 reps
  3. metagenome: the ~60 real metagenomes (run2)

For every extrachromosomal element (ECE = plasmid / virus / novel) it records length,
motif density per kb, the linked/unlinked call, and (for the ground-truth sets) whether
the confident link is to the correct host.

Definitions (fixed with the user):
  * Evaluable contig gate (isolate motif-sharing criteria, cal_jaccard):
        depth >= 10x AND length >= 5000 bp.
    We record both `passes_depth_gate` (depth>=10 only, for the length panel) and
    `passes_full_gate` (depth+length) so each figure can pick its universe.
  * Confident link = best host has final_score > 0.5 AND specificity < 0.01
    (recomputed from host_summary.csv directly -- host_summary.filter.csv can be stale/empty).
  * linked:
        isolate + simulated -> confident AND to the CORRECT true host (true positive)
        metagenome          -> any confident link (no ground truth)
  * Motif density = sequence occurrences per kb (motif_loci_num, fwd+rev) of a motif set,
    where the motif set is "filtered by MODIFI default": on the defining contig the motif has
    motif_modified_ratio > 0.4 (min_frac) AND motif_modified_num > 30 (min_sites).
        isolate + simulated -> count the HOST genome's motif set within the ECE
        metagenome          -> count the ECE's OWN motif set within the ECE

Read-only on all pipeline outputs. Writes only:
    /home/shuaiw/borg/revision/ece_linkability/ece_linkability_master.csv
Run on the terminal (I/O bound), multiprocessing capped at 64 workers.
"""
import os
import csv
import glob
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

# ---- thresholds (MODIFI defaults / isolate motif-sharing criteria) ----
MIN_FRAC = 0.4      # main.py --min_frac default; host/ECE motif must exceed this modified ratio
MIN_SITES = 30      # main.py --min_sites default; and exceed this number of methylated sites
S0, P0 = 0.5, 0.01  # filter_linkage confident-link thresholds
DEPTH_GATE = 10     # sample_object.depth_cutoff
LEN_GATE = 5000     # sample_object.length_cutoff

# ---- dataset roots ----
ISO_SUMMARY = "/home/shuaiw/borg/paper/isolation/GTDB_tree/anno/isolation_sample_summary.tsv"
ISO_BASE = "/home/shuaiw/borg/paper/isolation/batch2_results"
ISO_VARIANT = "_methylation4"                         # paper's Isolation_sample uses _methylation4
SIM_ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
SIM_REPS = ["bg_300", "bg_300_rep2", "bg_300_rep3", "bg_300_rep4", "bg_300_rep5"]
META_ROOT = "/home/shuaiw/borg/paper/run2"
META_VARIANTS = ["_methylation3", "_methylation4", "_methylation2", "_methylation"]

OUT_CSV = "/home/shuaiw/borg/revision/ece_linkability/ece_linkability_master.csv"

sra = lambda x: str(x).split("_")[0]

# optional strict isolate ECE whitelist {sample: {ece: (type, len, depth)}}, set from --iso_strict
ISO_STRICT_MAP = None


def load_iso_strict(path):
    """filterpass_isolate CSV -> {sample: {ece: (type, len, depth)}} for pass==1 rows."""
    df = pd.read_csv(path)
    df = df[df["pass"] == 1]
    m = {}
    for _, r in df.iterrows():
        s, e = str(r["sample"]), str(r["MGE"])
        try:
            ln = int(float(r["mge_len"]))
        except (ValueError, TypeError):
            ln = None
        try:
            dep = float(r["mge_depth"])
        except (ValueError, TypeError):
            dep = np.nan
        m.setdefault(s, {})[e] = (str(r["MGE_type"]), ln, dep)
    return m


# ----------------------------------------------------------------------
# small readers
# ----------------------------------------------------------------------
def read_profile(path):
    """Read a per-contig profile -> dict keyed by (motifString, centerPos) with
    (motif_loci_num, motif_modified_num, motif_modified_ratio). Returns {} if absent."""
    if not os.path.exists(path):
        return {}
    d = {}
    try:
        with open(path) as fh:
            r = csv.DictReader(fh)
            for row in r:
                try:
                    key = (row["motifString"], int(row["centerPos"]))
                    d[key] = (
                        int(float(row["motif_loci_num"])),
                        int(float(row["motif_modified_num"])),
                        float(row["motif_modified_ratio"]),
                    )
                except (KeyError, ValueError, TypeError):
                    continue
    except OSError:
        return {}
    return d


def filtered_set(profile):
    """MODIFI-default filtered motif set of a contig from its profile dict."""
    return {k for k, (loci, mod_n, mod_r) in profile.items()
            if mod_r > MIN_FRAC and mod_n > MIN_SITES}


def set_density(target_profile, motif_set, length_bp):
    """Occurrences (motif_loci_num) of motif_set within a target contig, and per-kb density."""
    if not motif_set or not target_profile or not length_bp or length_bp <= 0:
        return 0, np.nan
    n = 0
    for k in motif_set:
        v = target_profile.get(k)
        if v is not None:
            n += v[0]
    return n, n / (length_bp / 1000.0)


def gff_mod_count(path, min_score=30):
    """Total modified positions on a contig from its per-contig kinModCall GFF, counting
    `modified_base` lines whose score (column 6) is >= min_score. Motif-independent.
    Returns None if the GFF is absent."""
    if not os.path.exists(path):
        return None
    n = 0
    try:
        with open(path) as fh:
            for line in fh:
                if not line or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 6 or f[2] != "modified_base":
                    continue
                try:
                    if float(f[5]) >= min_score:
                        n += 1
                except ValueError:
                    continue
    except OSError:
        return None
    return n


def read_depth(path):
    """mean_depth.csv -> {contig: (depth, length_or_None)}."""
    out = {}
    if not os.path.exists(path):
        return out
    try:
        with open(path) as fh:
            r = csv.DictReader(fh)
            for row in r:
                c = row.get("contig")
                if c is None:
                    continue
                try:
                    dep = float(row["depth"])
                except (KeyError, ValueError, TypeError):
                    dep = np.nan
                ln = None
                if row.get("length") not in (None, ""):
                    try:
                        ln = int(float(row["length"]))
                    except ValueError:
                        ln = None
                out[c] = (dep, ln)
    except OSError:
        return out
    return out


def read_all_mge(path):
    """all_mge.tsv / mge_list.tsv -> {seq_name: (type, length)}."""
    out = {}
    if not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return out
    name_col = "seq_name" if "seq_name" in df.columns else df.columns[0]
    type_col = "type" if "type" in df.columns else ("mge_type" if "mge_type" in df.columns else None)
    len_col = "length" if "length" in df.columns else None
    for _, row in df.iterrows():
        nm = str(row[name_col])
        if nm in ("seq_name", "nan"):
            continue
        typ = str(row[type_col]) if type_col else "NA"
        try:
            ln = int(float(row[len_col])) if len_col and not pd.isna(row[len_col]) else None
        except (ValueError, TypeError):
            ln = None
        out[nm] = (typ, ln)
    return out


def read_host_candidates(path):
    """hosts/<ece>.host_prediction.csv -> list of (host, final_score, specificity) for every
    candidate host of one ECE (the full ranked list, not just the best). [] if absent."""
    if not os.path.exists(path):
        return []
    out = []
    try:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                try:
                    fs = float(row["final_score"])
                except (KeyError, ValueError, TypeError):
                    fs = np.nan
                try:
                    sp = float(row["specificity"])
                except (KeyError, ValueError, TypeError):
                    sp = np.nan
                out.append((row.get("host", ""), fs, sp))
    except OSError:
        return []
    return out


def read_host_summary(path):
    """host_summary.csv -> {MGE: (host, final_score, specificity, MGE_len)} keeping best score."""
    out = {}
    if not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path)
    except Exception:
        return out
    if "MGE" not in df.columns:
        return out
    for _, row in df.iterrows():
        mge = str(row["MGE"])
        try:
            fs = float(row["final_score"])
        except (ValueError, TypeError):
            fs = np.nan
        try:
            sp = float(row["specificity"])
        except (ValueError, TypeError):
            sp = np.nan
        ln = None
        if "MGE_len" in df.columns and not pd.isna(row["MGE_len"]):
            try:
                ln = int(float(row["MGE_len"]))
            except (ValueError, TypeError):
                ln = None
        prev = out.get(mge)
        if prev is None or (fs == fs and (prev[1] != prev[1] or fs > prev[1])):
            out[mge] = (str(row["host"]), fs, sp, ln)
    return out


# ----------------------------------------------------------------------
# generic per-sample assembler
# ----------------------------------------------------------------------
def build_records(dataset, sample, work_dir, mge_path, profiles_dir,
                  ground_truth, rep=None, ece_override=None, override_depth=None):
    """Return a list of per-ECE dict records for one sample/run.

    ground_truth: True for isolate/simulated (host known -> host motif set, correctness),
                  False for metagenome (ECE-self motif set, any confident link).
    ece_override: if given, {ece: (type, len)} to use as the ECE universe instead of all_mge
                  (still uses all_mge names to identify host contigs). override_depth: {ece: depth}
                  fallback when mean_depth.csv lacks the contig."""
    mge = read_all_mge(mge_path)
    depth = read_depth(os.path.join(work_dir, "mean_depth.csv"))
    hs = read_host_summary(os.path.join(work_dir, "host_summary.csv"))
    mge_names = set(mge.keys())
    if ece_override is not None:
        mge_names |= set(ece_override.keys())
        universe = ece_override
    else:
        universe = mge
    if not universe:
        return []

    # profile cache within this sample
    pcache = {}

    def prof(ctg):
        if ctg not in pcache:
            pcache[ctg] = read_profile(os.path.join(profiles_dir, f"{ctg}.motifs.profile.csv"))
        return pcache[ctg]

    # ---- host genome motif set(s) for ground-truth datasets ----
    # host contigs = non-MGE contigs with depth >= gate. Grouped by SRA prefix for the
    # simulated community; a single group ("__ISO__") for a pure isolate.
    gffs_dir = os.path.join(work_dir, "gffs")
    host_set_by_group = {}
    host_len_by_group = {}
    host_loci_by_group = {}
    if ground_truth:
        groups = {}
        for ctg, (dep, ln) in depth.items():
            if ctg in mge_names:
                continue
            if dep is None or dep != dep or dep < DEPTH_GATE:
                continue
            g = sra(ctg) if dataset == "simulated" else "__ISO__"
            groups.setdefault(g, []).append((ctg, ln))
        for g, ctgs in groups.items():
            uset = set()
            tot_len = 0
            per_contig = []
            for ctg, ln in ctgs:
                p = prof(ctg)
                uset |= filtered_set(p)
                per_contig.append((ctg, p, ln))
            # union set fixed; now host genome-wide occurrences and length
            tot_loci = 0
            for ctg, p, ln in per_contig:
                length = ln if ln else None
                nloci, _ = set_density(p, uset, length if length else 1)
                tot_loci += nloci
                if length:
                    tot_len += length
            host_set_by_group[g] = uset
            host_len_by_group[g] = tot_len
            host_loci_by_group[g] = tot_loci

    records = []
    for ece, (typ, mlen) in universe.items():
        dep, dlen = depth.get(ece, (np.nan, None))
        if (dep != dep) and override_depth:
            dep = override_depth.get(ece, np.nan)
        ece_len = mlen or dlen
        row = hs.get(ece)
        assigned_host = row[0] if row else ""
        final_score = row[1] if row else np.nan
        specificity = row[2] if row else np.nan
        if ece_len is None and row and row[3]:
            ece_len = row[3]

        confident = bool(row) and (final_score == final_score) and \
            (final_score > S0) and (specificity == specificity) and (specificity < P0)

        # correctness for ground-truth sets
        if ground_truth:
            if dataset == "simulated":
                correct = bool(assigned_host) and sra(assigned_host) == sra(ece)
            else:  # isolate: single strain -> correct if a confident candidate host is a
                   # chromosomal (non-MGE) contig. Scan the full ranked candidate list, not just
                   # the top host: the true host still counts even if it is not ranked #1.
                correct = bool(assigned_host) and (assigned_host not in mge_names) and confident
                if not correct:
                    for h, fs, sp in read_host_candidates(
                            os.path.join(work_dir, "hosts", f"{ece}.host_prediction.csv")):
                        if fs == fs and sp == sp and fs > S0 and sp < P0 \
                                and h and h not in mge_names:
                            correct = True
                            break
            linked = correct if dataset != "simulated" else (confident and correct)
            true_host_present = (sra(ece) in host_set_by_group) if dataset == "simulated" \
                else ("__ISO__" in host_set_by_group)
        else:
            correct = np.nan
            linked = confident
            true_host_present = np.nan

        # recognition-site density (host or ECE-self motif set)
        ece_prof = prof(ece)
        if ground_truth:
            motif_set_source = "host"
            g = sra(ece) if dataset == "simulated" else "__ISO__"
            hset = host_set_by_group.get(g, set())
            n_sites, dens = set_density(ece_prof, hset, ece_len)
            hlen = host_len_by_group.get(g, 0)
            host_own_density = (host_loci_by_group.get(g, 0) / (hlen / 1000.0)) if hlen else np.nan
        else:
            motif_set_source = "ece_self"
            eset = filtered_set(ece_prof)
            n_sites, dens = set_density(ece_prof, eset, ece_len)
            host_own_density = np.nan

        passes_depth = (dep == dep) and (dep >= DEPTH_GATE)
        passes_full = passes_depth and (ece_len is not None) and (ece_len >= LEN_GATE)

        # modification sites/density from the ECE's own GFF (motif-independent, score >= 30);
        # only read the GFF for depth-passing ECEs to bound I/O
        n_mod, mod_dens = np.nan, np.nan
        if passes_depth:
            nmod = gff_mod_count(os.path.join(gffs_dir, f"{ece}.gff"))
            if nmod is not None:
                n_mod = nmod
                mod_dens = nmod / (ece_len / 1000.0) if ece_len else np.nan

        records.append(dict(
            dataset=dataset, sample=sample, rep=rep if rep else sample,
            ece_id=ece, type=typ, ece_len=ece_len, depth=dep,
            passes_depth_gate=passes_depth, passes_full_gate=passes_full,
            motif_set_source=motif_set_source, n_motif_sites=n_sites,
            motif_density_per_kb=dens, host_own_density_per_kb=host_own_density,
            n_mod_sites=n_mod, mod_density_per_kb=mod_dens,
            assigned_host=assigned_host, final_score=final_score, specificity=specificity,
            confident_link=confident, correct_host=correct, linked=linked,
            true_host_present=true_host_present,
        ))
    return records


# ----------------------------------------------------------------------
# per-dataset workers (one call per sample -> list of records)
# ----------------------------------------------------------------------
def worker_isolate(sample):
    work_dir = os.path.join(ISO_BASE, sample, sample + ISO_VARIANT)
    mge_path = os.path.join(ISO_BASE, sample, "all_mge.tsv")
    profiles_dir = os.path.join(work_dir, "profiles")
    override = odepth = None
    if ISO_STRICT_MAP is not None:
        entry = ISO_STRICT_MAP.get(sample)
        if not entry:
            return []
        override = {e: (t, l) for e, (t, l, d) in entry.items()}
        odepth = {e: d for e, (t, l, d) in entry.items()}
    try:
        return build_records("isolate", sample, work_dir, mge_path, profiles_dir,
                             ground_truth=True, ece_override=override, override_depth=odepth)
    except Exception as e:  # never let one sample kill the pool
        return [dict(dataset="isolate", sample=sample, ece_id="__ERROR__", error=str(e))]


def worker_simulated(rep):
    work_dir = os.path.join(SIM_ROOT, rep, "modifi", rep)
    mge_path = os.path.join(SIM_ROOT, rep, f"{rep}.mge_list.tsv")
    profiles_dir = os.path.join(work_dir, "profiles")
    try:
        return build_records("simulated", rep, work_dir, mge_path, profiles_dir,
                             ground_truth=True, rep=rep)
    except Exception as e:
        return [dict(dataset="simulated", sample=rep, ece_id="__ERROR__", error=str(e))]


def worker_metagenome(sample):
    base = os.path.join(META_ROOT, sample)
    work_dir = None
    for v in META_VARIANTS:
        if os.path.exists(os.path.join(base, sample + v, "host_summary.csv")):
            work_dir = os.path.join(base, sample + v)
            break
    if work_dir is None:
        return []
    mge_path = os.path.join(base, "all_mge.tsv")
    profiles_dir = os.path.join(work_dir, "profiles")
    try:
        return build_records("metagenome", sample, work_dir, mge_path, profiles_dir,
                             ground_truth=False)
    except Exception as e:
        return [dict(dataset="metagenome", sample=sample, ece_id="__ERROR__", error=str(e))]


# ----------------------------------------------------------------------
def isolate_roster():
    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return -1.0
    with open(ISO_SUMMARY) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return [r["Sample"] for r in rows
            if r.get("Pure_anno") == "pure" and num(r.get("Average_DP")) >= 10]


def metagenome_samples():
    out = []
    for d in sorted(glob.glob(os.path.join(META_ROOT, "*"))):
        s = os.path.basename(d)
        for v in META_VARIANTS:
            if os.path.exists(os.path.join(d, s + v, "host_summary.csv")):
                out.append(s)
                break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--datasets", default="isolate,simulated,metagenome")
    ap.add_argument("--out", default=OUT_CSV)
    ap.add_argument("--iso_strict", default=None,
                    help="restrict isolate ECEs to pass==1 rows of this filterpass CSV")
    args = ap.parse_args()
    threads = min(args.threads, 64)
    want = set(args.datasets.split(","))

    global ISO_STRICT_MAP
    if args.iso_strict:
        ISO_STRICT_MAP = load_iso_strict(args.iso_strict)
        print(f"[isolate] strict ECE list: {sum(len(v) for v in ISO_STRICT_MAP.values())} "
              f"pass ECEs across {len(ISO_STRICT_MAP)} samples")

    all_recs = []

    if "isolate" in want:
        roster = sorted(ISO_STRICT_MAP) if ISO_STRICT_MAP else isolate_roster()
        print(f"[isolate] {len(roster)} isolates")
        with Pool(threads) as pool:
            for recs in pool.imap_unordered(worker_isolate, roster, chunksize=8):
                all_recs.extend(recs)
        print(f"[isolate] done, {len(all_recs)} ECE rows so far")

    if "simulated" in want:
        print(f"[simulated] {len(SIM_REPS)} bg300 reps")
        with Pool(min(threads, len(SIM_REPS))) as pool:
            for recs in pool.imap_unordered(worker_simulated, SIM_REPS):
                all_recs.extend(recs)
        print(f"[simulated] done, {len(all_recs)} ECE rows cumulative")

    if "metagenome" in want:
        samples = metagenome_samples()
        print(f"[metagenome] {len(samples)} samples with host_summary")
        with Pool(min(threads, len(samples) or 1)) as pool:
            for recs in pool.imap_unordered(worker_metagenome, samples, chunksize=1):
                all_recs.extend(recs)
        print(f"[metagenome] done, {len(all_recs)} ECE rows cumulative")

    df = pd.DataFrame(all_recs)
    errors = df[df.get("ece_id") == "__ERROR__"] if "ece_id" in df.columns else pd.DataFrame()
    if len(errors):
        print(f"WARNING: {len(errors)} sample-level errors; see error column")
        print(errors[["dataset", "sample", "error"]].head(20).to_string(index=False))
        df = df[df["ece_id"] != "__ERROR__"].copy()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\nwrote {args.out}  ({len(df)} ECE rows)")

    # quick summary
    for ds in ["isolate", "simulated", "metagenome"]:
        sub = df[df.dataset == ds]
        if not len(sub):
            continue
        pv = sub[sub.type.isin(["plasmid", "virus"])]
        full = sub[sub.passes_full_gate]
        print(f"[{ds}] ECEs={len(sub)}  plasmid/virus={len(pv)}  "
              f"full-gate={len(full)}  linked(full-gate)={int(full.linked.sum())}")


if __name__ == "__main__":
    main()
