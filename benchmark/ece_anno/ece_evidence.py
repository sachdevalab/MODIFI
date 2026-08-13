#!/usr/bin/env python3
"""
ECE validation evidence — single-sample engine.

For every predicted ECE (rows of the linkage CSV for one sample) gather six independent
lines of evidence and a combined confidence call, reusing geNomad's already-computed
proteins/annotations so no gene re-prediction or geNomad re-run is needed.

  1. Circularity          : hifiasm "C" suffix + geNomad topology (DTR/ITR)
  2. Plasmid markers       : ConjScan (macsyfinder) + geNomad conjugation genes
  3. Viral markers         : Pfam terminase/capsid/portal (+ optional VOGdb) HMMs
  4. Chromosomal SCMG (neg): bac71/arc76 single-copy core-gene HMMs
  5. Coverage pattern      : mge_depth/host_depth ratio (CSV) + per-contig CV (BAM)
  6. geNomad robustness    : score / fdr / n_hallmarks, survives default FDR threshold

Outputs <out_dir>/ece_evidence.tsv (+ intermediates). Read-only w.r.t. the sample dir.

Usage:
  python ece_evidence.py --sample ERR10042285 \
      --work_dir /home/shuaiw/borg/paper/isolation/batch2_results/ERR10042285 \
      --ece_csv  /home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv \
      --out_dir  /home/shuaiw/borg/revision/ece_anno/per_sample/ERR10042285 \
      --db_dir   /home/shuaiw/borg/revision/ece_anno/db --threads 16
"""
import argparse
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd

import anno_utils as au
import parse_conjscan as pcj

# ----------------------------------------------------------------------------
# Tunable thresholds
# ----------------------------------------------------------------------------
GENOMAD_FDR_MAX = 0.05      # default-preset FDR cutoff for robustness re-thresholding
COV_RATIO_LOW = 0.75        # coverage-distinctness band (P4)
COV_RATIO_HIGH = 1.33
CV_MAX = 0.75               # max CV for "uniform" coverage (P4)
CV_ARTIFACT = 1.5           # CV above this contributes to artifact flag
COV_MIN = 1.0               # mean depth below this contributes to artifact flag
SCMG_MIN = 5                # >= this many chromosomal single-copy markers -> chromosomal flag
SCMG_FRAC = 0.10            # or this fraction of a marker set
VOGDB_HMM = "/shared/db/vogdb/latest/hmm/VOG.all.hmm"

CIRCULAR_TOPOLOGY = {"DTR", "ITR"}


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------
def genomad_prefix(work_dir, sample):
    """Resolve the geNomad output prefix (usually <sample>.hifiasm.p_ctg.rename)."""
    cand = f"{sample}.hifiasm.p_ctg.rename"
    if os.path.isdir(os.path.join(work_dir, "Genomad", f"{cand}_summary")):
        return cand
    hits = glob.glob(os.path.join(work_dir, "Genomad", "*_summary"))
    if hits:
        return os.path.basename(hits[0])[: -len("_summary")]
    return cand


def load_eces(ece_csv, sample):
    """ECE rows for this sample from the linkage CSV, deduped by contig."""
    df = pd.read_csv(ece_csv)
    df = df[df["prefix"].astype(str) == str(sample)].copy()
    df = df[~df["mge_contig"].astype(str).str.contains(r"\|provirus", regex=True)]
    df = df.drop_duplicates("mge_contig")
    out = pd.DataFrame({
        "seq_name": df["mge_contig"].astype(str),
        "type": df["mge_type"].astype(str),
        "length": df["mge_length"],
        "mge_depth": df["mge_depth"],
        "host_depth": df["host_depth"],
        "host_lineage": df.get("host_lineage", ""),
    })
    return out.reset_index(drop=True)


def read_genomad_summaries(work_dir, gp):
    """Per-contig geNomad summary fields for both plasmid and virus tables."""
    sdir = os.path.join(work_dir, "Genomad", f"{gp}_summary")
    rows = {}
    for kind, score_col in (("plasmid", "plasmid_score"), ("virus", "virus_score")):
        path = os.path.join(sdir, f"{gp}_{kind}_summary.tsv")
        if not os.path.exists(path):
            continue
        t = pd.read_csv(path, sep="\t")
        for _, r in t.iterrows():
            name = str(r["seq_name"]).split("|")[0]  # drop provirus coord suffix
            cj = r.get("conjugation_genes")
            cj_n = 0
            if isinstance(cj, str) and cj.strip() and cj.strip() != "NA":
                cj_n = len(cj.split(";"))
            rows[name] = {
                "genomad_topology": r.get("topology", ""),
                "genomad_score": pd.to_numeric(r.get(score_col), errors="coerce"),
                "genomad_fdr": pd.to_numeric(r.get("fdr"), errors="coerce"),
                "genomad_n_hallmarks": pd.to_numeric(r.get("n_hallmarks"), errors="coerce"),
                "genomad_conjugation_genes": cj_n,
            }
    return rows


def load_genes_tsv(work_dir, gp):
    """geNomad per-gene annotation table (annotate stage). None if cleaned away."""
    path = os.path.join(work_dir, "Genomad", f"{gp}_annotate", f"{gp}_genes.tsv")
    if not os.path.exists(path):
        return None
    g = pd.read_csv(path, sep="\t")
    g["contig"] = g["gene"].astype(str).map(au.gene_to_contig)
    return g


# ----------------------------------------------------------------------------
# Protein subset (the reuse / speed-up step)
# ----------------------------------------------------------------------------
def _iter_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.rstrip("\n"))
    if name is not None:
        yield name, "".join(seq)


def subset_proteins(work_dir, gp, ece_set, ref_fa, out_faa, threads):
    """Write ECE-only proteins. Reuse geNomad proteins.faa if present, else pyrodigal."""
    faa = os.path.join(work_dir, "Genomad", f"{gp}_annotate", f"{gp}_proteins.faa")
    if os.path.exists(faa):
        n = 0
        with open(out_faa, "w") as o:
            for pid, seq in _iter_fasta(faa):
                if au.gene_to_contig(pid) in ece_set:
                    o.write(f">{pid}\n{seq}\n")
                    n += 1
        return "reused geNomad proteins", n
    # fallback: extract ECE contigs, run pyrodigal-gv
    sub_fna = out_faa + ".contigs.fna"
    with open(out_faa + ".ece.list", "w") as lf:
        lf.write("\n".join(sorted(ece_set)) + "\n")
    subprocess.run(["seqkit", "grep", "-f", out_faa + ".ece.list", "-o", sub_fna, ref_fa],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    au.run_pyrodigal(sub_fna, out_faa, threads=threads)
    n = sum(1 for _ in _iter_fasta(out_faa))
    return "ran pyrodigal (geNomad proteins absent)", n


# ----------------------------------------------------------------------------
# Evidence lines
# ----------------------------------------------------------------------------
def ev_circularity(eces, gmd):
    rows = []
    for name in eces["seq_name"]:
        topo = gmd.get(name, {}).get("genomad_topology", "") or ""
        circ_h = name.endswith("C")
        tr = str(topo) in CIRCULAR_TOPOLOGY
        rows.append({"seq_name": name, "circular_hifiasm": circ_h,
                     "genomad_topology": topo, "terminal_repeat": tr,
                     "circular_support": int(circ_h) + int(tr)})
    return pd.DataFrame(rows)


def ev_conjscan(out_dir, ece_faa, threads):
    """Run ConjScan and classify each ECE contig (replicon)."""
    cj_dir = os.path.join(out_dir, "conjscan")
    summary = os.path.join(cj_dir, "best_solution_summary.tsv")
    all_sys = os.path.join(cj_dir, "all_systems.tsv")
    macsy = os.environ.get("MACSYFINDER",
                           "/home/shuaiw/miniconda3/envs/conjscan/bin/macsyfinder")
    if os.path.getsize(ece_faa) == 0:
        return pd.DataFrame(columns=["seq_name", "conjscan_type", "conjscan_anno"])
    if os.path.isdir(cj_dir):
        subprocess.run(["rm", "-rf", cj_dir], check=False)
    try:
        subprocess.run([macsy, "--db-type", "gembase", "-o", cj_dir,
                        "--sequence-db", ece_faa, "-w", str(threads), "-m", "CONJScan"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"[conjscan] failed: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["seq_name", "conjscan_type", "conjscan_anno"])
    if not os.path.exists(summary):
        return pd.DataFrame(columns=["seq_name", "conjscan_type", "conjscan_anno"])
    parsed = pcj.parse_conjscan(summary, os.path.join(cj_dir, "conjscan_parsed.tsv"))
    rows = []
    for _, r in parsed.iterrows():
        if r["conjugative"]:
            ctype = "Conjugative"
        elif r["mobilizable"]:
            ctype = "Mobilizable"
        elif r["dCONJ"]:
            ctype = "dCONJ"
        else:
            ctype = "None"
        rows.append({"seq_name": r["contig"], "conjscan_type": ctype,
                     "conjscan_anno": r["conjscan_anno"]})
    return pd.DataFrame(rows)


def _count_hits_per_contig(tblout, colname):
    df = au.parse_hmmer_tblout(tblout)
    if df.empty:
        return pd.Series(dtype=int), df
    df["contig"] = df["gene_id"].map(au.gene_to_contig)
    return df.groupby("contig")["gene_id"].nunique(), df


VIRAL_COLS = ["seq_name", "vir_terminase_large", "vir_terminase_small",
              "vir_major_capsid", "vir_portal", "vir_marker_total", "vog_hallmark_hits"]
SCMG_COLS = ["seq_name", "scmg_bac_count", "scmg_arc_count", "scmg_count", "scmg_fraction"]


def ev_viral(out_dir, ece_faa, db_dir, threads, use_vogdb):
    vir_hmm = os.path.join(db_dir, "viral_markers", "viral_markers.hmm")
    mapf = os.path.join(db_dir, "viral_markers", "marker_map.tsv")
    mp = pd.read_csv(mapf, sep="\t")
    acc2class = dict(zip(mp["accession"], mp["class"]))
    tbl = os.path.join(out_dir, "viral_pfam.tblout")
    au.run_hmmsearch(ece_faa, vir_hmm, tbl, threads=threads, cut_ga=True)
    df = au.parse_hmmer_tblout(tbl)
    classes = ["terminase_large", "terminase_small", "major_capsid", "portal"]
    percontig = {}
    if not df.empty:
        df["contig"] = df["gene_id"].map(au.gene_to_contig)
        df["class"] = df["query_acc"].map(acc2class)
        # distinct genes per (contig,class)
        for (contig, cls), sub in df.groupby(["contig", "class"]):
            percontig.setdefault(contig, {}).setdefault(cls, set()).update(sub["gene_id"])
    # VOGdb backstop
    vog_counts = {}
    if use_vogdb and os.path.exists(VOGDB_HMM) and os.path.getsize(ece_faa) > 0:
        vtbl = os.path.join(out_dir, "viral_vogdb.tblout")
        au.run_hmmsearch(ece_faa, VOGDB_HMM, vtbl, threads=threads, evalue=1e-5)
        s, _ = _count_hits_per_contig(vtbl, "vog")
        vog_counts = s.to_dict()
    rows = []
    for name in set(percontig) | set(vog_counts):
        cc = percontig.get(name, {})
        counts = {c: len(cc.get(c, set())) for c in classes}
        rows.append({
            "seq_name": name,
            "vir_terminase_large": counts["terminase_large"],
            "vir_terminase_small": counts["terminase_small"],
            "vir_major_capsid": counts["major_capsid"],
            "vir_portal": counts["portal"],
            "vir_marker_total": sum(counts.values()),
            "vog_hallmark_hits": int(vog_counts.get(name, 0)),
        })
    return pd.DataFrame(rows, columns=VIRAL_COLS)


def ev_scmg(out_dir, ece_faa, db_dir, threads):
    res = {}
    sizes = {"bac71": 71, "arc76": 76}
    for base in ("bac71", "arc76"):
        hmm = os.path.join(db_dir, "scmg", f"{base}.hmm")
        tbl = os.path.join(out_dir, f"scmg_{base}.tblout")
        au.run_hmmsearch(ece_faa, hmm, tbl, threads=threads, cut_ga=True)
        df = au.parse_hmmer_tblout(tbl)
        if df.empty:
            continue
        df["contig"] = df["gene_id"].map(au.gene_to_contig)
        # distinct marker models hit per contig
        per = df.groupby("contig")["query_name"].nunique()
        for contig, cnt in per.items():
            res.setdefault(contig, {})[base] = int(cnt)
    rows = []
    for contig, d in res.items():
        bac = d.get("bac71", 0)
        arc = d.get("arc76", 0)
        cnt = max(bac, arc)
        frac = max(bac / sizes["bac71"], arc / sizes["arc76"])
        rows.append({"seq_name": contig, "scmg_bac_count": bac, "scmg_arc_count": arc,
                     "scmg_count": cnt, "scmg_fraction": round(frac, 3)})
    return pd.DataFrame(rows, columns=SCMG_COLS)


def ev_genes_crosscheck(genes, ece_set):
    """Fast geNomad genes.tsv cross-checks per contig."""
    if genes is None:
        return pd.DataFrame(columns=["seq_name", "genes_uscg", "genes_virus_hallmark",
                                     "genes_conjscan_hits"])
    g = genes[genes["contig"].isin(ece_set)].copy()
    rows = []
    for contig, sub in g.groupby("contig"):
        uscg = pd.to_numeric(sub.get("uscg"), errors="coerce").fillna(0).sum()
        vh = pd.to_numeric(sub.get("virus_hallmark"), errors="coerce").fillna(0).sum()
        cj = sub.get("annotation_conjscan")
        cj_hits = int((cj.notna() & (cj.astype(str) != "NA")).sum()) if cj is not None else 0
        rows.append({"seq_name": contig, "genes_uscg": int(uscg),
                     "genes_virus_hallmark": int(vh), "genes_conjscan_hits": cj_hits})
    return pd.DataFrame(rows)


def ev_coverage_cv(bam, contigs):
    """Per-contig coverage CV from the BAM (mean/variance via pysam pileup)."""
    out = {}
    if not os.path.exists(bam):
        return out
    import pysam
    b = pysam.AlignmentFile(bam, "rb")
    refset = set(b.references)
    for c in contigs:
        if c not in refset:
            continue
        L = b.get_reference_length(c)
        cov = np.zeros(L, dtype=np.int32)
        for col in b.pileup(c, truncate=True):
            cov[col.reference_pos] = col.nsegments
        m = float(cov.mean()) if L else 0.0
        cv = float(np.sqrt(cov.var()) / m) if m > 0 else float("nan")
        out[c] = cv
    b.close()
    return out


# ----------------------------------------------------------------------------
# Confidence
# ----------------------------------------------------------------------------
def decide(row):
    t = row["type"]
    p1 = bool(row["circular_hifiasm"]) or bool(row["terminal_repeat"])
    p2 = bool(row["survives_default_threshold"]) and (row["genomad_n_hallmarks"] or 0) >= 1
    if t == "virus":
        p3 = (row["vir_marker_total"] or 0) >= 1 or (row["vog_hallmark_hits"] or 0) >= 1
    else:  # plasmid / novel
        p3 = (row["conjscan_type"] not in (None, "None", "")) or \
             (row["genomad_conjugation_genes"] or 0) > 0
    ratio = row["cov_ratio"]
    cv = row["cov_cv"]
    p4 = (pd.notna(ratio) and (ratio < COV_RATIO_LOW or ratio > COV_RATIO_HIGH)
          and pd.notna(cv) and cv < CV_MAX)
    n_pos = int(p1) + int(p2) + int(p3) + int(p4)

    flag_chr = (row["scmg_count"] or 0) >= SCMG_MIN or (row["scmg_fraction"] or 0) > SCMG_FRAC
    cov_mean = row["cov_mean"]
    flag_art = ((pd.notna(cov_mean) and cov_mean < COV_MIN)
                or (pd.notna(cv) and cv > CV_ARTIFACT)) and n_pos == 0

    if flag_chr:
        conf, reason = "Low", "likely_chromosomal_fragment"
    elif flag_art:
        conf, reason = "Low", "likely_artifact"
    elif n_pos >= 2:
        conf, reason = "High", "well_supported"
    elif n_pos == 1:
        conf, reason = "Medium", "single_line_support"
    else:
        conf, reason = "Low", "weak_support"
    return pd.Series({"support_circular": p1, "support_genomad": p2,
                      "support_marker": p3, "support_coverage": p4,
                      "n_positive_lines": n_pos, "flag_chromosomal": flag_chr,
                      "flag_artifact": flag_art, "confidence": conf,
                      "confidence_reason": reason})


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--work_dir", required=True)
    ap.add_argument("--ece_csv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--db_dir", default="/home/shuaiw/borg/revision/ece_anno/db")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--vogdb", dest="vogdb", action="store_true", default=True)
    ap.add_argument("--no-vogdb", dest="vogdb", action="store_false")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    gp = genomad_prefix(args.work_dir, args.sample)
    ref_fa = os.path.join(args.work_dir, f"{args.sample}.hifiasm.p_ctg.rename.fa")
    bam = os.path.join(args.work_dir, f"{args.sample}.align.bam")

    eces = load_eces(args.ece_csv, args.sample)
    if eces.empty:
        print(f"[{args.sample}] no ECEs in CSV; nothing to do.")
        return
    ece_set = set(eces["seq_name"])
    print(f"[{args.sample}] {len(eces)} ECEs "
          f"({(eces['type']=='plasmid').sum()} plasmid, {(eces['type']=='virus').sum()} virus)")

    gmd = read_genomad_summaries(args.work_dir, gp)
    genes = load_genes_tsv(args.work_dir, gp)

    ece_faa = os.path.join(args.out_dir, "ece_proteins.faa")
    src, n_prot = subset_proteins(args.work_dir, gp, ece_set, ref_fa, ece_faa, args.threads)
    print(f"[{args.sample}] proteins: {src} ({n_prot} ORFs)")

    # evidence frames
    df = eces.copy()
    df = df.merge(ev_circularity(eces, gmd), on="seq_name", how="left")
    df = df.merge(ev_conjscan(args.out_dir, ece_faa, args.threads), on="seq_name", how="left")
    # VOGdb is the expensive scan; only run it when a virus-type ECE is present.
    run_vogdb = args.vogdb and bool((eces["type"] == "virus").any())
    df = df.merge(ev_viral(args.out_dir, ece_faa, args.db_dir, args.threads, run_vogdb),
                  on="seq_name", how="left")
    df = df.merge(ev_scmg(args.out_dir, ece_faa, args.db_dir, args.threads),
                  on="seq_name", how="left")
    df = df.merge(ev_genes_crosscheck(genes, ece_set), on="seq_name", how="left")

    # geNomad robustness columns (genomad_topology already provided by ev_circularity)
    gsum = pd.DataFrame([{"seq_name": k, "genomad_score": v["genomad_score"],
                          "genomad_fdr": v["genomad_fdr"],
                          "genomad_n_hallmarks": v["genomad_n_hallmarks"],
                          "genomad_conjugation_genes": v["genomad_conjugation_genes"]}
                         for k, v in gmd.items()])
    if not gsum.empty:
        df = df.merge(gsum, on="seq_name", how="left")
    for c in ["genomad_score", "genomad_fdr", "genomad_n_hallmarks",
              "genomad_conjugation_genes"]:
        if c not in df.columns:
            df[c] = np.nan
    df["survives_default_threshold"] = df["genomad_fdr"] <= GENOMAD_FDR_MAX

    # coverage
    df["cov_mean"] = pd.to_numeric(df["mge_depth"], errors="coerce")
    df["host_depth"] = pd.to_numeric(df["host_depth"], errors="coerce")
    df["cov_ratio"] = df["cov_mean"] / df["host_depth"]
    cv_map = ev_coverage_cv(bam, list(ece_set))
    df["cov_cv"] = df["seq_name"].map(cv_map)

    # fill / types
    int_cols = ["circular_support", "vir_terminase_large", "vir_terminase_small",
                "vir_major_capsid", "vir_portal", "vir_marker_total", "vog_hallmark_hits",
                "scmg_bac_count", "scmg_arc_count", "scmg_count", "genes_uscg",
                "genes_virus_hallmark", "genes_conjscan_hits"]
    for c in int_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df["scmg_fraction"] = pd.to_numeric(df.get("scmg_fraction"), errors="coerce").fillna(0.0)
    for c in ["circular_hifiasm", "terminal_repeat"]:
        df[c] = df[c].fillna(False).astype(bool)
    df["conjscan_type"] = df["conjscan_type"].fillna("None")
    df["conjscan_anno"] = df.get("conjscan_anno", "").fillna("")
    df["genomad_conjugation_genes"] = pd.to_numeric(
        df["genomad_conjugation_genes"], errors="coerce").fillna(0)
    df["genomad_n_hallmarks"] = pd.to_numeric(
        df["genomad_n_hallmarks"], errors="coerce").fillna(0)

    df = pd.concat([df, df.apply(decide, axis=1)], axis=1)

    cols = ["seq_name", "type", "length", "host_lineage",
            "circular_hifiasm", "genomad_topology", "terminal_repeat", "circular_support",
            "conjscan_type", "conjscan_anno", "genomad_conjugation_genes", "genes_conjscan_hits",
            "vir_terminase_large", "vir_terminase_small", "vir_major_capsid", "vir_portal",
            "vir_marker_total", "vog_hallmark_hits", "genes_virus_hallmark",
            "scmg_bac_count", "scmg_arc_count", "scmg_count", "scmg_fraction", "genes_uscg",
            "cov_mean", "cov_cv", "host_depth", "cov_ratio",
            "genomad_score", "genomad_fdr", "genomad_n_hallmarks", "survives_default_threshold",
            "support_circular", "support_genomad", "support_marker", "support_coverage",
            "n_positive_lines", "flag_chromosomal", "flag_artifact", "confidence",
            "confidence_reason"]
    df = df[[c for c in cols if c in df.columns]]
    out_tsv = os.path.join(args.out_dir, "ece_evidence.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[{args.sample}] wrote {out_tsv}")
    print(df["confidence"].value_counts().to_string())


if __name__ == "__main__":
    main()
