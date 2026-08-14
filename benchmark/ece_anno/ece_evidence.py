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
GENOMAD_SCORE_MIN = 0.8     # strict (conservative-like) score floor for P2
GENOMAD_FDR_MAX = 0.05      # strict FDR cutoff for P2
COV_RATIO_LOW = 0.75        # coverage-distinctness band (P4)
COV_RATIO_HIGH = 1.33
CV_MAX = 0.75               # max CV for "uniform" coverage (P4)
CV_ARTIFACT = 1.5           # CV above this contributes to artifact flag
COV_MIN = 1.0               # mean depth below this contributes to artifact flag
SCMG_MIN = 5                # >= this many chromosomal single-copy markers -> chromosomal flag
SCMG_FRAC = 0.10            # or this fraction of a marker set
VIR_MIN_CLASSES = 2         # virus P3 needs >= this many distinct structural classes
VOGDB_HMM = "/shared/db/vogdb/latest/hmm/VOG.all.hmm"

# Circularity = a closed circle only. DTR (direct terminal repeat) => circular;
# ITR (inverted terminal repeat) => a COMPLETE LINEAR genome (not circular).
CIRCULAR_TOPOLOGY = {"DTR"}
COMPLETE_LINEAR_TOPOLOGY = {"ITR"}


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


def load_eces(ece_csv, sample, cols):
    """ECE rows for this sample from the linkage CSV, deduped by contig.

    `cols` maps logical fields -> CSV column names (see --col_* args), so the same
    engine handles the isolate CSV (mge_contig/prefix/...) and the metagenome CSV
    (MGE/sample/MGE_cov/...).
    """
    df = pd.read_csv(ece_csv)
    df = df[df[cols["sample"]].astype(str) == str(sample)].copy()
    df = df[~df[cols["seqname"]].astype(str).str.contains(r"\|provirus", regex=True)]
    df = df.drop_duplicates(cols["seqname"])
    out = pd.DataFrame({
        "seq_name": df[cols["seqname"]].astype(str),
        "type": df[cols["type"]].astype(str),
        "length": df[cols["length"]],
        "mge_depth": df[cols["mgedepth"]],
        "host_depth": df[cols["hostdepth"]],
        "host_lineage": df.get(cols["lineage"], ""),
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
    """geNomad per-gene annotations, read from the ECE-only summary genes tables
    ({gp}_plasmid_genes.tsv / {gp}_virus_genes.tsv) so we never parse the whole
    metagenome's annotate/genes.tsv. None if absent."""
    sdir = os.path.join(work_dir, "Genomad", f"{gp}_summary")
    frames = []
    for kind in ("plasmid", "virus"):
        p = os.path.join(sdir, f"{gp}_{kind}_genes.tsv")
        if os.path.exists(p):
            frames.append(pd.read_csv(p, sep="\t"))
    if not frames:
        return None
    g = pd.concat(frames, ignore_index=True).drop_duplicates("gene")
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


def _extract_ece_contigs(work_dir, sample, gp, contigs, ref_fa, out_fna):
    """Collect nucleotide sequence for the given ECE contigs only, preferring the
    existing per-contig FASTAs (<sample>_methylation4/contigs/<contig>.fa) and
    falling back to a name-indexed seqkit grep from the assembly. Never scans the
    whole metagenome sequence linearly."""
    contig_dir = os.path.join(work_dir, f"{sample}_methylation4", "contigs")
    still = []
    n = 0
    with open(out_fna, "w") as o:
        for c in contigs:
            pcf = os.path.join(contig_dir, f"{c}.fa")
            if os.path.exists(pcf):
                for name, seq in _iter_fasta(pcf):
                    o.write(f">{c}\n{seq}\n")
                n += 1
            else:
                still.append(c)
    if still and os.path.exists(ref_fa):  # random-access fetch by name via .fai
        lst = out_fna + ".names"
        with open(lst, "w") as lf:
            lf.write("\n".join(still) + "\n")
        with open(out_fna, "a") as o:
            subprocess.run(["seqkit", "grep", "-f", lst, ref_fa],
                           check=True, stdout=o, stderr=subprocess.DEVNULL)
        n += len(still)
    return n


def subset_proteins(work_dir, sample, gp, ece_set, ref_fa, out_faa, threads):
    """Write proteins for the ECE contigs ONLY (never the whole metagenome).

    Reuses geNomad's ECE-only summary protein files
    ({gp}_summary/{gp}_{plasmid,virus}_proteins.faa); for any ECE contig not covered
    there (e.g. novel elements geNomad did not classify), predicts ORFs with
    pyrodigal-gv on that contig's own sequence only.
    """
    sdir = os.path.join(work_dir, "Genomad", f"{gp}_summary")
    src_files = [os.path.join(sdir, f"{gp}_plasmid_proteins.faa"),
                 os.path.join(sdir, f"{gp}_virus_proteins.faa")]
    covered, n = set(), 0
    with open(out_faa, "w") as o:
        for sf in src_files:
            if not os.path.exists(sf):
                continue
            for pid, seq in _iter_fasta(sf):
                c = au.gene_to_contig(pid)
                if c in ece_set:
                    o.write(f">{pid}\n{seq}\n")
                    n += 1
                    covered.add(c)
    src = f"reused geNomad ECE proteins ({len(covered)} contigs)"
    missing = ece_set - covered
    if missing:
        sub_fna = out_faa + ".missing.fna"
        tmp_faa = out_faa + ".missing.faa"
        if _extract_ece_contigs(work_dir, sample, gp, missing, ref_fa, sub_fna):
            au.run_pyrodigal(sub_fna, tmp_faa, threads=threads)
            with open(out_faa, "a") as o, open(tmp_faa) as fin:
                for line in fin:
                    o.write(line)
            n += sum(1 for _ in _iter_fasta(tmp_faa))
            src += f" + pyrodigal for {len(missing)} uncovered contig(s)"
    return src, n


# ----------------------------------------------------------------------------
# Evidence lines
# ----------------------------------------------------------------------------
def ev_circularity(eces, gmd):
    """Completeness signals. circular = hifiasm graph-cycle OR DTR; complete_linear = ITR
    (a complete linear genome, NOT circular). Either indicates a complete, non-fragment
    element and satisfies the P1 completeness line."""
    rows = []
    for name in eces["seq_name"]:
        topo = str(gmd.get(name, {}).get("genomad_topology", "") or "")
        circ_h = name.endswith("C")
        circular = circ_h or (topo in CIRCULAR_TOPOLOGY)
        complete_linear = topo in COMPLETE_LINEAR_TOPOLOGY
        rows.append({"seq_name": name, "circular_hifiasm": circ_h,
                     "genomad_topology": topo, "circular": circular,
                     "complete_linear": complete_linear,
                     "complete_support": int(circular) + int(complete_linear)})
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
    if not os.path.exists(summary):  # reuse existing ConjScan run if present
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
              "vir_major_capsid", "vir_portal", "vir_marker_total", "vir_n_classes",
              "vog_hallmark_hits"]
SCMG_COLS = ["seq_name", "scmg_bac_count", "scmg_arc_count", "scmg_count", "scmg_fraction"]
PLAS_COLS = ["seq_name", "plas_rep_initiation", "plas_partition", "plas_marker_total"]


def ev_viral(out_dir, ece_faa, db_dir, threads, use_vogdb):
    """Detect the four phage structural classes (terminase L/S, major capsid, portal)
    per contig. A class is present if EITHER the independent curated Pfam profiles OR
    VOGdb (mapped by consensus function) detect it. vir_n_classes = size of that union;
    this drives the virus P3 line (>= VIR_MIN_CLASSES)."""
    classes = ["terminase_large", "terminase_small", "major_capsid", "portal"]
    vir_hmm = os.path.join(db_dir, "viral_markers", "viral_markers.hmm")
    mp = pd.read_csv(os.path.join(db_dir, "viral_markers", "marker_map.tsv"), sep="\t")
    acc2class = dict(zip(mp["accession"], mp["class"]))
    tbl = os.path.join(out_dir, "viral_pfam.tblout")
    au.run_hmmsearch(ece_faa, vir_hmm, tbl, threads=threads, cut_ga=True, reuse=True)
    df = au.parse_hmmer_tblout(tbl)
    pfam = {}  # contig -> {class: set(genes)}
    if not df.empty:
        df["contig"] = df["gene_id"].map(au.gene_to_contig)
        df["class"] = df["query_acc"].map(acc2class)
        for (contig, cls), sub in df.groupby(["contig", "class"]):
            pfam.setdefault(contig, {}).setdefault(cls, set()).update(sub["gene_id"])

    # VOGdb structural classes (reuse tblout if present; map VOG -> class by function)
    vog_cls = {}       # contig -> set(class)
    vog_counts = {}    # contig -> total VOG hits (info)
    vmapf = os.path.join(db_dir, "viral_markers", "vog_structural_map.tsv")
    vtbl = os.path.join(out_dir, "viral_vogdb.tblout")
    if (use_vogdb and os.path.exists(vmapf)
            and (os.path.exists(vtbl) or os.path.getsize(ece_faa) > 0)):
        if os.path.exists(VOGDB_HMM):
            au.run_hmmsearch(ece_faa, VOGDB_HMM, vtbl, threads=threads, evalue=1e-5, reuse=True)
        vog2c = dict(zip(*[pd.read_csv(vmapf, sep="\t")[c] for c in ("vog", "class")]))
        vh = au.parse_hmmer_tblout(vtbl)
        if not vh.empty:
            vh["contig"] = vh["gene_id"].map(au.gene_to_contig)
            for contig, sub in vh.groupby("contig"):
                vog_counts[contig] = int(sub["gene_id"].nunique())
                cls = {vog2c[q] for q in sub["query_name"] if q in vog2c}
                if cls:
                    vog_cls[contig] = cls

    rows = []
    for name in set(pfam) | set(vog_counts):
        pf = {c: len(pfam.get(name, {}).get(c, set())) for c in classes}
        present = {c for c in classes if pf[c] > 0} | vog_cls.get(name, set())
        rows.append({
            "seq_name": name,
            "vir_terminase_large": pf["terminase_large"],
            "vir_terminase_small": pf["terminase_small"],
            "vir_major_capsid": pf["major_capsid"],
            "vir_portal": pf["portal"],
            "vir_marker_total": sum(pf.values()),
            "vir_n_classes": len(present),  # union of Pfam + VOGdb structural classes
            "vog_hallmark_hits": int(vog_counts.get(name, 0)),
        })
    return pd.DataFrame(rows, columns=VIRAL_COLS)


def ev_plasmid_markers(out_dir, ece_faa, db_dir, threads):
    """Plasmid hallmark markers: replication-initiation (Rep) + partition (ParA/ParB)."""
    hmm = os.path.join(db_dir, "plasmid_markers", "plasmid_markers.hmm")
    mapf = os.path.join(db_dir, "plasmid_markers", "marker_map.tsv")
    if not os.path.exists(hmm):
        return pd.DataFrame(columns=PLAS_COLS)
    acc2class = dict(zip(*[pd.read_csv(mapf, sep="\t")[c] for c in ("accession", "class")]))
    tbl = os.path.join(out_dir, "plasmid_pfam.tblout")
    au.run_hmmsearch(ece_faa, hmm, tbl, threads=threads, cut_ga=True, reuse=True)
    df = au.parse_hmmer_tblout(tbl)
    rows = []
    if not df.empty:
        df["contig"] = df["gene_id"].map(au.gene_to_contig)
        df["class"] = df["query_acc"].map(acc2class)
        for contig, sub in df.groupby("contig"):
            rep = sub.loc[sub["class"] == "rep_initiation", "gene_id"].nunique()
            par = sub.loc[sub["class"] == "partition", "gene_id"].nunique()
            rows.append({"seq_name": contig, "plas_rep_initiation": int(rep),
                         "plas_partition": int(par), "plas_marker_total": int(rep + par)})
    return pd.DataFrame(rows, columns=PLAS_COLS)


def ev_rrna(out_dir, ece_fna, db_dir, threads):
    """Detect rRNA genes (16S/23S) on the ECE contigs via cmsearch vs Rfam CMs.
    Presence of rRNA is chromosomal evidence (ECEs are depleted of rRNA)."""
    cm = os.path.join(db_dir, "rrna", "rrna.cm")
    if not os.path.exists(cm) or not os.path.exists(ece_fna) or os.path.getsize(ece_fna) == 0:
        return {}
    tbl = os.path.join(out_dir, "rrna.tblout")
    if not (os.path.exists(tbl) and os.path.getsize(tbl) > 0):  # reuse if present
        try:
            subprocess.run(["cmsearch", "--cut_ga", "--noali", "--cpu", str(threads),
                            "--tblout", tbl, cm, ece_fna],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"[rrna] cmsearch failed: {e}", file=sys.stderr)
            return {}
    counts = {}
    if os.path.exists(tbl):
        with open(tbl) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                # tblout: col0 = target (contig) name
                counts[p[0]] = counts.get(p[0], 0) + 1
    return counts


def ev_scmg(out_dir, ece_faa, db_dir, threads):
    res = {}
    sizes = {"bac71": 71, "arc76": 76}
    for base in ("bac71", "arc76"):
        hmm = os.path.join(db_dir, "scmg", f"{base}.hmm")
        tbl = os.path.join(out_dir, f"scmg_{base}.tblout")
        au.run_hmmsearch(ece_faa, hmm, tbl, threads=threads, cut_ga=True, reuse=True)
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


def _cv_from_pileup(bam_file, contig, region=None):
    import pysam
    b = pysam.AlignmentFile(bam_file, "rb")
    try:
        L = b.get_reference_length(contig)
    except (KeyError, ValueError):
        L = max(b.lengths) if b.lengths else 0
    if not L:
        b.close()
        return float("nan")
    cov = np.zeros(L, dtype=np.int32)
    it = b.pileup(contig, truncate=True) if region else b.pileup(truncate=True)
    for col in it:
        if 0 <= col.reference_pos < L:
            cov[col.reference_pos] = col.nsegments
    b.close()
    m = float(cov.mean())
    return float(np.sqrt(cov.var()) / m) if m > 0 else float("nan")


def ev_coverage_cv(work_dir, sample, contigs, fallback_bam):
    """Per-contig coverage CV (sqrt(var)/mean of per-base depth).

    Prefers the small per-contig BAMs under <work_dir>/<sample>_methylation4/bams/,
    falling back to the whole-sample alignment only if a per-contig BAM is absent.
    """
    bam_dir = os.path.join(work_dir, f"{sample}_methylation4", "bams")
    out = {}
    for c in contigs:
        pcb = os.path.join(bam_dir, f"{c}.bam")
        if os.path.exists(pcb):
            out[c] = _cv_from_pileup(pcb, c, region=False)
        elif os.path.exists(fallback_bam):
            out[c] = _cv_from_pileup(fallback_bam, c, region=True)
        else:
            out[c] = float("nan")
    return out


# ----------------------------------------------------------------------------
# Confidence
# ----------------------------------------------------------------------------
def decide(row):
    t = row["type"]
    # P1 completeness: circular (hifiasm cycle or DTR) OR complete linear (ITR).
    p1 = bool(row["circular"]) or bool(row["complete_linear"])
    # P2 geNomad robustness (strict): calibrated score >= 0.8 AND FDR <= 0.05.
    p2 = (row["genomad_score"] or 0) >= GENOMAD_SCORE_MIN and \
         pd.notna(row["genomad_fdr"]) and row["genomad_fdr"] <= GENOMAD_FDR_MAX
    # P3 type-specific marker (stringent).
    if t == "virus":
        # >= 2 distinct structural classes (terminase L/S, capsid, portal) — avoids
        # single phage-tail hallmark false positives.
        p3 = (row["vir_n_classes"] or 0) >= VIR_MIN_CLASSES
    else:  # plasmid / novel
        p3 = (row["conjscan_type"] not in (None, "None", "")) or \
             (row["genomad_conjugation_genes"] or 0) > 0 or \
             (row["plas_marker_total"] or 0) > 0
    ratio = row["cov_ratio"]
    cv = row["cov_cv"]
    p4 = (pd.notna(ratio) and (ratio < COV_RATIO_LOW or ratio > COV_RATIO_HIGH)
          and pd.notna(cv) and cv < CV_MAX)
    n_pos = int(p1) + int(p2) + int(p3) + int(p4)

    # chromosomal fragment: many single-copy core genes OR presence of rRNA genes.
    flag_chr = ((row["scmg_count"] or 0) >= SCMG_MIN or (row["scmg_fraction"] or 0) > SCMG_FRAC
                or (row["rrna_count"] or 0) > 0)
    cov_mean = row["cov_mean"]
    flag_art = ((pd.notna(cov_mean) and cov_mean < COV_MIN)
                or (pd.notna(cv) and cv > CV_ARTIFACT)) and n_pos == 0

    # very-high-confidence ECE = ALL four independent lines pass and no negative flag.
    # This is a deliberately stringent, high-precision set for validating the
    # ECE-host linkage method (recall is not the objective here).
    very_high = p1 and p2 and p3 and p4 and (not flag_chr) and (not flag_art)
    return pd.Series({"support_circular": p1, "support_genomad": p2,
                      "support_marker": p3, "support_coverage": p4,
                      "n_positive_lines": n_pos, "flag_chromosomal": flag_chr,
                      "flag_artifact": flag_art, "very_high_confidence": very_high})


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--work_dir", required=True)
    ap.add_argument("--ece_csv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--db_dir", default="/home/shuaiw/borg/revision/ece_anno/db")
    ap.add_argument("--threads", type=int, default=8)
    # VOGdb contributes to viral structural-class detection (union with Pfam); the
    # existing tblouts are reused, so this is cheap. Disable with --no-vogdb.
    ap.add_argument("--vogdb", dest="vogdb", action="store_true", default=True)
    ap.add_argument("--no-vogdb", dest="vogdb", action="store_false")
    # CSV column mapping (defaults = isolate jaccard_same_sample.csv schema)
    ap.add_argument("--col_sample", default="prefix")
    ap.add_argument("--col_seqname", default="mge_contig")
    ap.add_argument("--col_type", default="mge_type")
    ap.add_argument("--col_length", default="mge_length")
    ap.add_argument("--col_mgedepth", default="mge_depth")
    ap.add_argument("--col_hostdepth", default="host_depth")
    ap.add_argument("--col_lineage", default="host_lineage")
    args = ap.parse_args()
    cols = {"sample": args.col_sample, "seqname": args.col_seqname, "type": args.col_type,
            "length": args.col_length, "mgedepth": args.col_mgedepth,
            "hostdepth": args.col_hostdepth, "lineage": args.col_lineage}

    os.makedirs(args.out_dir, exist_ok=True)
    gp = genomad_prefix(args.work_dir, args.sample)
    ref_fa = os.path.join(args.work_dir, f"{args.sample}.hifiasm.p_ctg.rename.fa")
    bam = os.path.join(args.work_dir, f"{args.sample}.align.bam")

    eces = load_eces(args.ece_csv, args.sample, cols)
    if eces.empty:
        print(f"[{args.sample}] no ECEs in CSV; nothing to do.")
        return
    ece_set = set(eces["seq_name"])
    print(f"[{args.sample}] {len(eces)} ECEs "
          f"({(eces['type']=='plasmid').sum()} plasmid, {(eces['type']=='virus').sum()} virus)")

    gmd = read_genomad_summaries(args.work_dir, gp)
    genes = load_genes_tsv(args.work_dir, gp)

    ece_faa = os.path.join(args.out_dir, "ece_proteins.faa")
    src, n_prot = subset_proteins(args.work_dir, args.sample, gp, ece_set, ref_fa,
                                  ece_faa, args.threads)
    print(f"[{args.sample}] proteins: {src} ({n_prot} ORFs)")
    # ECE contig nucleotides (for rRNA cmsearch) — ECE contigs only, never the whole assembly
    ece_fna = os.path.join(args.out_dir, "ece_contigs.fna")
    _extract_ece_contigs(args.work_dir, args.sample, gp, ece_set, ref_fa, ece_fna)

    # evidence frames
    df = eces.copy()
    df = df.merge(ev_circularity(eces, gmd), on="seq_name", how="left")
    df = df.merge(ev_conjscan(args.out_dir, ece_faa, args.threads), on="seq_name", how="left")
    # VOGdb is off by default (info only); only scan when a virus ECE is present.
    run_vogdb = args.vogdb and bool((eces["type"] == "virus").any())
    df = df.merge(ev_viral(args.out_dir, ece_faa, args.db_dir, args.threads, run_vogdb),
                  on="seq_name", how="left")
    df = df.merge(ev_plasmid_markers(args.out_dir, ece_faa, args.db_dir, args.threads),
                  on="seq_name", how="left")
    df = df.merge(ev_scmg(args.out_dir, ece_faa, args.db_dir, args.threads),
                  on="seq_name", how="left")
    df = df.merge(ev_genes_crosscheck(genes, ece_set), on="seq_name", how="left")
    # rRNA presence (chromosomal signal)
    rrna_map = ev_rrna(args.out_dir, ece_fna, args.db_dir, args.threads)
    df["rrna_count"] = df["seq_name"].map(rrna_map)

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
    df["genomad_score"] = pd.to_numeric(df["genomad_score"], errors="coerce")
    df["genomad_fdr"] = pd.to_numeric(df["genomad_fdr"], errors="coerce")
    df["survives_strict"] = (df["genomad_score"] >= GENOMAD_SCORE_MIN) & \
                            (df["genomad_fdr"] <= GENOMAD_FDR_MAX)

    # coverage
    df["cov_mean"] = pd.to_numeric(df["mge_depth"], errors="coerce")
    df["host_depth"] = pd.to_numeric(df["host_depth"], errors="coerce")
    df["cov_ratio"] = df["cov_mean"] / df["host_depth"]
    cv_map = ev_coverage_cv(args.work_dir, args.sample, list(ece_set), bam)
    df["cov_cv"] = df["seq_name"].map(cv_map)

    # fill / types
    int_cols = ["complete_support", "vir_terminase_large", "vir_terminase_small",
                "vir_major_capsid", "vir_portal", "vir_marker_total", "vir_n_classes",
                "vog_hallmark_hits", "plas_rep_initiation", "plas_partition",
                "plas_marker_total", "scmg_bac_count", "scmg_arc_count", "scmg_count",
                "genes_uscg", "genes_virus_hallmark", "genes_conjscan_hits", "rrna_count"]
    for c in int_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df["scmg_fraction"] = pd.to_numeric(df.get("scmg_fraction"), errors="coerce").fillna(0.0)
    for c in ["circular_hifiasm", "circular", "complete_linear"]:
        df[c] = df[c].fillna(False).astype(bool)
    df["conjscan_type"] = df["conjscan_type"].fillna("None")
    df["conjscan_anno"] = df.get("conjscan_anno", "").fillna("")
    df["genomad_conjugation_genes"] = pd.to_numeric(
        df["genomad_conjugation_genes"], errors="coerce").fillna(0)
    df["genomad_n_hallmarks"] = pd.to_numeric(
        df["genomad_n_hallmarks"], errors="coerce").fillna(0)

    df = pd.concat([df, df.apply(decide, axis=1)], axis=1)

    cols = ["seq_name", "type", "length", "host_lineage",
            "circular_hifiasm", "genomad_topology", "circular", "complete_linear",
            "complete_support",
            "conjscan_type", "conjscan_anno", "genomad_conjugation_genes", "genes_conjscan_hits",
            "plas_rep_initiation", "plas_partition", "plas_marker_total",
            "vir_terminase_large", "vir_terminase_small", "vir_major_capsid", "vir_portal",
            "vir_marker_total", "vir_n_classes", "vog_hallmark_hits", "genes_virus_hallmark",
            "scmg_bac_count", "scmg_arc_count", "scmg_count", "scmg_fraction", "genes_uscg",
            "rrna_count",
            "cov_mean", "cov_cv", "host_depth", "cov_ratio",
            "genomad_score", "genomad_fdr", "genomad_n_hallmarks", "survives_strict",
            "support_circular", "support_genomad", "support_marker", "support_coverage",
            "n_positive_lines", "flag_chromosomal", "flag_artifact", "very_high_confidence"]
    df = df[[c for c in cols if c in df.columns]]
    out_tsv = os.path.join(args.out_dir, "ece_evidence.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[{args.sample}] wrote {out_tsv}  "
          f"(very_high_confidence: {int(df['very_high_confidence'].sum())}/{len(df)})")


if __name__ == "__main__":
    main()
