#!/usr/bin/env python
"""
Build MTase presence/absence and activity matrices for the Klebsiella pneumoniae
248_1 strain cluster (15 MAGs from the 59-metagenome dataset).

For each MAG we parse the MicrobeMod RM-gene annotation table, keep only the
methyltransferase genes (Gene type MT / IIG), and group them by their REBASE
homolog (the annotated recognition motif follows the homolog). We then build:

  - presence/absence matrix  (MAG x MTase, 1 if >=1 copy)
  - copy-number matrix       (MAG x MTase, number of gene copies)
  - activity matrix          (MAG x MTase, 1 if the MTase recognition motif is
                              detected as modified in that MAG, else 0, NA if the
                              MTase has no REBASE motif to test)

The activity call cross-references each MTase recognition motif against the
detected-methylation motifs of the same MAG (the same data underlying Fig 3e),
using an IUPAC-degenerate containment test.

Outputs (source-data CSVs) go to tmp/rev_figs/Klevsuekka_mtase/.
"""

import os
import re
import subprocess
import pandas as pd

# ----------------------------------------------------------------------------
# Paths / configuration
# ----------------------------------------------------------------------------
RUN2 = "/home/shuaiw/borg/paper/run2"
METHY = "methylation4"            # canonical run
OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/Klevsuekka_mtase"

# mmseqs work dir (extracted FASTA + clustering outputs). NOT on the small /home SSD.
WORK_DIR = "/home/shuaiw/borg/revision/MTase"
MMSEQS = "/shared/software/bin/mmseqs"
MIN_SEQ_ID = 0.9                  # ortholog identity threshold (~99% ANI genomes)
MIN_COV = 0.8

# The 15 MAGs of dRep secondary cluster 248_1 (Klebsiella pneumoniae).
MAGS = [
    "96plex_11_C", "infant_1_1_C", "infant_2_6_C", "infant_2_11_C",
    "infant_3_4_C", "infant_3_25_C", "infant_4_8_C", "infant_8_9_C",
    "infant_14_14_C", "infant_15_6_C", "infant_16_9_L", "infant_19_5_L",
    "infant_20_6_C", "infant_25_32_C", "infant_27_13_C",
]

# A detected motif counts as "active" at/above this modified fraction.
ACTIVITY_FRACTION = 0.25

# ----------------------------------------------------------------------------
# IUPAC-aware motif compatibility
# ----------------------------------------------------------------------------
IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"),
    "H": set("ACT"), "V": set("ACG"), "N": set("ACGT"),
}


def bases_compatible(a, b):
    """True if two IUPAC codes can represent a common concrete base."""
    return len(IUPAC.get(a, set()) & IUPAC.get(b, set())) > 0


COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def revcomp(motif):
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(motif))


def _equiv_same_len(a, b):
    """Same length and every position IUPAC-compatible."""
    if len(a) != len(b):
        return False
    return all(bases_compatible(a[i], b[i]) for i in range(len(a)))


def motifs_equivalent(mtase_motif, detected):
    """
    True only when the detected motif describes the SAME recognition site as the
    MTase's recognition motif: identical length and IUPAC-compatible at every
    position (checked on both strands). This deliberately does NOT accept a
    longer, more-specific detected motif that merely CONTAINS the MTase motif -
    e.g. Dam's 'GATC' is called active only when bare 'GATC' is detected, not
    when the methylation is confined to an extended motif such as 'GATCWADWD'
    (which reflects a different / more restricted specificity). This keeps the
    activity call consistent with the main-text motif-fraction figure.
    """
    m = str(mtase_motif).upper().strip()
    d = str(detected).upper().strip()
    if not m or not d:
        return False
    return _equiv_same_len(m, d) or _equiv_same_len(m, revcomp(d))


# ----------------------------------------------------------------------------
# Parsing helpers (idiom from scripts/analyze_RM.py)
# ----------------------------------------------------------------------------
def sample_of(mag):
    """MAG name -> sample dir (strip the last two '_'-fields)."""
    return "_".join(mag.split("_")[:-2])


def rm_table_path(mag):
    s = sample_of(mag)
    return os.path.join(RUN2, s, f"{s}_{METHY}", "RM_systems", "all_ctgs_RM.rm.genes.tsv")


def motifs_path(mag):
    s = sample_of(mag)
    return os.path.join(RUN2, s, f"{s}_{METHY}", "motifs", f"{mag}.motifs.csv")


def faa_path(mag):
    s = sample_of(mag)
    return os.path.join(RUN2, s, f"{s}_{METHY}", "RM_systems", "all_ctgs_RM.rm.genes.faa")


def load_mtases(mag):
    """Return the MTase (MT/IIG) genes for one MAG, deduplicated by Operon."""
    df = pd.read_csv(rm_table_path(mag), sep="\t", dtype=str)
    df = df[df["Gene"].str.startswith(mag + "_", na=False)]
    df = df[df["Gene type"].isin(["MT", "IIG"])]
    # One MTase gene can produce several domain hits sharing an Operon; keep one.
    df = df.drop_duplicates(subset=["Operon"]).copy()
    df["Homolog motif"] = df["Homolog motif"].fillna("").str.strip()
    df["REBASE homolog"] = df["REBASE homolog"].fillna("").str.strip()
    df["Predicted methylation"] = df["Predicted methylation"].fillna("").str.strip()
    df["MAG"] = mag
    return df


def read_faa(path):
    """Parse a protein FASTA -> {gene_id: sequence}. Gene id = 1st header token."""
    seqs, gid, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if gid is not None:
                    seqs[gid] = "".join(buf)
                gid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if gid is not None:
        seqs[gid] = "".join(buf)
    return seqs


def cluster_orthologs(gene_ids):
    """
    Extract protein sequences for `gene_ids` from the per-sample .faa files, write
    a combined FASTA to WORK_DIR, run mmseqs easy-cluster, and return a dict
    gene_id -> cluster representative gene_id.
    """
    os.makedirs(WORK_DIR, exist_ok=True)
    # cache faa per sample so each file is read once
    faa_cache = {}
    fasta = os.path.join(WORK_DIR, "kp248_mtase.faa")
    missing = []
    with open(fasta, "w") as out:
        for g in sorted(gene_ids):
            mag = "_".join(g.split("_")[:-1])         # contig = gene minus last field
            fp = faa_path(mag)
            if fp not in faa_cache:
                faa_cache[fp] = read_faa(fp)
            seq = faa_cache[fp].get(g)
            if seq is None:
                missing.append(g)
                continue
            out.write(f">{g}\n{seq}\n")
    if missing:
        print(f"WARNING: {len(missing)} gene(s) had no protein sequence: {missing}")

    prefix = os.path.join(WORK_DIR, "kp248_mtase_clu")
    tmp = os.path.join(WORK_DIR, "mmseqs_tmp")
    cmd = [MMSEQS, "easy-cluster", fasta, prefix, tmp,
           "--min-seq-id", str(MIN_SEQ_ID), "-c", str(MIN_COV),
           "--cov-mode", "0", "-v", "1"]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    gene2rep = {}
    with open(prefix + "_cluster.tsv") as f:
        for line in f:
            rep, member = line.rstrip("\n").split("\t")
            gene2rep[member] = rep
    return gene2rep


def load_detected_motifs(mag):
    """Return list of (motifString, fraction) detected as modified in the MAG."""
    p = motifs_path(mag)
    if not os.path.exists(p):
        return []
    m = pd.read_csv(p)
    out = []
    for _, r in m.iterrows():
        try:
            frac = float(r["fraction"])
        except (ValueError, TypeError):
            frac = 0.0
        out.append((str(r["motifString"]).upper().strip(), frac))
    return out


def clean_systype(s):
    return str(s).replace("RM_Type_", "Type ") if isinstance(s, str) else s


def first_nonempty(series):
    for v in series:
        if isinstance(v, str) and v.strip():
            return v.strip()
    return ""


def most_common_nonempty(series):
    vals = [v.strip() for v in series if isinstance(v, str) and v.strip()]
    if not vals:
        return ""
    return pd.Series(vals).mode().iloc[0]


def cluster_label(genes):
    """Human label for an ortholog cluster: representative REBASE homolog + motif."""
    homolog = most_common_nonempty(genes["REBASE homolog"])
    motif = first_nonempty(genes["Homolog motif"])
    if not homolog:
        homolog = f"[HMM] {most_common_nonempty(genes['HMM'])}"
    base = re.sub(r"^\[HMM\] ", "", homolog)
    return f"{base} ({motif})" if motif else base


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # ---- 1. collect one record per MTase gene across all 15 MAGs ----
    genes = pd.concat([load_mtases(m) for m in MAGS], ignore_index=True)
    genes["system_type"] = genes["System Type"].apply(clean_systype)
    detected_by_mag = {m: load_detected_motifs(m) for m in MAGS}

    # ---- 2. cluster the protein sequences into ortholog groups ----
    gene2rep = cluster_orthologs(set(genes["Gene"]))
    genes["cluster_rep"] = genes["Gene"].map(gene2rep)
    if genes["cluster_rep"].isna().any():
        # any gene without a sequence becomes its own singleton cluster
        genes["cluster_rep"] = genes["cluster_rep"].fillna(genes["Gene"])

    # ---- 3. assign each cluster a stable key + consensus annotation ----
    cl_info = {}
    for rep, grp in genes.groupby("cluster_rep"):
        cl_info[rep] = {
            "label": cluster_label(grp),
            "rebase_homolog": most_common_nonempty(grp["REBASE homolog"]),
            "hmm": most_common_nonempty(grp["HMM"]),
            "motif": first_nonempty(grp["Homolog motif"]),
            "mod_type": first_nonempty(grp["Predicted methylation"]),
            "system_type": most_common_nonempty(grp["system_type"]),
            "rebase_names": ";".join(sorted(set(
                v for v in grp["REBASE homolog"] if isinstance(v, str) and v))),
            "n_genes": len(grp),
            "genes": ";".join(sorted(grp["Gene"])),
        }
    # de-duplicate labels (append motif/hmm suffix if two clusters share a label)
    seen = {}
    for rep in cl_info:
        lab = cl_info[rep]["label"]
        if lab in seen.values():
            lab = f"{lab} #{sum(1 for v in seen.values() if v.startswith(cl_info[rep]['label'])) + 1}"
        seen[rep] = lab
        cl_info[rep]["label"] = lab

    # order columns: by system type then label, for a tidy default
    reps = sorted(cl_info, key=lambda r: (cl_info[r]["system_type"], cl_info[r]["label"]))
    keys = [cl_info[r]["label"] for r in reps]
    rep_of_key = {cl_info[r]["label"]: r for r in reps}

    # ---- 4. build matrices (rows = MAGs, cols = cluster labels) ----
    presence = {m: {} for m in MAGS}
    copies = {m: {} for m in MAGS}
    activity = {m: {} for m in MAGS}
    catalog_rows = []

    for rep, grp in genes.groupby("cluster_rep"):
        key = cl_info[rep]["label"]
        motif = cl_info[rep]["motif"]
        for mag, mag_grp in grp.groupby("MAG"):
            n = len(mag_grp)
            presence[mag][key] = 1
            copies[mag][key] = n
            if motif:
                hits = sorted(
                    {f"{dm}({frac:.2f})" for dm, frac in detected_by_mag[mag]
                     if frac >= ACTIVITY_FRACTION and motifs_equivalent(motif, dm)}
                )
                active = int(len(hits) > 0)
                matched = ";".join(hits)
            else:
                active = None
                matched = ""
            activity[mag][key] = active
            for _, r in mag_grp.iterrows():
                catalog_rows.append({
                    "MAG": mag, "cluster": key, "cluster_rep": rep,
                    "gene": r["Gene"], "operon": r["Operon"],
                    "rebase_homolog": r["REBASE homolog"], "hmm": r["HMM"],
                    "recognition_motif": motif, "mod_type": cl_info[rep]["mod_type"],
                    "system_type": r["system_type"],
                    "homolog_identity": r.get("Homolog identity(%)", ""),
                    "copies_in_mag": n, "active": active,
                    "matched_detected_motif": matched,
                })

    def as_matrix(d, fill=0):
        return pd.DataFrame([[d[m].get(k, fill) for k in keys] for m in MAGS],
                            index=MAGS, columns=keys)

    pres_df = as_matrix(presence, 0).astype(int)
    copy_df = as_matrix(copies, 0).astype(int)
    act_df = pd.DataFrame([[activity[m].get(k, None) for k in keys] for m in MAGS],
                          index=MAGS, columns=keys)

    annot_df = pd.DataFrame([{
        "cluster": cl_info[r]["label"],
        "rebase_homolog": cl_info[r]["rebase_homolog"],
        "hmm": cl_info[r]["hmm"], "motif": cl_info[r]["motif"],
        "mod_type": cl_info[r]["mod_type"], "system_type": cl_info[r]["system_type"],
    } for r in reps]).set_index("cluster").loc[keys]

    members_df = pd.DataFrame([{
        "cluster": cl_info[r]["label"], "cluster_rep": r,
        "n_genes": cl_info[r]["n_genes"], "rebase_names": cl_info[r]["rebase_names"],
        "genes": cl_info[r]["genes"],
    } for r in reps])

    catalog_df = pd.DataFrame(catalog_rows)

    # ---- write source data ----
    pres_df.to_csv(os.path.join(OUT_DIR, "mtase_presence_absence.csv"))
    copy_df.to_csv(os.path.join(OUT_DIR, "mtase_copy_number.csv"))
    act_df.to_csv(os.path.join(OUT_DIR, "mtase_activity.csv"))
    annot_df.to_csv(os.path.join(OUT_DIR, "mtase_column_annotation.csv"))
    members_df.to_csv(os.path.join(OUT_DIR, "mtase_cluster_members.csv"), index=False)
    catalog_df.to_csv(os.path.join(OUT_DIR, "mtase_catalog.csv"), index=False)

    # ---- screen summary ----
    print(f"\nKp 248_1 MTase repertoire: {len(MAGS)} MAGs x {len(keys)} ortholog clusters "
          f"(mmseqs {MIN_SEQ_ID} id / {MIN_COV} cov)\n")
    print("Per-MAG MTase gene count (MT/IIG):")
    for m in MAGS:
        print(f"  {m:16s} {int(copy_df.loc[m].sum()):3d}")
    print("\nPer-cluster presence / active counts across the 15 MAGs:")
    print(f"  {'cluster (homolog + motif)':40s} {'mod':5s} {'sys':9s} "
          f"{'#genes':>6s} {'pres':>4s} {'act':>4s}  rebase_names")
    for r in reps:
        k = cl_info[r]["label"]
        n_pres = int(pres_df[k].sum())
        n_act = int(act_df[k].fillna(0).sum())
        print(f"  {k[:40]:40s} {cl_info[r]['mod_type']:5s} "
              f"{str(cl_info[r]['system_type']):9s} {cl_info[r]['n_genes']:6d} "
              f"{n_pres:4d} {n_act:4d}  {cl_info[r]['rebase_names']}")

    print(f"\nWrote source data to {OUT_DIR}")
    return pres_df, act_df, annot_df


if __name__ == "__main__":
    main()
