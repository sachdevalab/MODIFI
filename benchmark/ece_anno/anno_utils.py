#!/usr/bin/env python3
"""
Small annotation helpers for the ECE-evidence module.

Trimmed, pandas-based reimplementation of the patterns in
/home/mingy/MING/helper_scripts/contigs_anno.py (which uses polars). ECE protein
sets are small (hundreds to a few thousand ORFs), so a plain `hmmsearch --cpu N`
is used instead of the seqkit-split + GNU-parallel fan-out.
"""
import os
import subprocess

import pandas as pd


def run_hmmsearch(input_faa, hmm_db, out_tblout, threads=8,
                  cut_ga=False, cut_tc=False, evalue=1e-5, reuse=False):
    """Run hmmsearch of `input_faa` against a pressed `hmm_db`, writing --tblout.

    Uses profile gathering/trusted cutoffs when requested, else an E-value cutoff.
    Returns out_tblout. If input has no sequences, writes an empty tblout and returns.
    If `reuse` and a non-empty tblout already exists, it is kept (skips the search).
    """
    if reuse and os.path.exists(out_tblout) and os.path.getsize(out_tblout) > 0:
        return out_tblout
    if not os.path.exists(input_faa) or os.path.getsize(input_faa) == 0:
        open(out_tblout, "w").close()
        return out_tblout
    cmd = ["hmmsearch", "--cpu", str(threads), "--tblout", out_tblout]
    if cut_ga:
        cmd.append("--cut_ga")
    elif cut_tc:
        cmd.append("--cut_tc")
    else:
        cmd += ["-E", str(evalue)]
    cmd += [hmm_db, input_faa]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return out_tblout


def parse_hmmer_tblout(tblout_path):
    """Parse a HMMER --tblout into a DataFrame.

    Columns: gene_id, query_name, query_acc, evalue, score.
    query_acc has its Pfam version suffix stripped (PF03354.14 -> PF03354).
    """
    rows = []
    if os.path.exists(tblout_path):
        with open(tblout_path) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.split()
                if len(p) < 6:
                    continue
                rows.append({
                    "gene_id": p[0],
                    "query_name": p[2],
                    "query_acc": p[3].split(".")[0],
                    "evalue": float(p[4]),
                    "score": float(p[5]),
                })
    return pd.DataFrame(
        rows, columns=["gene_id", "query_name", "query_acc", "evalue", "score"]
    )


def best_hit_per_gene(df, score_col="score"):
    """Keep the single highest-scoring row per gene_id."""
    if df.empty:
        return df
    return df.sort_values(score_col, ascending=False).drop_duplicates("gene_id")


def gene_to_contig(gene_id):
    """geNomad/pyrodigal protein id -> parent contig (strip the trailing _<geneidx>)."""
    return gene_id.rsplit("_", 1)[0]


def run_pyrodigal(input_fna, out_faa, threads=8):
    """Fallback ORF calling when geNomad proteins are unavailable.

    pyrodigal-gv meta mode; sanitizes '#' in headers for macsyfinder gembase.
    """
    gff = out_faa + ".gff"
    fna = out_faa + ".fna"
    tmp = out_faa + ".raw"
    subprocess.run(
        ["pyrodigal-gv", "-i", input_fna, "-a", tmp, "-d", fna, "-o", gff,
         "-p", "meta", "-j", str(threads)],
        check=True,
    )
    with open(tmp) as fin, open(out_faa, "w") as fout:
        for line in fin:
            fout.write(line.replace("#", "_") if line.startswith(">") else line)
    os.remove(tmp)
    return out_faa
