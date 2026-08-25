#!/usr/bin/env python
"""
Dereplicate MTase genes per genome (contig) with mmseqs2 and redraw the
MTase-vs-motif correlation for the metagenome high-quality contigs.

Motivation
----------
The MTase count that goes into the MTase-vs-motif correlation has been based on
HMM family groups (distinct `HMM` labels among `MT`-type genes per contig). HMM
grouping can misrepresent the true number of distinct MTase enzymes in a genome
(near-identical gene copies split across families, or several distinct enzymes
collapsed into one family). Here we instead dereplicate the MTase protein
sequences of each genome with mmseqs2 easy-cluster (>=90% amino-acid identity,
>=80% coverage) and count the resulting clusters as the number of unique MTases.

Per-genome dereplication
------------------------
Clustering is run independently for each contig (genome): two MTase copies on the
same genome are collapsed only when they themselves meet the 90/80 threshold. A
contig with 0 or 1 MTase gene needs no clustering (unique count = 0 or 1).

We reuse the existing per-contig motif counts (motif_num_all_samples.csv), which
already restrict to the high-quality contigs (present in 59 of the 64
metagenomes), and only recompute the MTase count. We also carry the old
HMM-group count (MT_hmm_num) and the RM-system count (RM_num) for comparison.

Outputs
-------
  data  : /home/shuaiw/borg/revision/MTase_drep/mtase_gene_clusters.tsv
  source: /home/shuaiw/MODIFI/tmp/rev_figs/MTase_drep/mtase_drep_motif_corr.sourcedata.csv
  figure: /home/shuaiw/MODIFI/tmp/rev_figs/MTase_drep/MTase_drep_vs_motif_num.{pdf,png}
"""

import os
import shutil
import subprocess
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# ----------------------------------------------------------------------------
# Paths / configuration
# ----------------------------------------------------------------------------
RUN2 = "/home/shuaiw/borg/paper/run2"
METHY = "methylation4"                       # canonical run
MOTIF_CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/motif_num_all_samples.csv"

DATA_DIR = "/home/shuaiw/borg/revision/MTase_drep"      # off the small /home SSD
# per-contig mmseqs scratch on fast node-local disk (transient; consolidated TSV -> DATA_DIR).
# Override with MTASE_DREP_WORK (the SLURM wrapper points this at node-local /tmp).
WORK_DIR = os.environ.get("MTASE_DREP_WORK", os.path.join(DATA_DIR, "work"))
OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/MTase_drep"

MMSEQS = "/shared/software/bin/mmseqs"
# amino-acid identity / coverage thresholds; override for a cutoff sweep
MIN_SEQ_ID = float(os.environ.get("MTASE_MIN_SEQ_ID", "0.9"))
MIN_COV = float(os.environ.get("MTASE_MIN_COV", "0.8"))
TAG = f"id{int(round(MIN_SEQ_ID * 100))}_cov{int(round(MIN_COV * 100))}"
# per-contig mmseqs pool size; full node under SLURM (each mmseqs job is --threads 1)
N_WORKERS = int(os.environ.get("SLURM_CPUS_ON_NODE", "24"))

MTASE_TYPES = {"MT", "IIG"}                  # methyltransferase gene types


# ----------------------------------------------------------------------------
# Parsing helpers (idiom from benchmark/Klevsuekka_mtase/build_mtase_matrix.py)
# ----------------------------------------------------------------------------
def rm_table_path(sample):
    return os.path.join(RUN2, sample, f"{sample}_{METHY}", "RM_systems",
                        "all_ctgs_RM.rm.genes.tsv")


def faa_path(sample):
    return os.path.join(RUN2, sample, f"{sample}_{METHY}", "RM_systems",
                        "all_ctgs_RM.rm.genes.faa")


def contig_of(gene):
    """Contig id = gene id minus its last '_'-field (e.g. soil_1_3_C_149 -> soil_1_3_C)."""
    return "_".join(str(gene).split("_")[:-1])


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


def load_sample_mtases(sample):
    """
    Return, for one metagenome sample:
      contig_genes : {contig: [gene_id, ...]}   (MT/IIG genes, deduped by Operon)
      contig_hmms  : {contig: set(HMM labels)}  (distinct HMM among MT genes; old count)
    """
    tsv = rm_table_path(sample)
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    df = df[df["Gene type"].isin(MTASE_TYPES)].copy()
    # one MTase gene can produce several domain hits sharing an Operon; keep one
    df = df.drop_duplicates(subset=["Operon"])
    df["contig"] = df["Gene"].apply(contig_of)

    contig_genes = defaultdict(list)
    contig_hmms = defaultdict(set)
    for _, r in df.iterrows():
        contig_genes[r["contig"]].append(r["Gene"])
        if isinstance(r["HMM"], str) and r["HMM"].strip():
            contig_hmms[r["contig"]].add(r["HMM"].strip())
    return contig_genes, contig_hmms


# ----------------------------------------------------------------------------
# Per-contig mmseqs dereplication
# ----------------------------------------------------------------------------
def cluster_contig(args):
    """
    Cluster one contig's MTase proteins with mmseqs easy-cluster.
    args = (contig, {gene_id: sequence})
    Returns (contig, n_clusters, [(gene_id, cluster_rep), ...]).
    Contigs with <=1 sequence skip mmseqs.
    """
    contig, seqs = args
    genes = sorted(seqs)
    if len(genes) <= 1:
        rep = genes[0] if genes else None
        pairs = [(g, g) for g in genes]
        return contig, len(genes), pairs

    safe = contig.replace("/", "_")
    fasta = os.path.join(WORK_DIR, f"{safe}.faa")
    prefix = os.path.join(WORK_DIR, f"{safe}_clu")
    tmp = os.path.join(WORK_DIR, f"{safe}_tmp")
    with open(fasta, "w") as out:
        for g in genes:
            out.write(f">{g}\n{seqs[g]}\n")

    cmd = [MMSEQS, "easy-cluster", fasta, prefix, tmp,
           "--min-seq-id", str(MIN_SEQ_ID), "-c", str(MIN_COV),
           "--cov-mode", "0", "--threads", "1", "-v", "1"]
    subprocess.run(cmd, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    pairs = []
    reps = set()
    with open(prefix + "_cluster.tsv") as f:
        for line in f:
            rep, member = line.rstrip("\n").split("\t")
            pairs.append((member, rep))
            reps.add(rep)

    # clean up mmseqs scratch, keep nothing per-contig
    for p in (fasta, prefix + "_cluster.tsv", prefix + "_rep_seq.fasta",
              prefix + "_all_seqs.fasta"):
        if os.path.exists(p):
            os.remove(p)
    shutil.rmtree(tmp, ignore_errors=True)
    return contig, len(reps), pairs


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(WORK_DIR, exist_ok=True)

    motif = pd.read_csv(MOTIF_CSV)
    samples = sorted(motif["sample"].unique())
    # only score the high-quality contigs that are actually in the correlation
    keep_contigs = set(motif["contig"])
    print(f"Loaded {len(motif)} high-quality contigs across {len(samples)} metagenomes")

    # ---- 1. gather MTase genes + sequences per contig (high-quality contigs only) ----
    contig_genes = {}          # contig -> [gene, ...]
    contig_hmms = {}           # contig -> set(HMM)
    contig_seqs = {}           # contig -> {gene: seq}
    for s in samples:
        cg, ch = load_sample_mtases(s)
        cg = {c: genes for c, genes in cg.items() if c in keep_contigs}
        if not cg:
            continue
        need = {g for genes in cg.values() for g in genes}
        faa = read_faa(faa_path(s))
        for contig, genes in cg.items():
            contig_genes[contig] = genes
            contig_hmms[contig] = ch.get(contig, set())
            contig_seqs[contig] = {g: faa[g] for g in genes if g in faa}
        missing = [g for g in need if g not in faa]
        if missing:
            print(f"  {s}: WARNING {len(missing)} MTase gene(s) had no protein sequence")
    print(f"Contigs carrying >=1 MTase gene: {len(contig_genes)}")
    multi = sum(1 for c in contig_seqs if len(contig_seqs[c]) >= 2)
    print(f"Contigs needing mmseqs (>=2 MTase genes): {multi}")

    # ---- 2. per-genome (per-contig) mmseqs dereplication ----
    jobs = [(c, contig_seqs[c]) for c in contig_seqs]
    drep_num = {}
    cluster_rows = []
    with Pool(N_WORKERS) as pool:
        for contig, n_clusters, pairs in pool.imap_unordered(cluster_contig, jobs, chunksize=8):
            drep_num[contig] = n_clusters
            sample = "_".join(contig.split("_")[:-2]) if contig.count("_") >= 2 else contig
            for gene, rep in pairs:
                cluster_rows.append({"sample": sample, "contig": contig,
                                     "gene": gene, "cluster_rep": rep})
    # contigs present in the motif table but with no MTase gene at all -> 0
    for contig in motif["contig"]:
        drep_num.setdefault(contig, 0)

    clusters_df = pd.DataFrame(cluster_rows).sort_values(["sample", "contig", "gene"])
    clusters_df.to_csv(os.path.join(DATA_DIR, f"mtase_gene_clusters.{TAG}.tsv"),
                       sep="\t", index=False)

    # ---- 3. assemble the per-contig source table (join onto motif counts) ----
    df = motif[["sample", "environment", "contig", "motif_num", "RM_num",
                "ctg_len"]].copy()
    df["MT_hmm_num"] = df["contig"].map(lambda c: len(contig_hmms.get(c, set())))
    df["MT_drep_num"] = df["contig"].map(lambda c: drep_num.get(c, 0)).astype(int)
    df = df[["sample", "environment", "contig", "motif_num", "RM_num",
             "MT_hmm_num", "MT_drep_num", "ctg_len"]]
    src = os.path.join(OUT_DIR, f"mtase_drep_motif_corr.{TAG}.sourcedata.csv")
    df.to_csv(src, index=False)

    # ---- 4. correlation + figure ----
    r, p = pearsonr(df["MT_drep_num"], df["motif_num"])
    r_hmm, p_hmm = pearsonr(df["MT_hmm_num"], df["motif_num"])
    r_rm, p_rm = pearsonr(df["RM_num"], df["motif_num"])

    envs = sorted(df["environment"].dropna().unique())
    palette = dict(zip(envs, sns.color_palette("tab10", len(envs))))

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    rng = np.random.default_rng(0)
    jx = df["MT_drep_num"] + rng.uniform(-0.18, 0.18, len(df))
    jy = df["motif_num"] + rng.uniform(-0.18, 0.18, len(df))
    for env in envs:
        m = df["environment"] == env
        ax.scatter(jx[m], jy[m], s=16, alpha=0.65, color=palette[env],
                   edgecolor="none", label=env)
    sns.regplot(x="MT_drep_num", y="motif_num", data=df, scatter=False,
                ax=ax, color="black", line_kws={"linewidth": 1.4},
                truncate=False)
    ax.set_xlabel(f"No. of unique MTases (mmseqs, >={int(round(MIN_SEQ_ID*100))}% id, "
                  f">={int(round(MIN_COV*100))}% cov)")
    ax.set_ylabel("No. of motifs")
    ax.set_title(f"Pearson r = {r:.2f}, p = {p:.2e}  (n = {len(df)})")
    ax.legend(title="Environment", fontsize=7, title_fontsize=8,
              markerscale=1.4, frameon=False, loc="best")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"MTase_drep_vs_motif_num.{TAG}.pdf"))
    fig.savefig(os.path.join(OUT_DIR, f"MTase_drep_vs_motif_num.{TAG}.png"), dpi=300)

    # ---- 5. screen summary (for copy-paste) ----
    print(f"\n=== cutoff {TAG} (>= {MIN_SEQ_ID} id, >= {MIN_COV} cov) ===")
    print("=== Correlation with No. of motifs (n = %d contigs) ===" % len(df))
    print(f"  mmseqs unique MTase (MT_drep_num): Pearson r = {r:.3f}, p = {p:.3e}")
    print(f"  HMM-group MTase     (MT_hmm_num) : Pearson r = {r_hmm:.3f}, p = {p_hmm:.3e}")
    print(f"  RM systems          (RM_num)     : Pearson r = {r_rm:.3f}, p = {p_rm:.3e}")
    delta = (df["MT_hmm_num"] - df["MT_drep_num"])
    print(f"\nMT_hmm_num vs MT_drep_num: mean {df['MT_hmm_num'].mean():.2f} vs "
          f"{df['MT_drep_num'].mean():.2f}; "
          f"contigs where drep < hmm: {int((delta > 0).sum())}, "
          f"drep > hmm: {int((delta < 0).sum())}, equal: {int((delta == 0).sum())}")
    print("\nSource data head:")
    print(df.head(10).to_string(index=False))
    print(f"\nWrote:\n  {src}"
          f"\n  {os.path.join(OUT_DIR, f'MTase_drep_vs_motif_num.{TAG}.pdf')}"
          f"\n  {os.path.join(DATA_DIR, f'mtase_gene_clusters.{TAG}.tsv')}")


if __name__ == "__main__":
    main()
