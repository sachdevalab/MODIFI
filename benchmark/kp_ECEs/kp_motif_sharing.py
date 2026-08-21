#!/usr/bin/env python3
"""
Step 5 - Does Klebsiella just dominate, or do its ECEs genuinely fail to share modification motifs
with co-resident Enterobacteriaceae (so they cannot transfer)?

Two complementary read-outs:

(A) Network-based host reach of the 22 Klebsiella-associated ECE clusters: how many are
    Klebsiella-exclusive vs. shared across other Enterobacteriaceae genera (from kp_22_clusters.tsv).

(B) Motif-based: in each infant gut that co-harbours a Klebsiella host chromosome and >=1 other
    Enterobacteriaceae host chromosome, quantify how similar their modification-motif profiles are
    (cosine on the full per-contig methylation-fraction vector, matching the pipeline's cos_sim;
    plus Jaccard on the binarised "modified motif set"). Same computation for each Klebsiella-linked
    ECE vs. the co-resident non-Klebsiella Enterobacteriaceae host. Low similarity => motif-blocked
    (supports host restriction of transfer); high similarity => motif-permissive (transfer expected,
    and indeed the shared multi-genus clusters realise it).

Inputs (all read-only):
  - network linkage table (host contigs + host_taxa)
  - per-sample cross-sample motif profiles (contig x motif methylation fraction)
  - kp_22_clusters.tsv (from build_kp_clusters.py)
Outputs to /home/shuaiw/borg/revision/kp_eces and a copy of source data to the figure dir.
"""
import os
import re
import glob
import numpy as np
import pandas as pd

LINKAGE_CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
PROFILE_GLOB = "/groups/banfield/projects/multienv/methylation/infant_cross_sample_meth/{sample}_motif_profile_cross.csv"
CLUSTERS_TSV = "/home/shuaiw/borg/revision/kp_eces/kp_22_clusters.tsv"
OUT_DIR = "/home/shuaiw/borg/revision/kp_eces"
FIG_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"

ENTERO_GENERA = {"Klebsiella", "Citrobacter", "Enterobacter", "Escherichia", "Serratia",
                 "Salmonella", "Pluralibacter", "Leclercia", "Cronobacter", "Shigella",
                 "Kluyvera", "Raoultella", "Kosakonia"}
MOD_THRESHOLD = 0.5  # a motif is "modified" on a contig if methylation fraction >= this


def species_of(taxa):
    return re.sub(r"^s__", "", str(taxa)).strip()


def genus_of(taxa):
    sp = species_of(taxa)
    if not sp:
        return "NA"
    return re.sub(r"_[A-Z]$", "", sp.split(" ")[0])


def cosine(a, b):
    a = np.nan_to_num(a.astype(float)); b = np.nan_to_num(b.astype(float))
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return np.nan
    return float(np.dot(a, b) / (na * nb))


def jaccard(a, b, thr=MOD_THRESHOLD):
    sa = set(np.where(np.nan_to_num(a.astype(float)) >= thr)[0])
    sb = set(np.where(np.nan_to_num(b.astype(float)) >= thr)[0])
    if not sa and not sb:
        return np.nan
    return len(sa & sb) / len(sa | sb)


def load_profile(sample):
    path = PROFILE_GLOB.format(sample=sample)
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, index_col=0)  # rows=motifs, cols=contigs
    return df


def main():
    gc = pd.read_csv(LINKAGE_CSV)
    gc["host_genus"] = gc["host_taxa"].map(genus_of)
    gc["host_species"] = gc["host_taxa"].map(species_of)
    gc["ece_genus"] = gc["host_genus"]  # the ECE's linked host genus

    # ---- (A) network host reach of the 22 clusters ----
    clu = pd.read_csv(CLUSTERS_TSV, sep="\t")
    clu["is_multi_genus"] = clu["n_distinct_host_genera"] > 1
    n_excl = int((~clu["is_multi_genus"]).sum())
    n_multi = int(clu["is_multi_genus"].sum())
    multi = clu[clu["is_multi_genus"]][
        ["cluster_rep", "n_cluster_members", "n_distinct_host_genera", "linked_host_genera",
         "linked_host_species"]]
    with open(os.path.join(OUT_DIR, "kp_cluster_hostreach_summary.txt"), "w") as fh:
        fh.write(f"22 Klebsiella pneumoniae-associated ECE clusters:\n")
        fh.write(f"  Klebsiella-exclusive (single host genus in network): {n_excl}\n")
        fh.write(f"  span >1 host genus (cross-Enterobacteriaceae): {n_multi}\n\n")
        fh.write("Multi-genus clusters:\n")
        fh.write(multi.to_string(index=False))
        fh.write("\n")
    print(f"[A] host reach: {n_excl} Klebsiella-exclusive, {n_multi} multi-genus of 22")
    print(multi.to_string(index=False))

    # ---- (B) motif-profile similarity among co-resident Enterobacteriaceae hosts ----
    # per infant sample, the set of Enterobacteriaceae host chromosomes in the network
    infant = gc[gc["sample"].str.startswith("infant") | gc["sample"].str.startswith("asthma")]
    entero_hosts = (infant[infant["host_genus"].isin(ENTERO_GENERA)]
                    [["sample", "host", "host_species", "host_genus"]]
                    .drop_duplicates())

    host_rows = []       # host-vs-host similarity
    ece_rows = []        # klebsiella ECE-vs-coresident host similarity
    for sample, sub in entero_hosts.groupby("sample"):
        genera = set(sub["host_genus"])
        if "Klebsiella" not in genera or len(sub) < 2:
            continue
        prof = load_profile(sample)
        if prof is None:
            continue
        cols = set(prof.columns)
        kleb = sub[sub["host_genus"] == "Klebsiella"]
        others = sub[sub["host_genus"] != "Klebsiella"]
        for _, k in kleb.iterrows():
            if k["host"] not in cols:
                continue
            kvec = prof[k["host"]].values
            for _, o in others.iterrows():
                if o["host"] not in cols:
                    continue
                ovec = prof[o["host"]].values
                host_rows.append({
                    "sample": sample,
                    "kleb_host": k["host"], "kleb_species": k["host_species"],
                    "other_host": o["host"], "other_species": o["host_species"],
                    "other_genus": o["host_genus"],
                    "cosine": round(cosine(kvec, ovec), 4),
                    "jaccard_modset": round(jaccard(kvec, ovec), 4),
                })
        # ECE-vs-coresident-host: Klebsiella-linked ECEs in this sample
        sample_eces = gc[(gc["sample"] == sample) &
                         (gc["host_species"] == "Klebsiella pneumoniae")]["MGE"].unique()
        for ece in sample_eces:
            if ece not in cols:
                continue
            evec = prof[ece].values
            for _, o in others.iterrows():
                if o["host"] not in cols:
                    continue
                ece_rows.append({
                    "sample": sample, "ece": ece,
                    "other_host": o["host"], "other_species": o["host_species"],
                    "other_genus": o["host_genus"],
                    "cosine": round(cosine(evec, prof[o["host"]].values), 4),
                    "jaccard_modset": round(jaccard(evec, prof[o["host"]].values), 4),
                })

    host_df = pd.DataFrame(host_rows)
    ece_df = pd.DataFrame(ece_rows)
    host_df.to_csv(os.path.join(OUT_DIR, "kp_host_motif_similarity.tsv"), sep="\t", index=False)
    ece_df.to_csv(os.path.join(OUT_DIR, "kp_ece_host_motif_similarity.tsv"), sep="\t", index=False)
    os.makedirs(FIG_DIR, exist_ok=True)
    host_df.to_csv(os.path.join(FIG_DIR, "kp_host_motif_similarity_sourcedata.csv"), index=False)

    print(f"\n[B] host-host Enterobacteriaceae motif comparisons: {len(host_df)} pairs "
          f"across {host_df['sample'].nunique() if len(host_df) else 0} guts")
    if len(host_df):
        print(host_df.to_string(index=False))
        print(f"\n  Klebsiella-vs-other-Enterobacteriaceae cosine: "
              f"median={host_df['cosine'].median():.3f}, "
              f"range={host_df['cosine'].min():.3f}-{host_df['cosine'].max():.3f}")
        print(f"  Klebsiella-vs-other-Enterobacteriaceae Jaccard(modset): "
              f"median={host_df['jaccard_modset'].median():.3f}")
    if len(ece_df):
        print(f"\n  Klebsiella-ECE vs co-resident non-Klebsiella host cosine: "
              f"median={ece_df['cosine'].median():.3f}, n={len(ece_df)}")


if __name__ == "__main__":
    main()
