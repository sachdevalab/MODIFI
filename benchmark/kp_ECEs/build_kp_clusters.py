#!/usr/bin/env python3
"""
Build the roster of ECE (plasmid/virus) clusters that associate with a target host species
(default: Klebsiella pneumoniae), reproducing the manuscript "22 clusters" number.

Logic mirrors benchmark/network/plot_linkage_data.py::analyze_MGEs (lines 614-622):
  - read the host-ECE linkage table (mge_host_gc_cov.csv)
  - filter rows whose host_taxa contains the target species
  - map each linked ECE (MGE) to its 95% ANI cluster (megablast.cluster.95ani.tsv)
  - the unique cluster set is the answer

Also emits, for every one of those clusters, the FULL membership across all samples with per-member
host species/genus (from the linkage table), so cross-genus clusters (e.g. infant_15_35_C) are flagged.

Read-only on inputs. Writes TSVs to the data dir and copies a summary to the figure dir.
"""
import os
import sys
import re
import argparse
from collections import defaultdict
import pandas as pd

LINKAGE_CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
CLUSTER_TSV = "/home/shuaiw/borg/paper/MGE/cluster/megablast.cluster.95ani.tsv"
OUT_DIR = "/home/shuaiw/borg/revision/kp_eces"
FIG_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"


def read_mge_cluster(mge_clu_file):
    """Reuse of plot_linkage_data.read_mge_cluster: rep<TAB>comma-joined members."""
    mge_clu_dict = {}
    clu_mge_dict = defaultdict(list)
    df = pd.read_csv(mge_clu_file, sep="\t", header=None)
    for _, row in df.iterrows():
        cluster = row[0]
        mges = str(row[1]).split(",")
        for mge in mges:
            mge_clu_dict[mge] = cluster
            clu_mge_dict[cluster].append(mge)
    return mge_clu_dict, clu_mge_dict


def species_from_taxa(host_taxa):
    """host_taxa is like 's__Klebsiella pneumoniae' -> 'Klebsiella pneumoniae'."""
    if not isinstance(host_taxa, str):
        return "NA"
    return re.sub(r"^s__", "", host_taxa).strip()


def genus_from_taxa(host_taxa):
    sp = species_from_taxa(host_taxa)
    if sp in ("NA", ""):
        return "NA"
    # GTDB alphabetic suffixes like 'Citrobacter_A' collapse to the base genus
    genus = sp.split(" ")[0]
    return re.sub(r"_[A-Z]$", "", genus)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxon", default="s__Klebsiella pneumoniae",
                    help="host_taxa substring to select (default: s__Klebsiella pneumoniae)")
    ap.add_argument("--linkage", default=LINKAGE_CSV)
    ap.add_argument("--clusters", default=CLUSTER_TSV)
    ap.add_argument("--outdir", default=OUT_DIR)
    ap.add_argument("--figdir", default=FIG_DIR)
    args = ap.parse_args()

    gc_df = pd.read_csv(args.linkage)
    mge_clu_dict, clu_mge_dict = read_mge_cluster(args.clusters)

    # --- reproduce analyze_MGEs ---
    target = gc_df[gc_df["host_taxa"].astype(str).str.contains(args.taxon, regex=False)].copy()
    target["MGE_cluster"] = target["MGE"].map(mge_clu_dict)
    unique_clusters = sorted(target["MGE_cluster"].dropna().unique().tolist())
    n_member_eces = target["MGE"].nunique()

    print(f"[build_kp_clusters] target taxon: {args.taxon}")
    print(f"[build_kp_clusters] member ECE sequences linked to target: {n_member_eces}")
    print(f"[build_kp_clusters] unique 95%-ANI ECE clusters: {len(unique_clusters)}")
    missing = target[target["MGE_cluster"].isna()]["MGE"].tolist()
    if missing:
        print(f"[build_kp_clusters] WARNING: {len(missing)} linked ECEs not found in cluster file: {missing}")

    # --- per-member linkage annotation for the target-linked ECEs ---
    keep_cols = ["MGE", "host", "MGE_type", "environment", "sample", "mge_len",
                 "host_taxa", "host_phylum", "cos_sim", "MGE_cov", "host_cov"]
    kp_linked = target[keep_cols].copy()
    kp_linked["MGE_cluster"] = target["MGE_cluster"]
    kp_linked["host_species"] = kp_linked["host_taxa"].map(species_from_taxa)
    kp_linked["host_genus"] = kp_linked["host_taxa"].map(genus_from_taxa)
    kp_linked = kp_linked.sort_values(["MGE_cluster", "sample"])
    kp_linked_path = os.path.join(args.outdir, "kp_linked_eces.tsv")
    kp_linked.to_csv(kp_linked_path, sep="\t", index=False)
    print(f"[build_kp_clusters] wrote {kp_linked_path} ({len(kp_linked)} target-linked ECE rows)")

    # --- full cluster membership + which host genera/species each cluster links to (network-wide) ---
    # For every ECE in the whole linkage table, we know its host species; use that to describe the
    # cross-genus reach of each of the target clusters.
    all_link = gc_df.copy()
    all_link["MGE_cluster"] = all_link["MGE"].map(mge_clu_dict)
    all_link["host_species"] = all_link["host_taxa"].map(species_from_taxa)
    all_link["host_genus"] = all_link["host_taxa"].map(genus_from_taxa)

    rows = []
    for clu in unique_clusters:
        members = clu_mge_dict.get(clu, [])
        # linkage rows anywhere in the network for members of this cluster
        clu_links = all_link[all_link["MGE_cluster"] == clu]
        linked_species = sorted(set(clu_links["host_species"].dropna()) - {"NA", ""})
        linked_genera = sorted(set(clu_links["host_genus"].dropna()) - {"NA", ""})
        # representative info (the cluster id is the representative seq name)
        rep = clu
        rep_rows = all_link[all_link["MGE"] == rep]
        rep_len = int(rep_rows["mge_len"].iloc[0]) if len(rep_rows) else ""
        rep_type = rep_rows["MGE_type"].iloc[0] if len(rep_rows) else ""
        rep_host = rep_rows["host_species"].iloc[0] if len(rep_rows) else ""
        n_samples = clu_links["sample"].nunique()
        rows.append({
            "cluster_rep": clu,
            "rep_type": rep_type,
            "rep_len_bp": rep_len,
            "rep_linked_host": rep_host,
            "n_cluster_members": len(members),
            "n_members_linked_in_network": clu_links["MGE"].nunique(),
            "n_samples_linked": n_samples,
            "n_distinct_host_genera": len(linked_genera),
            "n_distinct_host_species": len(linked_species),
            "linked_host_genera": ";".join(linked_genera),
            "linked_host_species": ";".join(linked_species),
            "all_cluster_members": ";".join(members),
        })
    clu_df = pd.DataFrame(rows).sort_values("n_cluster_members", ascending=False)
    clu_path = os.path.join(args.outdir, "kp_22_clusters.tsv")
    clu_df.to_csv(clu_path, sep="\t", index=False)
    print(f"[build_kp_clusters] wrote {clu_path} ({len(clu_df)} clusters)")

    # compact summary copy for the figure/response dir
    os.makedirs(args.figdir, exist_ok=True)
    summary_path = os.path.join(args.figdir, "kp_22_clusters_summary.tsv")
    clu_df.to_csv(summary_path, sep="\t", index=False)

    # cross-genus flag summary
    cross = clu_df[clu_df["n_distinct_host_genera"] > 1]
    print(f"[build_kp_clusters] clusters linking >1 host genus: {len(cross)} / {len(clu_df)}")
    if "infant_15_35_C" in set(clu_df["cluster_rep"]):
        r = clu_df[clu_df["cluster_rep"] == "infant_15_35_C"].iloc[0]
        print(f"[build_kp_clusters] infant_15_35_C -> genera=[{r['linked_host_genera']}] "
              f"species=[{r['linked_host_species']}] members={r['n_cluster_members']}")

    # write the plain member-ECE list for sequence extraction downstream
    all_member_eces = sorted(set(
        [m for clu in unique_clusters for m in clu_mge_dict.get(clu, [])]
    ))
    kp_linked_eces = sorted(set(kp_linked["MGE"]))
    with open(os.path.join(args.outdir, "kp_all_cluster_member_eces.txt"), "w") as fh:
        fh.write("\n".join(all_member_eces) + "\n")
    with open(os.path.join(args.outdir, "kp_linked_eces.txt"), "w") as fh:
        fh.write("\n".join(kp_linked_eces) + "\n")
    print(f"[build_kp_clusters] {len(all_member_eces)} total cluster-member ECEs, "
          f"{len(kp_linked_eces)} directly Klebsiella-linked ECEs (for sequence extraction)")


if __name__ == "__main__":
    main()
