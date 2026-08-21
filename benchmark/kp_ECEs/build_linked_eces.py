#!/usr/bin/env python3
"""Generalise the Klebsiella-ECE analysis to ALL confidently modification-linked ECEs in the network.
Builds the metadata table (env, sample, type, length, inferred host taxonomy, 95%-ANI cluster) and the
ID list used to pull sequences for mob_typer / AMRFinderPlus."""
import os
import re
from collections import defaultdict
import pandas as pd

LINKAGE = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
CLUSTER_TSV = "/home/shuaiw/borg/paper/MGE/cluster/megablast.cluster.95ani.tsv"
OUT = "/home/shuaiw/borg/revision/linked_eces"


def read_mge_cluster(f):
    m = {}
    df = pd.read_csv(f, sep="\t", header=None)
    for _, r in df.iterrows():
        for mge in str(r[1]).split(","):
            m[mge] = r[0]
    return m


def sp(t):
    return re.sub(r"^s__", "", str(t)).strip()


def genus(t):
    s = sp(t)
    return re.sub(r"_[A-Z]+$", "", s.split(" ")[0]) if s else "NA"


def main():
    gc = pd.read_csv(LINKAGE)
    clu = read_mge_cluster(CLUSTER_TSV)
    gc["cluster"] = gc["MGE"].map(clu)
    gc["host_species"] = gc["host_taxa"].map(sp)
    gc["host_genus"] = gc["host_taxa"].map(genus)
    keep = ["MGE", "MGE_type", "environment", "sample", "mge_len",
            "host", "host_taxa", "host_species", "host_genus", "host_phylum", "cluster"]
    meta = gc[keep].drop_duplicates("MGE").sort_values(["MGE_type", "environment", "mge_len"],
                                                       ascending=[True, True, False])
    meta.to_csv(os.path.join(OUT, "linked_eces_meta.tsv"), sep="\t", index=False)
    with open(os.path.join(OUT, "linked_eces.txt"), "w") as fh:
        fh.write("\n".join(meta["MGE"]) + "\n")
    print(f"[linked] {len(meta)} linked ECEs "
          f"({(meta.MGE_type=='plasmid').sum()} plasmid, {(meta.MGE_type=='virus').sum()} virus, "
          f"{(meta.MGE_type=='novel').sum()} novel)")
    print(f"[linked] environments: {meta['environment'].value_counts().to_dict()}")
    print(f"[linked] distinct host species: {meta['host_species'].nunique()}, "
          f"genera: {meta['host_genus'].nunique()}")


if __name__ == "__main__":
    main()
