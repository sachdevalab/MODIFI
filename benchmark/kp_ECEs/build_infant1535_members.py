#!/usr/bin/env python3
"""
Step 3 finalizer - cross-genus member table for the infant_15_35_C plasmid cluster.
Each of the 24 members with: sample, length, network-linked host species/genus (if any), and
ANI (%identity) + coverage to the representative. Shows the cluster spans multiple Enterobacteriaceae
genera - direct evidence on cross-Enterobacteriaceae transfer.
"""
import os
import re
import pandas as pd

BASE = "/home/shuaiw/borg/revision/kp_eces"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
CLUSTER_TSV = "/home/shuaiw/borg/paper/MGE/cluster/megablast.cluster.95ani.tsv"
FAI = "/home/shuaiw/borg/revision/kp_eces/kp_cluster_member_lengths.tsv"
LINKAGE = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
ANI_ROWS = os.path.join(BASE, "deepdive/infant_15_35_C.ani_rows.tsv")
REP = "infant_15_35_C"


def clean_species(s):
    s = re.sub(r"^s__", "", str(s)).strip()
    return re.sub(r"_[A-Z]+ ", " ", re.sub(r"_[A-Z]+$", "", s))


def main():
    # members of the representative cluster
    members = []
    with open(CLUSTER_TSV) as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            if rep == REP:
                members = mem.split(",")
                break
    # lengths
    fai = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["seq", "len"])
    len_map = dict(zip(fai["seq"], fai["len"]))
    # network host species per ECE (if linked)
    gc = pd.read_csv(LINKAGE)
    host_map = {}
    for _, r in gc.iterrows():
        if r["MGE"] in members:
            host_map.setdefault(r["MGE"], set()).add(clean_species(r["host_taxa"]))
    # ANI to representative
    ani = {}
    if os.path.exists(ANI_ROWS) and os.path.getsize(ANI_ROWS) > 0:
        adf = pd.read_csv(ANI_ROWS, sep="\t", header=None,
                          names=["q", "t", "num_alns", "pid", "qcov", "tcov"])
        for _, r in adf.iterrows():
            other = r["t"] if r["q"] == REP else (r["q"] if r["t"] == REP else None)
            if other is None:
                continue
            # keep the strongest (max pid*cov) record per member
            cov = max(r["qcov"], r["tcov"])
            prev = ani.get(other)
            if prev is None or r["pid"] > prev[0]:
                ani[other] = (r["pid"], cov)

    rows = []
    for m in members:
        sp = ";".join(sorted(host_map.get(m, []))) if m in host_map else ""
        genus = ";".join(sorted({s.split(" ")[0] for s in host_map.get(m, [])})) if m in host_map else ""
        pid, cov = ani.get(m, ("", ""))
        rows.append({
            "member": m,
            "is_representative": (m == REP),
            "length_bp": len_map.get(m, ""),
            "network_linked_host_species": sp,
            "network_linked_host_genus": genus,
            "pid_to_representative": pid,
            "aln_cov_to_representative": cov,
        })
    df = pd.DataFrame(rows)
    df["length_bp"] = pd.to_numeric(df["length_bp"], errors="coerce")
    df = df.sort_values("length_bp", ascending=False)
    out = os.path.join(BASE, "infant_15_35_C_members.tsv")
    df.to_csv(out, sep="\t", index=False)
    df.to_csv(os.path.join(FIG, "infant_15_35_C_members.tsv"), sep="\t", index=False)
    print(f"[members] wrote {out} ({len(df)} members)")
    print(df.to_string(index=False))
    genera = sorted({g for gs in df["network_linked_host_genus"] for g in str(gs).split(";") if g})
    print(f"\n[members] distinct network-linked host genera in cluster: {genera}")


if __name__ == "__main__":
    main()
