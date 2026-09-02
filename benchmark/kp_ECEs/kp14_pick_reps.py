#!/usr/bin/env python3
"""Select one representative genome per K. pneumoniae-related ECE cluster (revised linkage).
14 clusters (Kp-linked ECEs -> 95% ANI clusters). Representative rule: prefer circular; else longest;
among circular pick longest. Circularity/length from the union of ece_evidence_all.tsv tables."""
import os, glob, csv, re
import pandas as pd

BASE = "/home/shuaiw/borg/revision"
NR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/mge_host_gc_cov.csv"
CMAP = f"{BASE}/ece_anno/high_conf_linkage/linkage_ece_cluster.tsv"
RCLUST = f"{BASE}/network/MGE_cluster/megablast.cluster.95ani.tsv"
EV_GLOB = [f"{BASE}/ece_anno/ece_evidence_all.tsv"] + glob.glob(f"{BASE}/ece_anno/*/ece_evidence_all.tsv")
OUT = f"{BASE}/kp_eces/gene_profile"
os.makedirs(OUT, exist_ok=True)


def genus(taxa):
    s = re.sub(r"^s__", "", str(taxa)).strip()
    return re.sub(r"_[A-Z]+$", "", s.split(" ")[0]) if s else ""


def main():
    nr = pd.read_csv(NR)
    # Kp-linked MGEs -> their clusters
    kp_mges = set(nr[nr["host_taxa"].astype(str).str.contains("Klebsiella pneumoniae", regex=False)]["MGE"])
    cmap = pd.read_csv(CMAP, sep="\t")           # contig, cluster_id
    mge2clu = dict(zip(cmap["contig"], cmap["cluster_id"]))
    kp_clusters = sorted({mge2clu[m] for m in kp_mges if m in mge2clu})
    print(f"[kp14] Kp-linked ECEs: {len(kp_mges)}; Kp clusters: {len(kp_clusters)}")

    # full membership of each cluster (revised 95% ANI cluster file: rep<TAB>comma-members)
    clu_members = {}
    with open(RCLUST) as fh:
        for line in fh:
            rep, mem = line.rstrip("\n").split("\t")
            if rep in kp_clusters:
                clu_members[rep] = mem.split(",")

    # evidence map: seq_name -> (circular_bool, length, type); union of all evidence tables
    ev = {}
    for f in EV_GLOB:
        if not os.path.exists(f):
            continue
        with open(f) as fh:
            rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
            i_s, i_c, i_l, i_t = (hdr.index("seq_name"), hdr.index("circular"),
                                  hdr.index("length"), hdr.index("type"))
            for row in rd:
                if len(row) <= max(i_s, i_c, i_l, i_t):
                    continue
                s = row[i_s]
                circ = str(row[i_c]).strip().lower() in ("true", "1", "yes")
                try:
                    ln = float(row[i_l])
                except ValueError:
                    ln = 0.0
                prev = ev.get(s)
                # keep the record with the longer length / any circular=True
                if prev is None or ln > prev[1] or (circ and not prev[0]):
                    ev[s] = (circ, ln, row[i_t])

    # network host info for members that are linked (host species + genus)
    host_sp = {}; host_gen = {}
    for _, r in nr.iterrows():
        host_sp.setdefault(r["MGE"], set()).add(re.sub(r"^s__", "", str(r["host_taxa"])))
        host_gen.setdefault(r["MGE"], set()).add(genus(r["host_taxa"]))

    rows = []
    for clu in kp_clusters:
        members = clu_members.get(clu, [clu])
        best = None  # (circ, length, seq)
        for m in members:
            circ, ln, _ = ev.get(m, (False, 0.0, ""))
            key = (1 if circ else 0, ln)
            if best is None or key > best[0]:
                best = (key, m, circ, ln)
        rep = best[1]; rep_circ = best[2]; rep_len = best[3]
        rep_type = ev.get(rep, (None, None, ""))[2]
        kp_linked = [m for m in members if m in kp_mges]
        linked_genera = sorted({g for m in members for g in host_gen.get(m, set()) if g})
        linked_species = sorted({s for m in members for s in host_sp.get(m, set()) if s and s != "nan"})
        n_circ = sum(1 for m in members if ev.get(m, (False,))[0])
        rows.append({
            "cluster": clu, "representative": rep, "rep_type": rep_type,
            "rep_circular": rep_circ, "rep_length": int(rep_len), "n_members": len(members),
            "n_circular_members": n_circ, "rep_sample": re.sub(r"_[0-9]+_[CL]$", "", rep),
            "kp_linked_members": ";".join(kp_linked),
            "linked_host_genera": ";".join(linked_genera),
            "linked_host_species": ";".join(linked_species),
        })
    df = pd.DataFrame(rows).sort_values("n_members", ascending=False)
    out = f"{OUT}/kp14_representatives.tsv"
    df.to_csv(out, sep="\t", index=False)
    # member list for downstream (also write full membership map)
    with open(f"{OUT}/kp14_reps.txt", "w") as fh:
        fh.write("\n".join(df["representative"]) + "\n")
    print(f"[kp14] wrote {out}")
    print(df[["cluster","representative","rep_type","rep_circular","rep_length","n_members","n_circular_members"]].to_string(index=False))
    print(f"\n[kp14] representatives circular: {int(df['rep_circular'].sum())}/{len(df)}")


if __name__ == "__main__":
    main()
