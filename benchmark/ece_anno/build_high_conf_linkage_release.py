#!/usr/bin/env python3
"""High-confidence ECE-host linkage release (for Zenodo / reviewer upload).
From the 317 final strict high-confidence linkages, write:
  high_conf_linkage_table.csv   (per-linkage table with cluster columns; -> data dir + figure dir)
  linkage_ece.fasta             (the 317 linked ECE contigs)
  linkage_host.fasta            (the linked host contigs only, deduplicated)
  linkage_ece_cluster.tsv       (ECE contig -> 95%-ANI MGE cluster id)   [2 cols]
  linkage_host_cluster.tsv      (host contig -> 99% dRep host cluster id) [2 cols]
  README.md                     (dataset description)
"""
import csv, os, sys
import pysam
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/network")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
import plot_linkage_data as pld
from sample_object import classify_taxa, get_detail_taxa_name

FP = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"
N = "/home/shuaiw/borg/revision/network"
RUN2 = "/home/shuaiw/borg/paper/run2"
ECE_FA = f"{N}/all_mge_revised.fa"
DATA = "/home/shuaiw/borg/revision/ece_anno/high_conf_linkage"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
os.makedirs(DATA, exist_ok=True)

lk = list(csv.DictReader(open(FP)))
ctg_taxa = pld.get_ctg_taxa(f"{RUN2}/")
mge_clu, _ = pld.read_mge_cluster(f"{N}/MGE_cluster/megablast.cluster.95ani.tsv")
_, host_clu, _ = pld.read_drep_cluster(f"{N}/rescued_network/dRep_99_out/data_tables/Cdb.csv")

# --- fasta handles ---
ece_fh = pysam.FastaFile(ECE_FA)              # builds .fai if absent
host_fh = {}
def host_handle(sample):
    if sample not in host_fh:
        p = f"{RUN2}/{sample}/{sample}.hifiasm.p_ctg.rename.fa"
        host_fh[sample] = pysam.FastaFile(p) if os.path.exists(p) else None
    return host_fh[sample]

def sp(host):
    s = classify_taxa(ctg_taxa.get(host, "NA"), "species")
    return s.replace("s__", "") if s not in ("Unknown", "s__", "", None) and len(s) > 3 else get_detail_taxa_name(ctg_taxa.get(host, "NA"))

# --- table ---
rows = []
for r in lk:
    mge, host, s = r["MGE"], r["host"], r["sample"]
    hh = host_handle(s)
    mlen = ece_fh.get_reference_length(mge) if mge in ece_fh.references else ""
    hlen = hh.get_reference_length(host) if (hh and host in hh.references) else ""
    rows.append({
        "sample": s, "environment": r["environment"],
        "MGE": mge, "MGE_cluster": mge_clu.get(mge, mge), "MGE_type": r["type"],
        "MGE_len": mlen, "MGE_gc": r["MGE_gc"], "MGE_cov": r["MGE_cov"],
        "host": host, "host_cluster": host_clu.get(host, host),
        "host_species": sp(host), "host_phylum": r["host_phylum"],
        "host_taxa": get_detail_taxa_name(ctg_taxa.get(host, "NA")), "host_len": hlen,
        "host_gc": r["host_gc"], "host_cov": r["host_cov"],
        "cos_sim": r["cos_sim"]})
cols = list(rows[0].keys())
for outp in (f"{DATA}/high_conf_linkage_table.csv", f"{FIG}/high_conf_linkage_table.csv"):
    with open(outp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

# --- sequences ---
def write_fa(handle_of, names_samples, path):
    n = 0
    with open(path, "w") as out:
        for name, hnd in names_samples:
            if hnd and name in hnd.references:
                out.write(f">{name}\n{hnd.fetch(name)}\n"); n += 1
            else:
                print("  MISSING seq:", name)
    return n

eces = sorted(set(r["MGE"] for r in lk))
n_ece = write_fa(None, [(m, ece_fh) for m in eces], f"{DATA}/linkage_ece.fasta")
hosts = sorted(set((r["host"], r["sample"]) for r in lk))
n_host = write_fa(None, [(h, host_handle(s)) for h, s in hosts], f"{DATA}/linkage_host.fasta")

# --- cluster maps (2 cols each) ---
with open(f"{DATA}/linkage_ece_cluster.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["contig", "cluster_id"])
    for m in eces: w.writerow([m, mge_clu.get(m, m)])
with open(f"{DATA}/linkage_host_cluster.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["contig", "cluster_id"])
    for h, _s in hosts: w.writerow([h, host_clu.get(h, h)])

print(f"linkages: {len(rows)}")
print(f"ECE contigs: {n_ece}/{len(eces)}   host contigs: {n_host}/{len(hosts)}")
print(f"MGE clusters: {len(set(mge_clu.get(m, m) for m in eces))}   "
      f"host clusters: {len(set(host_clu.get(h, h) for h, _ in hosts))}")
print(f"table -> {DATA}/high_conf_linkage_table.csv  (+ figure dir)")
print(f"seqs  -> linkage_ece.fasta, linkage_host.fasta; maps -> linkage_{{ece,host}}_cluster.tsv")
