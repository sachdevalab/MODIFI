#!/usr/bin/env python3
"""
Step 2 finalizer - merge mob_typer + AMRFinderPlus + existing CONJscan/geNomad evidence onto the
Klebsiella-associated ECE clusters, producing the master answer to reviewer parts 1-2.

Outputs:
  kp_22_clusters_annotated.tsv  - one row per cluster: rep, size, host reach, mobility, replicon /
                                  relaxase MOB type, closest known plasmid + mash distance +
                                  predicted host range, AMR/virulence gene count and list.
  kp_cluster_members_annotated.tsv - per-member ECE annotation (mob_typer + AMR).
"""
import os
import re
import pandas as pd

BASE = "/home/shuaiw/borg/revision/kp_eces"
ANN = os.path.join(BASE, "annotate")
FIG_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
EVID = "/home/shuaiw/borg/revision/ece_anno/metagenome/ece_evidence_all.tsv"


def safe_read_tsv(path, **kw):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return pd.read_csv(path, sep="\t", **kw)
    return pd.DataFrame()


def main():
    clu = pd.read_csv(os.path.join(BASE, "kp_22_clusters.tsv"), sep="\t")
    mob = safe_read_tsv(os.path.join(ANN, "mob_typer.tsv"))
    amr = safe_read_tsv(os.path.join(ANN, "amrfinder.tsv"))
    evid = safe_read_tsv(EVID)

    # ---- per-member mob_typer ----
    # mob_typer 'sample_id' column holds the record id (the ECE contig name)
    mob_id_col = next((c for c in ["sample_id", "id", "sample"] if c in mob.columns), None)
    mob_keep = {}
    if mob_id_col:
        mcols = [c for c in ["num_contigs", "size", "rep_type(s)", "relaxase_type(s)",
                             "mpf_type", "orit_type(s)", "predicted_mobility",
                             "mash_nearest_neighbor", "mash_neighbor_distance",
                             "mash_neighbor_identification",
                             "primary_cluster_id", "secondary_cluster_id",
                             "predicted_host_range_overall_rank",
                             "predicted_host_range_overall_name"] if c in mob.columns]
        mob_keep = mob.set_index(mob_id_col)[mcols].to_dict("index")

    # ---- per-member AMR / stress (metal+biocide) ----
    # AMRFinderPlus 'Element type' separates true antibiotic AMR from STRESS (metal/biocide) and
    # VIRULENCE; the reviewer asks specifically about AMR, so track them separately.
    amr_by_ece = {}      # ece -> {AMR:[genes], STRESS:[genes], VIRULENCE:[genes]}
    if len(amr):
        # prefer the contig column (nucleotide mode fills 'Contig id'; 'Protein identifier' is NA)
        ctg_col = next((c for c in amr.columns if c.lower() in ("contig id", "contig")), None) \
            or next((c for c in amr.columns if c.lower() in ("protein identifier", "protein id")), None)
        gene_col = next((c for c in amr.columns if c.lower() in
                         ("gene symbol", "element symbol", "gene")), None)
        type_col = next((c for c in amr.columns if c.lower() == "element type"), None)
        if ctg_col and gene_col:
            for _, r in amr.iterrows():
                ece = str(r[ctg_col])
                etype = str(r[type_col]) if type_col else "AMR"
                amr_by_ece.setdefault(ece, {}).setdefault(etype, []).append(str(r[gene_col]))

    # existing evidence (mobility, replicon, partition) keyed by seq_name
    evid_map = {}
    if len(evid):
        for _, r in evid.iterrows():
            evid_map[r["seq_name"]] = r

    # ---- per-member annotated table ----
    members = pd.read_csv(os.path.join(BASE, "kp_all_cluster_member_eces.txt"),
                          header=None, names=["ece"])
    mrows = []
    for ece in members["ece"]:
        m = mob_keep.get(ece, {})
        ev = evid_map.get(ece, {})
        by_type = amr_by_ece.get(ece, {})
        amr_genes = sorted(set(by_type.get("AMR", [])))
        stress_genes = sorted(set(by_type.get("STRESS", [])))
        vir_genes = sorted(set(by_type.get("VIRULENCE", [])))
        mrows.append({
            "ece": ece,
            "mob_size": m.get("size", ""),
            "rep_type": m.get("rep_type(s)", ""),
            "relaxase_type": m.get("relaxase_type(s)", ""),
            "mpf_type": m.get("mpf_type", ""),
            "predicted_mobility": m.get("predicted_mobility", ""),
            "mash_nearest_neighbor": m.get("mash_nearest_neighbor", ""),
            "mash_neighbor_distance": m.get("mash_neighbor_distance", ""),
            "mash_neighbor_species": m.get("mash_neighbor_identification", ""),
            "predicted_host_range_rank": m.get("predicted_host_range_overall_rank", ""),
            "predicted_host_range_name": m.get("predicted_host_range_overall_name", ""),
            "conjscan_type": ev.get("conjscan_type", "") if len(ev) else "",
            "plas_rep_initiation": ev.get("plas_rep_initiation", "") if len(ev) else "",
            "plas_partition": ev.get("plas_partition", "") if len(ev) else "",
            "n_amr_genes": len(amr_genes),
            "amr_genes": ";".join(amr_genes),
            "n_stress_genes": len(stress_genes),
            "stress_genes": ";".join(stress_genes),
            "n_virulence_genes": len(vir_genes),
            "virulence_genes": ";".join(vir_genes),
        })
    mem_df = pd.DataFrame(mrows)
    mem_df.to_csv(os.path.join(BASE, "kp_cluster_members_annotated.tsv"), sep="\t", index=False)

    # ---- cluster-level rollup ----
    mem_df["_ece"] = mem_df["ece"]
    reps = clu["cluster_rep"].tolist()
    crows = []
    # map each member to its cluster via the cluster membership file
    memb_to_clu = {}
    for _, r in clu.iterrows():
        for mm in str(r["all_cluster_members"]).split(";"):
            memb_to_clu[mm] = r["cluster_rep"]
    mem_df["cluster_rep"] = mem_df["ece"].map(memb_to_clu)

    for _, r in clu.iterrows():
        rep = r["cluster_rep"]
        sub = mem_df[mem_df["cluster_rep"] == rep]
        rep_row = mem_df[mem_df["ece"] == rep]
        def uniq(col, sep=";"):
            return sorted(set(g for gs in sub[col] for g in str(gs).split(sep) if g))
        all_amr = uniq("amr_genes")
        all_stress = uniq("stress_genes")
        all_vir = uniq("virulence_genes")
        mobil = sub["predicted_mobility"].replace("", pd.NA).dropna()
        reptypes = sorted(set(x for xs in sub["rep_type"] for x in str(xs).split(",") if x and x != "-"))
        relaxases = sorted(set(x for xs in sub["relaxase_type"] for x in str(xs).split(",") if x and x != "-"))
        crows.append({
            "cluster_rep": rep,
            "rep_type": r["rep_type"],
            "rep_len_bp": r["rep_len_bp"],
            "n_cluster_members": r["n_cluster_members"],
            "n_distinct_host_genera": r["n_distinct_host_genera"],
            "linked_host_genera": r["linked_host_genera"],
            "linked_host_species": r["linked_host_species"],
            "predicted_mobility": ";".join(sorted(set(mobil))) if len(mobil) else "",
            "replicon_types": ";".join(reptypes),
            "relaxase_MOB_types": ";".join(relaxases),
            "rep_mash_nearest_neighbor": rep_row["mash_nearest_neighbor"].iloc[0] if len(rep_row) else "",
            "rep_mash_distance": rep_row["mash_neighbor_distance"].iloc[0] if len(rep_row) else "",
            "rep_mash_neighbor_species": rep_row["mash_neighbor_species"].iloc[0] if len(rep_row) else "",
            "rep_predicted_host_range": rep_row["predicted_host_range_name"].iloc[0] if len(rep_row) else "",
            "n_amr_genes_cluster": len(all_amr),
            "amr_genes_cluster": ";".join(all_amr),
            "n_stress_genes_cluster": len(all_stress),
            "stress_genes_cluster": ";".join(all_stress),
            "n_virulence_genes_cluster": len(all_vir),
            "virulence_genes_cluster": ";".join(all_vir),
        })
    cdf = pd.DataFrame(crows).sort_values("n_cluster_members", ascending=False)
    out = os.path.join(BASE, "kp_22_clusters_annotated.tsv")
    cdf.to_csv(out, sep="\t", index=False)
    cdf.to_csv(os.path.join(FIG_DIR, "kp_22_clusters_annotated.tsv"), sep="\t", index=False)
    print(f"[merge] wrote {out}")
    print(cdf[["cluster_rep", "n_cluster_members", "predicted_mobility",
               "n_amr_genes_cluster", "n_stress_genes_cluster", "rep_mash_distance",
               "rep_mash_neighbor_species"]].to_string(index=False))
    print(f"\n[merge] clusters with >=1 antibiotic-AMR gene: {(cdf['n_amr_genes_cluster']>0).sum()}/{len(cdf)}")
    print(f"[merge] antibiotic-AMR genes across clusters: "
          f"{sorted(set(g for gs in cdf['amr_genes_cluster'] for g in str(gs).split(';') if g))}")
    print(f"[merge] clusters with >=1 metal/biocide-STRESS gene: {(cdf['n_stress_genes_cluster']>0).sum()}/{len(cdf)}")
    print(f"[merge] clusters that are known plasmids (rep mash dist < 0.05): "
          f"{(pd.to_numeric(cdf['rep_mash_distance'], errors='coerce') < 0.05).sum()}/{len(cdf)}")


if __name__ == "__main__":
    main()
