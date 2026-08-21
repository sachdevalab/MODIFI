#!/usr/bin/env python3
"""Merge mob_typer + AMRFinderPlus onto all linked ECEs, and validate modification-inferred hosts
against the documented host of the closest known plasmid (the generalisation of the Klebsiella 8/9
result to the whole network)."""
import os
import re
import pandas as pd

BASE = "/home/shuaiw/borg/revision/linked_eces"
ANN = os.path.join(BASE, "annotate")
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
KNOWN_MAX = 0.05   # Mash distance threshold for "known plasmid"


def genus_of(s):
    s = str(s)
    if s in ("", "nan", "NA"):
        return ""
    return re.sub(r"_[A-Z]+$", "", s.split(" ")[0])


def main():
    meta = pd.read_csv(os.path.join(BASE, "linked_eces_meta.tsv"), sep="\t")
    mob = pd.read_csv(os.path.join(ANN, "mob_typer.tsv"), sep="\t")
    amr = pd.read_csv(os.path.join(ANN, "amrfinder.tsv"), sep="\t")

    mob = mob.rename(columns={"sample_id": "MGE"})
    mcols = ["MGE", "size", "rep_type(s)", "relaxase_type(s)", "mpf_type", "predicted_mobility",
             "mash_nearest_neighbor", "mash_neighbor_distance", "mash_neighbor_identification",
             "predicted_host_range_overall_name"]
    m = meta.merge(mob[[c for c in mcols if c in mob.columns]], on="MGE", how="left")

    # AMR / stress counts per ECE
    ctg = "Contig id" if "Contig id" in amr.columns else amr.columns[2]
    et = "Element type"
    amr_c = amr[amr[et] == "AMR"].groupby(ctg)["Gene symbol"].apply(lambda s: ";".join(sorted(set(s))))
    str_c = amr[amr[et] == "STRESS"].groupby(ctg)["Gene symbol"].apply(lambda s: ";".join(sorted(set(s))))
    m["amr_genes"] = m["MGE"].map(amr_c).fillna("")
    m["stress_genes"] = m["MGE"].map(str_c).fillna("")
    m["n_amr"] = m["amr_genes"].apply(lambda s: 0 if s == "" else len(s.split(";")))
    m["n_stress"] = m["stress_genes"].apply(lambda s: 0 if s == "" else len(s.split(";")))

    # known-plasmid + inferred-vs-reference host match (plasmids only; mash DB is plasmids)
    m["mash_dist"] = pd.to_numeric(m["mash_neighbor_distance"], errors="coerce")
    m["is_plasmid"] = m["MGE_type"] == "plasmid"
    m["known_plasmid"] = m["is_plasmid"] & (m["mash_dist"] < KNOWN_MAX)
    m["ref_host_genus"] = m["mash_neighbor_identification"].map(genus_of)
    m["inferred_host_genus"] = m["host_genus"]
    def match(r):
        if not r["known_plasmid"]:
            return ""
        return "yes" if r["ref_host_genus"] and r["ref_host_genus"] == r["inferred_host_genus"] else "no"
    m["ref_matches_inferred_genus"] = m.apply(match, axis=1)

    out = os.path.join(BASE, "all_linked_eces_annotated.tsv")
    m.to_csv(out, sep="\t", index=False)
    m.to_csv(os.path.join(FIG, "all_linked_eces_annotated.tsv"), sep="\t", index=False)

    # ---- summaries ----
    pl = m[m.is_plasmid]
    print(f"[merge-linked] {len(m)} linked ECEs ({len(pl)} plasmids, {(~m.is_plasmid).sum()} viruses)")
    print(f"[merge-linked] plasmids matching a known plasmid (mash<{KNOWN_MAX}): "
          f"{m['known_plasmid'].sum()}/{len(pl)} ({100*m['known_plasmid'].sum()/len(pl):.0f}%)")
    print(f"[merge-linked] mobility (plasmids): {pl['predicted_mobility'].fillna('NA').apply(lambda s: 'conjugative' if 'conjugative' in str(s) else ('mobilizable' if re.search(r'(^|;)mobilizable', str(s)) else 'non-mobilizable')).value_counts().to_dict()}")
    print(f"[merge-linked] ECEs with >=1 antibiotic-AMR gene: {(m.n_amr>0).sum()}/{len(m)}; "
          f"with metal/biocide: {(m.n_stress>0).sum()}/{len(m)}")
    kn = m[m.known_plasmid]
    vc = kn["ref_matches_inferred_genus"].value_counts().to_dict()
    tot = vc.get("yes", 0) + vc.get("no", 0)
    print(f"\n[merge-linked] KNOWN-vs-INFERRED HOST VALIDATION (known plasmids, genus level): "
          f"{vc.get('yes',0)}/{tot} agree ({100*vc.get('yes',0)/tot:.0f}%)")
    # agreement by environment
    print("[merge-linked] agreement by environment:")
    for env, g in kn.groupby("environment"):
        y = (g.ref_matches_inferred_genus == "yes").sum(); n = (g.ref_matches_inferred_genus == "no").sum()
        if y + n:
            print(f"    {env}: {y}/{y+n}")
    # AMR gene repertoire
    allamr = sorted(set(x for gs in m.amr_genes for x in str(gs).split(";") if x))
    print(f"\n[merge-linked] AMR gene repertoire across all linked ECEs ({len(allamr)}): {allamr}")


if __name__ == "__main__":
    main()
