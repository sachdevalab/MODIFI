#!/usr/bin/env python3
"""Prepare tidy CSVs for the result figures (cluster x host-genus matrix; gene-cargo categories)."""
import os
import re
import pandas as pd

BASE = "/home/shuaiw/borg/revision/kp_eces"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
ENTERO_ORDER = ["Klebsiella", "Citrobacter", "Enterobacter", "Escherichia",
                "Serratia", "Pluralibacter", "Leclercia", "Salmonella"]


def main():
    # ---- cluster x host-genus presence (long) ----
    clu = pd.read_csv(os.path.join(BASE, "kp_22_clusters_annotated.tsv"), sep="\t")
    rows = []
    for _, r in clu.iterrows():
        genera = [g for g in str(r["linked_host_genera"]).split(";") if g]
        for g in genera:
            gg = re.sub(r"_[A-Z]+$", "", g)
            rows.append({"cluster_rep": r["cluster_rep"],
                         "n_members": r["n_cluster_members"],
                         "genus": gg})
    pd.DataFrame(rows).to_csv(os.path.join(FIG, "fig_cluster_genus_long.csv"), index=False)

    # ---- gene-cargo functional categories for the two deep-dive plasmids ----
    # case-SENSITIVE regexes (a case-insensitive [A-Z] wrongly matches trans*, etc.). AMR and
    # metal/biocide counts come from AMRFinderPlus (precise), not fuzzy product strings.
    amr = pd.read_csv(os.path.join(BASE, "annotate/amrfinder.tsv"), sep="\t")
    ctg_col = "Contig id" if "Contig id" in amr.columns else amr.columns[2]

    cats = [
        ("Conjugation / transfer (T4SS)",
         r"conjug|[Tt]ype IV secretion|[Tt]ype IV conjugativ|conjugal transfer|\bTra[A-Z]\b|\bTrb[A-Z]|relaxase|pilus assembly"),
        ("Replication", r"[Rr]eplication init|RepA|\bRep\b|origin of replication"),
        ("Partition / stability", r"[Pp]artition|ParA|ParB|StbA|plasmid segregation"),
        ("Toxin-antitoxin", r"toxin|antitoxin|CcdB|RelE|ParE|AbiE"),
        ("IS / transposon / recombinase",
         r"transposase|[Ii]nsertion sequence|\bIS\d|\bTn\d|resolvase|recombinase|integrase|Mobile element"),
        ("DNA modification (MTase / RM)",
         r"methyltransferase|methylase|restriction endonuclease|modification methylase"),
    ]
    for seq in ["infant_15_35_C", "infant_8_26_C"]:
        f = os.path.join(BASE, f"deepdive/{seq}_bakta/{seq}.tsv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, sep="\t", comment="#", header=None,
                         names=["seqid", "type", "start", "stop", "strand",
                                "locus", "gene", "product", "dbxref"])
        prods = df["product"].fillna("").astype(str)
        out = []
        assigned = pd.Series(False, index=df.index)
        for cat, pat in cats:
            hit = prods.str.contains(pat, case=True, regex=True)
            out.append({"plasmid": seq, "category": cat, "n_genes": int(hit.sum())})
            assigned = assigned | hit
        # precise resistance counts from AMRFinderPlus for this contig
        a = amr[amr[ctg_col] == seq]
        n_amr = int((a["Element type"] == "AMR").sum()) if "Element type" in a.columns else 0
        n_str = int((a["Element type"] == "STRESS").sum()) if "Element type" in a.columns else 0
        out.append({"plasmid": seq, "category": "Antibiotic resistance", "n_genes": n_amr})
        out.append({"plasmid": seq, "category": "Metal / biocide resistance", "n_genes": n_str})
        out.append({"plasmid": seq, "category": "Other / hypothetical",
                    "n_genes": int((~assigned).sum())})
        pd.DataFrame(out).to_csv(os.path.join(FIG, f"{seq}_gene_categories.csv"), index=False)
        print(f"[prep] {seq}: {len(df)} CDS; AMR={n_amr} STRESS={n_str} -> {seq}_gene_categories.csv")

    print("[prep] wrote fig_cluster_genus_long.csv")


if __name__ == "__main__":
    main()
