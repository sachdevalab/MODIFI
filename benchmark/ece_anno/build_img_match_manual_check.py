#!/usr/bin/env python3
"""Per-linkage manual-inspection table for the IMG host cross-check (panel j of ece_profile_final).
For every 'comparable' linkage (an ECE with a strict IMG hit -- anicalc pid>=95 & qcov>=85 -- so a host
comparison is possible), list our ECE + host contig + our host genus/full GTDB taxonomy, the matched IMG
reference (db ECE), the reference's host genus/full taxonomy, the alignment ANI (pid) and query fraction
(qcov), and whether our host genus agrees with the reference (is_match). Output CSV -> figure dir."""
import csv, re, pickle
import pandas as pd

D = "/home/shuaiw/borg/revision/ece_anno/expanded/img_map"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
LK = "/home/shuaiw/borg/revision/ece_anno/high_conf_linkage/high_conf_linkage_table.csv"
TAXA_PKL = "/home/shuaiw/borg/paper/gene_anno/meta_ctg_taxa_dict.pkl"
ANI, QCOV = 95.0, 85.0

IMGVR_HOST = "/shared/db/imgvr/v4_IMG_VR_2022-09-20_6/metadata/IMGVR_all_Host_information.tsv"
IMGPR_HOST = "/shared/db/imgpr/2023-08-08_1/IMGPR_plasmid_data.tsv"


def genus_full(lin):
    """genus token as GTDB writes it, keeping any _A/_B/_E subgenus suffix."""
    for t in str(lin).split(";"):
        t = t.strip()
        if t.startswith("g__"):
            return t[3:]
    return ""


def base_genus(g):
    """strip GTDB alphabet subgenus suffix so Ruminococcus_B == Ruminococcus."""
    return re.sub(r"_[A-Z]+$", "", str(g))


def genus(lin):
    return base_genus(genus_full(lin))


def uvig(s):
    return str(s).split("|")[0].split()[0]


def load_host_table(path, key_col, tax_col, want):
    """refid -> full host taxonomy string, for referenced ids only."""
    hmap = {}
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        ik, it = hdr.index(key_col), hdr.index(tax_col)
        for row in rd:
            if len(row) > it and row[ik] in want:
                hmap[row[ik]] = row[it].strip()
    return hmap


def process(ani_tsv, host_path, key_col, tax_col, mge_type, db_label, lk, taxa):
    b = pd.read_csv(ani_tsv, sep="\t")                 # qname,tname,num_alns,pid,qcov,tcov
    b["refid"] = b["tname"].map(uvig)
    hmap = load_host_table(host_path, key_col, tax_col, set(b["refid"]))
    b["ref_tax"] = b["refid"].map(hmap).fillna("")
    b["ref_g"] = b["ref_tax"].map(genus)
    b["strict"] = (b.pid >= ANI) & (b.qcov >= QCOV)

    rows = []
    for ece, g in b.groupby("qname"):
        gs = g[g.strict]
        if not len(gs):
            continue                                   # not comparable (no strict IMG hit)
        info = lk.get(ece)
        if info is None:
            continue
        host_ctg = info["host"]
        our_tax = taxa.get(host_ctg, "Unknown")
        our_g_full = genus_full(our_tax)                # GTDB token, keeps _A/_B suffix
        our_g = base_genus(our_g_full)                  # base genus for the match test
        gsr = gs[gs.ref_g != ""]                         # strict hits with a genus-resolved ref host
        resolved = len(gsr) > 0
        if resolved:
            match_hits = gsr[gsr.ref_g == our_g]
            pick_pool = match_hits if len(match_hits) else gsr
            best = pick_pool.sort_values("pid", ascending=False).iloc[0]
            is_match = bool(len(match_hits) > 0) and bool(our_g)
        else:
            best = gs.sort_values("pid", ascending=False).iloc[0]
            is_match = ""
        # flag agreements that hold only after stripping the GTDB subgenus suffix
        note = ""
        if is_match is True and our_g_full != best["ref_g"]:
            note = f"base-genus match (GTDB subgenus {our_g_full} vs {best['ref_g']})"
        rows.append({
            "sample": info["sample"], "environment": info["environment"],
            "ece": ece, "ece_type": mge_type, "db": db_label,
            "our_host_contig": host_ctg, "our_host_genus": our_g_full,
            "our_host_full_taxonomy": our_tax,
            "db_ece": best["refid"], "db_ref_full_id": best["tname"],
            "db_host_genus": best["ref_g"],
            "db_host_full_taxonomy": best["ref_tax"],
            "ANI_pid": round(float(best["pid"]), 2),
            "fraction_qcov": round(float(best["qcov"]), 2),
            "target_cov_tcov": round(float(best["tcov"]), 2),
            "is_match": is_match, "match_note": note,
            "db_host_genus_resolved": "yes" if resolved else "no",
            "n_strict_hits": int(len(gs)),
            "all_strict_ref_genera": ";".join(sorted({x for x in gsr.ref_g if x})),
        })
    return rows


def main():
    lkdf = pd.read_csv(LK)
    lk = {r["MGE"]: {"sample": r["sample"], "environment": r["environment"],
                     "host": r["host"], "host_genus_col": ""} for _, r in lkdf.iterrows()}
    taxa = pickle.load(open(TAXA_PKL, "rb"))

    rows = process(f"{D}/plasmid_imgpr_ani.tsv", IMGPR_HOST, "plasmid_id", "host_taxonomy",
                   "plasmid", "IMG/PR", lk, taxa)
    rows += process(f"{D}/virus_imgvr_ani.tsv", IMGVR_HOST, "UVIG", "Host taxonomy prediction",
                    "virus", "IMG/VR", lk, taxa)
    out = pd.DataFrame(rows)
    # comparable = a strict IMG hit whose reference host is genus-resolved, so a host comparison is possible
    out = out[out.db_host_genus_resolved == "yes"].copy()
    out = out.sort_values(["ece_type", "is_match", "ece"], ascending=[True, False, True])
    cols = ["sample", "environment", "ece", "ece_type", "db",
            "our_host_contig", "our_host_genus", "our_host_full_taxonomy",
            "db_ece", "db_ref_full_id", "db_host_genus", "db_host_full_taxonomy",
            "ANI_pid", "fraction_qcov", "target_cov_tcov", "is_match", "match_note",
            "n_strict_hits", "all_strict_ref_genera"]
    out = out[cols]
    path = f"{FIG}/img_match_manual_check.csv"
    out.to_csv(path, index=False)
    print(f"wrote {path}")
    print(f"comparable (genus-resolved) linkages={len(out)}  "
          f"match={int((out.is_match==True).sum())} mismatch={int((out.is_match==False).sum())}")
    print(f"  by type: plasmid={int((out.ece_type=='plasmid').sum())} "
          f"virus={int((out.ece_type=='virus').sum())}")


if __name__ == "__main__":
    main()
