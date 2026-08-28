#!/usr/bin/env python3
"""IMG-catalogue host cross-check of the FINAL strict linkage set (317 linked ECEs), adapted from
benchmark/kp_ECEs/recompute_strict.py. KNOWN is defined at the SAME stringency as the ECE clusters
(cluster_MGE.py): a strict hit = anicalc weighted-ANI (pid) >= 95 AND query coverage (qcov) >= 85.
Funnel per MGE type: total linked -> known(strict) -> host genus-resolved (strict hit whose reference
carries a genus host) -> host-supported (inferred genus matches ANY such reference genus).
Writes *_strict_validation.tsv + ece_validation_summary_strict.csv under img_map/, and the summary
into the figure dir."""
import csv, re
import pandas as pd

D = "/home/shuaiw/borg/revision/ece_anno/expanded/img_map"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
ANI, QCOV = 95.0, 85.0   # cluster_MGE.py stringency (ANI>=95, query cov>=85)

def genus(lin):
    for t in str(lin).split(";"):
        t = t.strip()
        if t.startswith("g__"):
            return re.sub(r"_[A-Z]+$", "", t[3:])
    return ""

def uvig(s):
    return str(s).split("|")[0].split()[0]

def run(ani_tsv, hostfile, host_key, host_col, mge_type, ntot, label):
    meta = pd.read_csv(f"{D}/linked_eces_meta.tsv", sep="\t")
    sub = meta[meta.MGE_type == mge_type].copy()
    inf = dict(zip(sub.MGE, sub.host_genus.fillna("")))
    b = pd.read_csv(ani_tsv, sep="\t")            # cols: qname,tname,num_alns,pid,qcov,tcov
    b["refid"] = b["tname"].map(uvig)
    # stream host table for the referenced ids only
    want = set(b["refid"]); hmap = {}
    with open(hostfile) as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        ik = hdr.index(host_key); ih = hdr.index(host_col)
        for row in rd:
            if len(row) > ih and row[ik] in want:
                hmap[row[ik]] = row[ih]
    b["ref_g"] = b["refid"].map(hmap).map(genus)
    b["strict"] = (b.pid >= ANI) & (b.qcov >= QCOV)
    rows = []
    for q, g in b.groupby("qname"):
        if not bool(g.strict.any()):
            continue                              # -> counted as "no strict match"
        gs = g[g.strict & (g.ref_g != "")]        # strict hits carrying a genus host
        if len(gs):
            ig = inf.get(q, ""); genera = set(gs.ref_g)
            supp = bool(ig) and ig in genera
            rg = ig if supp else gs.iloc[0].ref_g
            cat = "host-supported (agrees)" if supp else "host mismatch"
            rows.append((q, True, rg, ig, cat))
        else:
            rows.append((q, True, "", inf.get(q, ""), "known (strict), host not genus-resolved"))
    per = pd.DataFrame(rows, columns=["MGE", "known_strict", "ref_host_genus", "inferred_host_genus", "category"])
    per.to_csv(f"{D}/{label}_strict_validation.tsv", sep="\t", index=False)
    n_known = int(per.known_strict.sum())
    res = per[per.category.isin(["host-supported (agrees)", "host mismatch"])]
    supp = int((res.category == "host-supported (agrees)").sum())
    n_nores = int((per.category == "known (strict), host not genus-resolved").sum())
    counts = {"host-supported (agrees)": supp, "host mismatch": len(res) - supp,
              "known (strict), host not genus-resolved": n_nores, "no strict match": ntot - n_known}
    rr = f"{100*supp/len(res):.0f}%" if len(res) else "NA"
    print(f"{label}: total={ntot} known(strict)={n_known} ({100*n_known/ntot:.0f}%) "
          f"host-resolved={len(res)} supported={supp} ({rr} of resolved)")
    return counts

def main():
    meta = pd.read_csv(f"{D}/linked_eces_meta.tsv", sep="\t")
    n_pl = int((meta.MGE_type == "plasmid").sum())
    n_vi = int((meta.MGE_type == "virus").sum())
    vc = run(f"{D}/virus_imgvr_ani.tsv",
             "/shared/db/imgvr/v4_IMG_VR_2022-09-20_6/metadata/IMGVR_all_Host_information.tsv",
             "UVIG", "Host taxonomy prediction", "virus", n_vi, "virus_imgvr")
    pc = run(f"{D}/plasmid_imgpr_ani.tsv",
             "/shared/db/imgpr/2023-08-08_1/IMGPR_plasmid_data.tsv",
             "plasmid_id", "host_taxonomy", "plasmid", n_pl, "plasmid_imgpr")
    order = ["host-supported (agrees)", "host mismatch",
             "known (strict), host not genus-resolved", "no strict match"]
    rows = [("Viruses (IMG/VR)", c, vc[c]) for c in order] + \
           [("Plasmids (IMG/PR)", c, pc[c]) for c in order]
    out = pd.DataFrame(rows, columns=["type", "category", "n"])
    out.to_csv(f"{D}/ece_validation_summary_strict.csv", index=False)
    out.to_csv(f"{FIG}/ece_validation_summary_strict.csv", index=False)
    print(f"wrote ece_validation_summary_strict.csv (plasmid n={n_pl}, virus n={n_vi})")

if __name__ == "__main__":
    main()
