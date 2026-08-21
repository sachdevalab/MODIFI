#!/usr/bin/env python3
"""Recompute host-validation with KNOWN = strict best hit (>=90% id & >=50% qcov), consistent across
the funnel: total -> known(strict) -> host genus-resolved (a strict hit whose reference has a genus
host) -> host-supported (genus agrees). Updates the summary CSV and per-ECE validation TSVs."""
import csv, re
import pandas as pd

B="/home/shuaiw/borg/revision/linked_eces"; FIG="/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
PID,QCOV=90.0,50.0

def genus(lin):
    for t in str(lin).split(";"):
        t=t.strip()
        if t.startswith("g__"): return re.sub(r"_[A-Z]+$","",t[3:])
    return ""

def run(blast, hostfile, id_from_stitle, host_key, host_col, mge_type, ntot, label):
    meta=pd.read_csv(f"{B}/linked_eces_meta.tsv",sep="\t")
    sub=meta[meta.MGE_type==mge_type].copy()
    inf=dict(zip(sub.MGE, sub.host_genus))
    b=pd.read_csv(blast,sep="\t",header=None,
        names=["q","stitle","pident","length","qlen","slen","qcovs","evalue","bitscore"])
    b["refid"]=b["stitle"].map(id_from_stitle)
    # stream host table for the referenced ids
    want=set(b["refid"]); hmap={}
    with open(hostfile) as fh:
        rd=csv.reader(fh,delimiter="\t"); hdr=next(rd)
        ik=hdr.index(host_key); ih=hdr.index(host_col)
        for row in rd:
            if len(row)>ih and row[ik] in want: hmap[row[ik]]=row[ih]
    b["ref_g"]=b["refid"].map(hmap).map(genus)
    b["strict"]=(b.pident>=PID)&(b.qcovs>=QCOV)
    rows=[]
    for q,g in b.groupby("q"):
        g=g.sort_values("bitscore",ascending=False)
        known=bool(g.iloc[0].strict)                    # best hit is strict
        if not known:
            continue                                    # -> counted as "no strict match"
        gs=g[g.strict & (g.ref_g!="")]                  # strict hits carrying a genus host
        if len(gs):
            ig=inf.get(q,""); genera=set(gs.ref_g)
            supp = ig in genera                         # ANY strict hit's host genus matches inferred
            rg = ig if supp else gs.iloc[0].ref_g       # for plotting: matching genus, else best-hit genus
            cat=("host-supported (agrees)" if supp else "host mismatch")
            rows.append((q,known,rg,ig,cat))
        else:
            rows.append((q,known,"","","known (strict), host not genus-resolved"))
    seen={r[0] for r in rows}
    per=pd.DataFrame(rows,columns=["MGE","known_strict","ref_host_genus","inferred_host_genus","category"])
    per.to_csv(f"{B}/{label}_strict_validation.tsv",sep="\t",index=False)
    per.to_csv(f"{FIG}/{label}_strict_validation.tsv",sep="\t",index=False)
    n_known=int(per.known_strict.sum())
    res=per[per.category.isin(["host-supported (agrees)","host mismatch"])]
    supp=(res.category=="host-supported (agrees)").sum()
    n_nores=(per.category=="known (strict), host not genus-resolved").sum()
    n_nomatch=ntot-n_known
    counts={"host-supported (agrees)":supp,"host mismatch":len(res)-supp,
            "known (strict), host not genus-resolved":n_nores,"no strict match":n_nomatch}
    print(f"{label}: total={ntot} known(strict)={n_known} host-resolved={len(res)} supported={supp}"
          f" ({100*supp/len(res):.0f}% of resolved)")
    return counts

def uvig(s): return str(s).split("|")[0].split()[0]

vc=run(f"{B}/virus_imgvr_blast.tsv",
       "/shared/db/imgvr/v4_IMG_VR_2022-09-20_6/metadata/IMGVR_all_Host_information.tsv",
       uvig,"UVIG","Host taxonomy prediction","virus",116,"virus_imgvr")
pc=run(f"{B}/plasmid_imgpr_blast.tsv",
       "/shared/db/imgpr/2023-08-08_1/IMGPR_plasmid_data.tsv",
       uvig,"plasmid_id","host_taxonomy","plasmid",263,"plasmid_imgpr")
order=["host-supported (agrees)","host mismatch","known (strict), host not genus-resolved","no strict match"]
rows=[("Viruses (IMG/VR)",c,vc[c]) for c in order]+[("Plasmids (IMG/PR)",c,pc[c]) for c in order]
pd.DataFrame(rows,columns=["type","category","n"]).to_csv(f"{FIG}/ece_validation_summary_strict.csv",index=False)
print("wrote ece_validation_summary_strict.csv")
