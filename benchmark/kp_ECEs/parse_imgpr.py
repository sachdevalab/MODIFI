#!/usr/bin/env python3
"""Compare linked plasmids to IMG/PR (metagenomic plasmid catalog): best megablast hit -> IMG/PR
host_taxonomy -> genus, vs our modification-inferred host genus. Complements the clinical-biased
mob_suite/RefSeq known-plasmid result with environmental-plasmid coverage."""
import os, re, csv
import pandas as pd

BASE="/home/shuaiw/borg/revision/linked_eces"
FIG="/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
PDATA="/shared/db/imgpr/2023-08-08_1/IMGPR_plasmid_data.tsv"
STRICT=dict(pid=90.0,qcov=50.0); REL=dict(pid=70.0,qcov=30.0)

def genus(lin):
    for tok in str(lin).split(";"):
        t=tok.strip()
        if t.startswith("g__"): return re.sub(r"_[A-Z]+$","",t[3:])
    return ""
def pid_of(stitle): return str(stitle).split("|")[0].split()[0]

def main():
    meta=pd.read_csv(f"{BASE}/linked_eces_meta.tsv",sep="\t")
    pl=meta[meta.MGE_type=="plasmid"].copy(); pl["inferred_host_genus"]=pl["host_genus"]
    bl=pd.read_csv(f"{BASE}/plasmid_imgpr_blast.tsv",sep="\t",header=None,
                   names=["q","stitle","pident","length","qlen","slen","qcovs","evalue","bitscore"])
    bl["pid_ref"]=bl["stitle"].map(pid_of)
    # stream IMG/PR metadata, keep only hit plasmid_ids
    want=set(bl["pid_ref"]); hmap={}
    with open(PDATA) as fh:
        rd=csv.reader(fh,delimiter="\t"); hdr=next(rd)
        ip=hdr.index("plasmid_id"); ih=hdr.index("host_taxonomy")
        for row in rd:
            if len(row)>ih and row[ip] in want: hmap[row[ip]]=row[ih]
    bl["host_lin"]=bl["pid_ref"].map(hmap); bl["ref_host_genus"]=bl["host_lin"].map(genus)
    best=bl.sort_values("bitscore",ascending=False).drop_duplicates("q").set_index("q")
    hosted=bl[bl.ref_host_genus!=""].sort_values("bitscore",ascending=False).drop_duplicates("q").set_index("q")

    rows=[]
    for _,r in pl.iterrows():
        q=r["MGE"]; row={"MGE":q,"environment":r["environment"],"mge_len":r["mge_len"],
                         "inferred_host_genus":r["inferred_host_genus"],
                         "imgpr_hit": q in best.index}
        if q in hosted.index:
            h=hosted.loc[q]; strict=(h.pident>=STRICT["pid"])and(h.qcovs>=STRICT["qcov"])
            rel=(h.pident>=REL["pid"])and(h.qcovs>=REL["qcov"]); rg=h.ref_host_genus
            row.update(host_ref=h.pid_ref,host_pident=round(h.pident,1),host_qcov=round(h.qcovs,1),
                       ref_host_genus=rg,tier=("strict" if strict else("related" if rel else "weak")),
                       ref_matches_inferred=("yes" if rg==r["inferred_host_genus"] else "no"))
        else: row["tier"]="none"
        rows.append(row)
    out=pd.DataFrame(rows)
    out.to_csv(f"{BASE}/plasmid_imgpr_validation.tsv",sep="\t",index=False)
    out.to_csv(f"{FIG}/plasmid_imgpr_validation.tsv",sep="\t",index=False)
    print(f"[imgpr] {len(out)} linked plasmids; with IMG/PR hit: {int(out.imgpr_hit.sum())}")
    for tier in ("strict","related"):
        sub=out[out.tier.isin(["strict"] if tier=="strict" else ["strict","related"])]
        sub=sub[sub.ref_matches_inferred.isin(["yes","no"])]
        if len(sub):
            y=(sub.ref_matches_inferred=="yes").sum()
            print(f"[imgpr] genus-resolved host check at {tier}: {len(sub)}; agreement {y}/{len(sub)} ({100*y/len(sub):.0f}%)")
    rel=out[out.tier.isin(["strict","related"]) & out.ref_matches_inferred.isin(["yes","no"])]
    if len(rel):
        print("[imgpr] agreement by environment (strict+related):")
        for env,g in rel.groupby("environment"):
            print(f"    {env}: {(g.ref_matches_inferred=='yes').sum()}/{len(g)}")

if __name__=="__main__": main()
