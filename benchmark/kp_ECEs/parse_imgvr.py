#!/usr/bin/env python3
"""Compare linked viruses to IMG/VR: best dc-megablast hit -> IMG/VR UViG host prediction ->
genus, compared to our modification-inferred host genus. IMG/VR covers uncultivated/environmental
phages that RefSeq lacks."""
import os, re
import pandas as pd

BASE="/home/shuaiw/borg/revision/linked_eces"
FIG="/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
HOSTINFO="/shared/db/imgvr/v4_IMG_VR_2022-09-20_6/metadata/IMGVR_all_Host_information.tsv"
# thresholds: "same/known virus" (strict) vs "related known virus" (looser)
STRICT=dict(pid=90.0,qcov=50.0)
REL=dict(pid=70.0,qcov=30.0)

def genus(lin):
    for tok in str(lin).split(";"):
        t=tok.strip()
        if t.startswith("g__"):
            return re.sub(r"_[A-Z]+$","",t[3:])
    return ""

def uvig(stitle):
    return str(stitle).split("|")[0].split()[0]

def main():
    meta=pd.read_csv(f"{BASE}/linked_eces_meta.tsv",sep="\t")
    vir=meta[meta.MGE_type=="virus"].copy()
    vir["inferred_host_genus"]=vir["host_genus"]

    bl=pd.read_csv(f"{BASE}/virus_imgvr_blast.tsv",sep="\t",header=None,
                   names=["q","stitle","pident","length","qlen","slen","qcovs","evalue","bitscore"])
    bl["uvig"]=bl["stitle"].map(uvig)
    # stream the huge IMG/VR host table, keeping only UViGs that are blast hits
    want=set(bl["uvig"]); hmap={}
    import csv as _csv
    with open(HOSTINFO) as _fh:
        rd=_csv.reader(_fh,delimiter="\t"); hdr=next(rd)
        iu=hdr.index("UVIG"); ih=hdr.index("Host taxonomy prediction")
        for row in rd:
            if len(row)>ih and row[iu] in want:
                hmap[row[iu]]=row[ih]
    bl["host_lin"]=bl["uvig"].map(hmap)
    bl["ref_host_genus"]=bl["host_lin"].map(genus)
    # best hit overall (by bitscore) and best hit that HAS a host prediction
    best=bl.sort_values("bitscore",ascending=False).drop_duplicates("q").set_index("q")
    hosted=bl[bl["ref_host_genus"]!=""].sort_values("bitscore",ascending=False).drop_duplicates("q").set_index("q")

    rows=[]
    for _,r in vir.iterrows():
        q=r["MGE"]; row={"MGE":q,"environment":r["environment"],"mge_len":r["mge_len"],
                         "inferred_host_genus":r["inferred_host_genus"]}
        if q in best.index:
            b=best.loc[q]
            row.update(imgvr_hit=b.uvig, pident=round(b.pident,1), qcov=round(b.qcovs,1),
                       imgvr_taxon_hit=True)
        else:
            row["imgvr_taxon_hit"]=False
        # host comparison uses the best HOST-annotated hit
        if q in hosted.index:
            h=hosted.loc[q]
            strict = (h.pident>=STRICT["pid"]) and (h.qcovs>=STRICT["qcov"])
            rel    = (h.pident>=REL["pid"]) and (h.qcovs>=REL["qcov"])
            rg=h.ref_host_genus
            row.update(host_uvig=h.uvig, host_pident=round(h.pident,1), host_qcov=round(h.qcovs,1),
                       ref_host_genus=rg, tier=("strict" if strict else ("related" if rel else "weak")),
                       ref_matches_inferred=("yes" if rg==r["inferred_host_genus"] else "no"))
        else:
            row["tier"]="none"
        rows.append(row)
    out=pd.DataFrame(rows)
    out.to_csv(f"{BASE}/virus_imgvr_validation.tsv",sep="\t",index=False)
    out.to_csv(f"{FIG}/virus_imgvr_validation.tsv",sep="\t",index=False)

    n=len(out)
    print(f"[imgvr] {n} linked viruses")
    print(f"[imgvr] with any IMG/VR hit: {int(out.imgvr_taxon_hit.sum())}")
    for tier in ("strict","related"):
        sub=out[out.tier.isin(["strict"] if tier=="strict" else ["strict","related"])]
        sub=sub[sub.ref_matches_inferred.isin(["yes","no"])]
        if len(sub):
            y=(sub.ref_matches_inferred=="yes").sum()
            print(f"[imgvr] host-annotated hit at {tier} tier: {len(sub)}; genus agreement {y}/{len(sub)} ({100*y/len(sub):.0f}%)")
    # breakdown by environment (related tier)
    rel=out[out.tier.isin(["strict","related"]) & out.ref_matches_inferred.isin(["yes","no"])]
    if len(rel):
        print("[imgvr] agreement by environment (strict+related):")
        for env,g in rel.groupby("environment"):
            y=(g.ref_matches_inferred=="yes").sum()
            print(f"    {env}: {y}/{len(g)}")

if __name__=="__main__":
    main()
