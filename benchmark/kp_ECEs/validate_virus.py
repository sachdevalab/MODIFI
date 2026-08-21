#!/usr/bin/env python3
"""Virus-linkage validation: match each linked virus to its closest host-annotated RefSeq viral
genome (blastn best hit), and compare that reference's documented host to our modification-inferred
host. Also folds in the CRISPR-spacer confirmation already computed in the repo."""
import os, re
import pandas as pd

BASE="/home/shuaiw/borg/revision/linked_eces"
FIG="/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
CRISPR="/home/shuaiw/MODIFI/benchmark/spacer/crispr_validation_breakdown.tsv"
# "known virus" thresholds (nucleotide): a real, close relative
MIN_PID=90.0; MIN_QCOV=50.0

def genus(s):
    s=re.sub(r"^(uncultured|Candidatus)\s+","",str(s)).strip()
    return re.sub(r"_[A-Z]+$","",s.split(" ")[0]) if s else ""

def main():
    meta=pd.read_csv(f"{BASE}/linked_eces_meta.tsv",sep="\t")
    vir=meta[meta.MGE_type=="virus"].copy()
    vir["inferred_host_genus"]=vir["host_genus"]

    bl=pd.read_csv(f"{BASE}/virus_blast.tsv",sep="\t",header=None,
                   names=["q","s","pident","length","qlen","slen","qcovs","evalue","bitscore"])
    host=pd.read_csv(f"{BASE}/viral_ref/viral_ref_host.tsv",sep="\t")
    hmap=dict(zip(host.accession,host.host)); hgmap=dict(zip(host.accession,host.host_genus))

    # best hit per virus by bitscore
    bl=bl.sort_values("bitscore",ascending=False).drop_duplicates("q")
    bl["ref_host"]=bl["s"].map(hmap); bl["ref_host_genus"]=bl["s"].map(hgmap)
    best=bl.set_index("q")

    rows=[]
    for _,r in vir.iterrows():
        q=r["MGE"]; b=best.loc[q] if q in best.index else None
        if b is None:
            rows.append({**base_row(r), "has_hit":False}); continue
        known = (b.pident>=MIN_PID) and (b.qcovs>=MIN_QCOV)
        rg=genus(b.ref_host_genus)
        match = "" if not known else ("yes" if rg and rg==r["inferred_host_genus"] else "no")
        rows.append({**base_row(r), "has_hit":True,
                     "closest_ref":b.s, "pident":round(b.pident,1), "qcov":round(b.qcovs,1),
                     "ref_host":b.ref_host, "ref_host_genus":rg,
                     "known_virus":known, "ref_matches_inferred_genus":match})
    out=pd.DataFrame(rows)

    # fold in CRISPR spacer confirmation
    cr=pd.read_csv(CRISPR,sep="\t"); cr=cr[cr.type=="virus"][["MGE","validated"]]
    out=out.merge(cr.rename(columns={"validated":"crispr_validated"}),on="MGE",how="left")
    out.to_csv(f"{BASE}/virus_validation.tsv",sep="\t",index=False)
    out.to_csv(f"{FIG}/virus_validation.tsv",sep="\t",index=False)

    nv=len(out); nhit=int(out.has_hit.sum()); nknown=int(out.get("known_virus",pd.Series(dtype=bool)).sum())
    kn=out[out.known_virus==True]
    y=(kn.ref_matches_inferred_genus=="yes").sum(); n=(kn.ref_matches_inferred_genus=="no").sum()
    print(f"[virus] {nv} linked viruses; {nhit} have a RefSeq-viral hit; "
          f"{nknown} match a KNOWN virus (pident>={MIN_PID}, qcov>={MIN_QCOV}%)")
    if y+n:
        print(f"[virus] known-vs-inferred host agreement (genus): {y}/{y+n} ({100*y/(y+n):.0f}%)")
    print(f"[virus] CRISPR-spacer confirmed links: {int(out.crispr_validated.fillna(0).sum())}/{nv}")
    if nknown:
        print("\n[virus] known viruses (ref host vs inferred host):")
        cols=["MGE","environment","inferred_host_genus","closest_ref","pident","qcov","ref_host_genus","ref_matches_inferred_genus"]
        print(kn[cols].to_string(index=False))

def base_row(r):
    return {"MGE":r["MGE"],"environment":r["environment"],"mge_len":r["mge_len"],
            "inferred_host_species":r["host_species"],"inferred_host_genus":r["inferred_host_genus"]}

if __name__=="__main__":
    main()
