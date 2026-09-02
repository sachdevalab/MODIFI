#!/usr/bin/env python3
"""Representative x every-Pfam-domain gene-count matrix for the 14 Kp cluster reps, from the eggNOG
PFAMs column. Rows = reps (ordered by length), columns = all distinct Pfam domains."""
import pandas as pd, re
GP="/home/shuaiw/borg/revision/kp_eces/gene_profile"
FIG="/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
A=f"{GP}/kp14.emapper.annotations"

rows=[l.rstrip("\n").split("\t") for l in open(A) if not l.startswith("##")]
hdr=[c.lstrip("#") for c in rows[0]]; df=pd.DataFrame(rows[1:],columns=hdr)
df=df[df["query"].str.len()>0]
df["rep"]=df["query"].map(lambda q: q.rsplit("_",1)[0])
rep=pd.read_csv(f"{GP}/kp14_representatives.tsv",sep="\t").sort_values("rep_length",ascending=False)
order=list(rep["representative"])
df=df[df["rep"].isin(order)]

recs=[]
for _,r in df.iterrows():
    pf=str(r.get("PFAMs","-") or "-")
    if pf in ("","-"): continue
    for p in pf.split(","):
        p=p.strip()
        if p and p!="-": recs.append((r["rep"], p))
long=pd.DataFrame(recs, columns=["rep","pfam"])
mat=long.groupby(["rep","pfam"]).size().unstack(fill_value=0).reindex(index=order, fill_value=0)
# order columns by total prevalence (desc) so the heatmap reads left->right by abundance
mat=mat.reindex(columns=mat.sum().sort_values(ascending=False).index)
lab={r["representative"]: f'{r["cluster"]}  ({r["rep_length"]/1000:.0f} kb)' for _,r in rep.iterrows()}
mat.index=[lab[r] for r in mat.index]
mat.to_csv(f"{FIG}/kp14_pfam_matrix.csv")
long.to_csv(f"{FIG}/kp14_pfam_long_sourcedata.csv", index=False)
print(f"[pfam] matrix {mat.shape[0]} reps x {mat.shape[1]} Pfam domains; nonzero cells={int((mat>0).sum().sum())}")
print(f"[pfam] Pfam-bearing ORFs total: {len(long)}")
print("[pfam] top 15 Pfams:", list(mat.columns[:15]))

# --- Pfam -> curated functional family map (majority family of the ORFs carrying each Pfam) ---
import importlib.util, collections
spec=importlib.util.spec_from_file_location("gm","/home/shuaiw/MODIFI/benchmark/kp_ECEs/kp14_gene_matrix.py")
gm=importlib.util.module_from_spec(spec)
# load only the FAMILIES + classify_family without running main()
src=open("/home/shuaiw/MODIFI/benchmark/kp_ECEs/kp14_gene_matrix.py").read()
ns={}; exec("import re\n"+src[src.index("FAMILIES = ["):src.index("\n\ndef main")], ns)
FAMILIES=ns["FAMILIES"]; classify=ns["classify_family"]
# re-read annotations to get per-ORF text + COG for family classification
rows2=[l.rstrip("\n").split("\t") for l in open(A) if not l.startswith("##")]
h2=[c.lstrip("#") for c in rows2[0]]; d2=pd.DataFrame(rows2[1:],columns=h2)
d2["text"]=(d2.get("Description","").astype(str)+" | "+d2.get("Preferred_name","").astype(str)
            +" | "+d2.get("PFAMs","").astype(str)+" | "+d2.get("KEGG_ko","").astype(str))
d2["family"]=[classify(t,c) for t,c in zip(d2["text"], d2.get("COG_category",""))]
pf_fam=collections.defaultdict(collections.Counter)
for _,r in d2.iterrows():
    pf=str(r.get("PFAMs","-") or "-")
    if pf in ("","-"): continue
    for p in pf.split(","):
        p=p.strip()
        if p and p!="-": pf_fam[p][r["family"]]+=1
pres=(mat>0).sum(0)
fam_rows=[{"pfam":p,"family":pf_fam[p].most_common(1)[0][0] if p in pf_fam else "Other/hypothetical",
           "n_reps":int(pres.get(p,0))} for p in mat.columns]
pd.DataFrame(fam_rows).to_csv(f"{FIG}/kp14_pfam_family_map.csv",index=False)
print(f"[pfam] wrote family map for {len(fam_rows)} Pfams")
