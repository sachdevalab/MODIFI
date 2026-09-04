#!/usr/bin/env python3
"""Summarise antiSMASH BGC regions for the 14 Kp ECE representatives.
Reads each antismash/<rep>/<rep>.json, extracts every detected region (product type(s), coordinates,
on-edge, best known-cluster (MIBiG) hit + similarity), joins cluster/type/length from
kp14_representatives.tsv, and writes a per-region TSV and a per-rep summary TSV.
"""
import os, json, glob
import pandas as pd

GP = "/home/shuaiw/borg/revision/kp_eces/gene_profile"
AS = "/home/shuaiw/borg/revision/kp_eces/antismash"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"

reps = pd.read_csv(f"{GP}/kp14_representatives.tsv", sep="\t")
meta = {r["representative"]: (r["cluster"], r["rep_type"], int(r["rep_length"]))
        for _, r in reps.iterrows()}

rows = []
for rep in reps["representative"]:
    jf = f"{AS}/{rep}/{rep}.json"
    if not os.path.exists(jf):
        rows.append(dict(representative=rep, region="NA", note="antismash not run/failed"))
        continue
    data = json.load(open(jf))
    for rec in data.get("records", []):
        contig = rec.get("id", rep)
        # regions live in areas / features of type 'region'
        for feat in rec.get("features", []):
            if feat.get("type") != "region":
                continue
            q = feat.get("qualifiers", {})
            products = q.get("product", [])
            loc = feat.get("location", "")
            # known-cluster comparison (clusterblast/knownclusterblast) is under modules; capture best hit
            rows.append(dict(
                representative=rep, contig=contig,
                region=(q.get("region_number", [""]) or [""])[0],
                products=";".join(products),
                location=loc,
                contig_edge=(q.get("contig_edge", [""]) or [""])[0],
            ))

# attach best knowncluster hit from the knownclusterblast module (if present)
def best_known(rep):
    jf = f"{AS}/{rep}/{rep}.json"
    if not os.path.exists(jf):
        return {}
    data = json.load(open(jf))
    out = {}
    for rec in data.get("records", []):
        mods = rec.get("modules", {})
        kcb = mods.get("antismash.modules.clusterblast", {}) or mods.get("antismash.modules.knownclusterblast", {})
        # knownclusterblast results: results -> [{ranking: [(ref, score), ...]}] per region
        kc = mods.get("antismash.modules.clusterblast", {})
        knn = kc.get("knowncluster", {}) if isinstance(kc, dict) else {}
        results = knn.get("results", []) if isinstance(knn, dict) else []
        for i, res in enumerate(results, 1):
            ranking = res.get("ranking", [])
            if ranking:
                ref = ranking[0][0]
                desc = ref.get("description", "") if isinstance(ref, dict) else str(ref)
                out[i] = desc
    return out

# best MIBiG known-cluster hit per (rep, region_number)
def known_map(rep):
    jf = f"{AS}/{rep}/{rep}.json"
    out = {}
    if not os.path.exists(jf):
        return out
    d = json.load(open(jf))
    for rec in d.get("records", []):
        kc = rec.get("modules", {}).get("antismash.modules.clusterblast", {})
        kn = kc.get("knowncluster", {}) if isinstance(kc, dict) else {}
        for i, res in enumerate(kn.get("results", []), 1):
            rk = res.get("ranking", [])
            if rk:
                ref, sc = rk[0][0], rk[0][1]
                desc = ref.get("description", "") if isinstance(ref, dict) else str(ref)
                acc = ref.get("accession", "") if isinstance(ref, dict) else ""
                hits = sc.get("hits", "") if isinstance(sc, dict) else ""
                out[str(i)] = f"{desc} ({acc}; {hits} gene hits)"
    return out

df = pd.DataFrame(rows)
if not df.empty and "products" in df.columns:
    df["cluster"] = df["representative"].map(lambda r: meta.get(r, ("", "", 0))[0])
    df["ece_type"] = df["representative"].map(lambda r: meta.get(r, ("", "", 0))[1])
    df["rep_length"] = df["representative"].map(lambda r: meta.get(r, ("", "", 0))[2])
    km = {rep: known_map(rep) for rep in reps["representative"]}
    df["best_known_cluster"] = df.apply(lambda x: km.get(x["representative"], {}).get(str(x["region"]), ""), axis=1)
    df = df[["cluster", "representative", "ece_type", "rep_length", "contig", "region",
             "products", "location", "contig_edge", "best_known_cluster"]]
    df.to_csv(f"{FIG}/kp14_antismash_bgc.tsv", sep="\t", index=False)
    df.to_csv(f"{AS}/kp14_antismash_bgc.tsv", sep="\t", index=False)

# per-rep summary
summ = []
for _, r in reps.iterrows():
    rep = r["representative"]
    sub = df[(df["representative"] == rep) & (df.get("products", "").astype(str) != "")] if not df.empty else pd.DataFrame()
    prods = sorted({p for row in sub.get("products", []) for p in str(row).split(";") if p})
    summ.append(dict(cluster=r["cluster"], representative=rep, ece_type=r["rep_type"],
                     rep_length=int(r["rep_length"]), n_BGC=len(sub),
                     BGC_types=";".join(prods)))
sdf = pd.DataFrame(summ).sort_values("rep_length", ascending=False)
sdf.to_csv(f"{FIG}/kp14_antismash_summary.tsv", sep="\t", index=False)
sdf.to_csv(f"{AS}/kp14_antismash_summary.tsv", sep="\t", index=False)

print(f"reps with >=1 BGC: {(sdf['n_BGC']>0).sum()}/{len(sdf)}; total regions: {int(sdf['n_BGC'].sum())}")
print(sdf.to_string(index=False))
