#!/usr/bin/env python3
"""Direct RefSeq mapping of the K. pneumoniae ECEs (new-network set), same criteria as the IMG search
(anicalc weighted-ANI pid>=95 & query-cov qcov>=85). Plasmids vs PLSDB (RefSeq plasmids), viruses vs
RefSeq viral. Writes a per-ECE hit table + a per-cluster summary, and cross-checks the mob_suite Mash
closest-reference. No figure."""
import csv
import pandas as pd

D = "/home/shuaiw/borg/revision/kp_eces/refseq_map"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/kp"
PLSDB_META = "/groups/diamond/databases/plasmid/PLSDB/2023_11_23_v2/plsdb.tsv"
VREF_META = "/home/shuaiw/borg/revision/linked_eces/viral_ref/viral_ref_host.tsv"
MOBTYPER = "/home/shuaiw/borg/revision/kp_eces/annotate/mob_typer.tsv"
ANI, QCOV = 95.0, 85.0

def acc(s):
    return str(s).split("|")[0].split()[0]

def stream_lookup(path, key_col, val_cols, want):
    """Return {key: {val_col: value}} for keys in `want`, streaming a big TSV."""
    out = {}
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
        ik = hdr.index(key_col); iv = {c: hdr.index(c) for c in val_cols}
        for row in rd:
            if len(row) > max(iv.values()) and row[ik] in want:
                out[row[ik]] = {c: row[iv[c]] for c in val_cols}
    return out

AMRFINDER = "/home/shuaiw/borg/revision/kp_eces/annotate/amrfinder.tsv"
PLASMIDFINDER = "/home/shuaiw/borg/revision/kp_eces/annotate/plasmidfinder.tsv"

def load_plasmidfinder():
    """abricate/PlasmidFinder per-contig incompatibility (Inc) / replicon groups."""
    import os, re
    from collections import defaultdict
    p = defaultdict(list)
    if not os.path.exists(PLASMIDFINDER):
        return p
    df = pd.read_csv(PLASMIDFINDER, sep="\t", dtype=str).fillna("")
    seqc = "SEQUENCE" if "SEQUENCE" in df.columns else df.columns[1]
    for _, r in df.iterrows():
        g = re.sub(r"_\d+.*$", "", r.get("GENE", ""))  # IncFIB(K)_1_Kpn3 -> IncFIB(K)
        if g and g not in p[r[seqc]]:
            p[r[seqc]].append(g)
    return p

def load_mob():
    """mob_typer per-contig: mobility + replicon/relaxase/MPF types + nearest RefSeq plasmid (Mash)."""
    import os
    m = {}
    if not os.path.exists(MOBTYPER):
        return m
    for r in csv.DictReader(open(MOBTYPER), delimiter="\t"):
        m[r["sample_id"]] = dict(
            mobility=r.get("predicted_mobility", ""),
            rep_types=r.get("rep_type(s)", ""), relaxase=r.get("relaxase_type(s)", ""),
            mpf=r.get("mpf_type", ""), orit=r.get("orit_type(s)", ""),
            mash_ref=r.get("mash_nearest_neighbor", ""), mash_dist=r.get("mash_neighbor_distance", ""),
            mash_id=r.get("mash_neighbor_identification", ""),
            host_range=r.get("predicted_host_range_overall_name", ""))
    return m

def load_amr():
    """AMRFinderPlus per-contig: gene symbols split by Element type (AMR / STRESS / VIRULENCE)."""
    import os
    from collections import defaultdict
    a = defaultdict(lambda: {"AMR": [], "STRESS": [], "VIRULENCE": []})
    if not os.path.exists(AMRFINDER):
        return a
    df = pd.read_csv(AMRFINDER, sep="\t", dtype=str).fillna("")
    cc = [c for c in df.columns if "ontig" in c][0]
    for _, r in df.iterrows():
        et = r.get("Element type", "").upper()
        g = r.get("Gene symbol", "")
        if et in ("AMR", "STRESS", "VIRULENCE") and g:
            a[r[cc]][et].append(g)
    return a

meta = pd.read_csv(f"{D}/kp_members_meta.tsv", sep="\t")
inf = dict(zip(meta.MGE, meta.inferred_host_genus.fillna("")))
clu = dict(zip(meta.MGE, meta.MGE_cluster))
linked = dict(zip(meta.MGE, meta.is_kp_linked))
mob = load_mob()
amr = load_amr()
pf = load_plasmidfinder()
def uniq(xs):
    seen = []
    for x in xs:
        if x and x not in seen: seen.append(x)
    return ";".join(seen)

def run(ani_tsv, db):
    """Return per-ECE best-strict-hit dict for one catalogue."""
    b = pd.read_csv(ani_tsv, sep="\t")           # qname,tname,num_alns,pid,qcov,tcov
    if not len(b):
        return {}
    b["refid"] = b["tname"].map(acc)
    want = set(b["refid"])
    if db == "plsdb":
        md = stream_lookup(PLSDB_META, "NUCCORE_ACC",
                           ["NUCCORE_Description", "TAXONOMY_genus", "TAXONOMY_species"], want)
        b["ref_desc"] = b["refid"].map(lambda a: md.get(a, {}).get("NUCCORE_Description", ""))
        b["ref_genus"] = b["refid"].map(lambda a: md.get(a, {}).get("TAXONOMY_genus", ""))
    else:
        md = stream_lookup(VREF_META, "accession", ["organism", "host_genus"], want)
        b["ref_desc"] = b["refid"].map(lambda a: md.get(a, {}).get("organism", ""))
        b["ref_genus"] = b["refid"].map(lambda a: md.get(a, {}).get("host_genus", ""))
    b["strict"] = (b.pid >= ANI) & (b.qcov >= QCOV)
    res = {}
    for q, g in b.groupby("qname"):
        gs = g[g.strict].sort_values(["pid", "qcov"], ascending=False)
        genera = sorted(set(x for x in gs.ref_genus if x))
        if len(gs):
            top = gs.iloc[0]
            res[q] = dict(known=1, ref_acc=top.refid, ref_desc=top.ref_desc, ani=top.pid,
                          qcov=top.qcov, tcov=top.tcov, ref_genus=top.ref_genus,
                          n_strict=len(gs), n_ref_genera=len(genera), ref_genera=";".join(genera))
        else:
            res[q] = dict(known=0)
    return res

pl = run(f"{D}/kp_plasmid_plsdb_ani.tsv", "plsdb")
vi = run(f"{D}/kp_virus_refseqviral_ani.tsv", "viral")

# ---- per-ECE table ----
rows = []
for _, m in meta.iterrows():
    q = m.MGE; r = (pl if m.MGE_type == "plasmid" else vi).get(q, {"known": 0})
    ig = inf.get(q, "")
    hs = int(bool(r.get("known")) and bool(ig) and ig in str(r.get("ref_genera", "")).split(";"))
    mn = mob.get(q, {})
    ag = amr.get(q, {"AMR": [], "STRESS": [], "VIRULENCE": []})
    rows.append(dict(MGE=q, cluster=m.MGE_cluster, type=m.MGE_type, is_kp_linked=m.is_kp_linked,
                     inferred_host_genus=ig, known=r.get("known", 0),
                     ref_acc=r.get("ref_acc", ""), ref_desc=r.get("ref_desc", ""),
                     ani=r.get("ani", ""), qcov=r.get("qcov", ""), tcov=r.get("tcov", ""),
                     ref_host_genus=r.get("ref_genus", ""), n_ref_genera=r.get("n_ref_genera", ""),
                     ref_genera=r.get("ref_genera", ""), host_supported=hs,
                     inc_groups=";".join(pf.get(q, [])),
                     mobility=mn.get("mobility", ""), rep_types=mn.get("rep_types", ""),
                     relaxase_MOB=mn.get("relaxase", ""), mpf_type=mn.get("mpf", ""),
                     amr_genes=uniq(ag["AMR"]), stress_genes=uniq(ag["STRESS"]),
                     virulence_genes=uniq(ag["VIRULENCE"]),
                     mob_host_range=mn.get("host_range", ""), mob_mash_ref=mn.get("mash_ref", ""),
                     mob_mash_dist=mn.get("mash_dist", ""), mob_mash_id=mn.get("mash_id", "")))
per = pd.DataFrame(rows).sort_values(["type", "cluster", "known"], ascending=[True, True, False])
per.to_csv(f"{D}/kp_refseq_hits.csv", index=False)
per.to_csv(f"{FIG}/kp_refseq_hits.csv", index=False)

# ---- per-cluster summary (14) ----
crows = []
for c, g in per.groupby("cluster"):
    known = g[g.known == 1]
    best = known.sort_values(["ani", "qcov"], ascending=False)
    b0 = best.iloc[0] if len(best) else None
    genera = sorted(set(x for s in known.ref_genera for x in str(s).split(";") if x))
    def uagg(col):
        return uniq([x for s in g[col] for x in str(s).split(";") if x])
    # cluster-level mobility: use the representative (cluster-id) member's call; fragmented members
    # lose T4SS/relaxase genes and get downgraded, so the representative is the most reliable.
    reprow = g[g.MGE == c]
    rep_mob = reprow.mobility.iloc[0] if len(reprow) else ""
    n_conj = int((g.mobility == "conjugative").sum())
    crows.append(dict(cluster=c, type=g.type.iloc[0], n_members=len(g),
                      n_kp_linked=int(g.is_kp_linked.sum()), n_known=len(known),
                      best_ref_acc=(b0.ref_acc if b0 is not None else ""),
                      best_ref_desc=(b0.ref_desc if b0 is not None else ""),
                      best_ani=(b0.ani if b0 is not None else ""),
                      best_qcov=(b0.qcov if b0 is not None else ""),
                      ref_host_genus=(b0.ref_host_genus if b0 is not None else ""),
                      n_ref_genera=len(genera), ref_genera=";".join(genera),
                      rep_mobility=rep_mob, n_conjugative_members=n_conj,
                      inc_groups=uagg("inc_groups"),
                      mobility=uagg("mobility"), rep_types=uagg("rep_types"),
                      relaxase_MOB=uagg("relaxase_MOB"), mpf_type=uagg("mpf_type"),
                      amr_genes=uagg("amr_genes"), stress_genes=uagg("stress_genes"),
                      virulence_genes=uagg("virulence_genes"),
                      mob_mash_ref=g.mob_mash_ref.iloc[0], mob_mash_dist=g.mob_mash_dist.iloc[0],
                      mob_mash_id=g.mob_mash_id.iloc[0], mob_host_range=g.mob_host_range.iloc[0]))
csum = pd.DataFrame(crows).sort_values(["n_known", "best_ani"], ascending=False)
csum.to_csv(f"{D}/kp_clusters_refseq_summary.csv", index=False)
csum.to_csv(f"{FIG}/kp_clusters_refseq_summary.csv", index=False)

# ---- report ----
plmem = set(meta[meta.MGE_type == "plasmid"].MGE); vimem = set(meta[meta.MGE_type == "virus"].MGE)
npl_known = sum(1 for q, r in pl.items() if q in plmem and r.get("known"))
nvi_known = sum(1 for q, r in vi.items() if q in vimem and r.get("known"))
print(f"K. pneumoniae ECE members: {len(meta)} (plasmid {int((meta.MGE_type=='plasmid').sum())}, "
      f"virus {int((meta.MGE_type=='virus').sum())})")
print(f"KNOWN (strict RefSeq hit): plasmids {npl_known}/{int((meta.MGE_type=='plasmid').sum())} vs PLSDB; "
      f"viruses {nvi_known}/{int((meta.MGE_type=='virus').sum())} vs RefSeq viral")
print(f"clusters with a known member: {int((csum.n_known>0).sum())}/{len(csum)}")
print("\n=== per-cluster summary ===")
cols = ["cluster", "type", "n_members", "n_known", "best_ani", "n_ref_genera",
        "mobility", "inc_groups", "relaxase_MOB", "amr_genes"]
with pd.option_context("display.max_columns", None, "display.width", 240, "display.max_colwidth", 42):
    print(csum[cols].to_string(index=False))
print(f"\nwrote kp_refseq_hits.csv ({len(per)} ECEs) + kp_clusters_refseq_summary.csv ({len(csum)} clusters)")
