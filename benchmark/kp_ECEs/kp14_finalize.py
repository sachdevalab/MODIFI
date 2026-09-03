#!/usr/bin/env python3
"""After run_kp14_reprofile.sh: (1) update kp14_ani_known.tsv for the changed representative(s),
(2) rebuild kp14_master_table.csv, all restricted to the current (strict) 14 representatives.
The 13 unchanged reps keep their existing ANI rows; only reps missing from the table (e.g. the new
infant_2_845_L) are recomputed from kp14_plsdb_ani.tsv + PLSDB metadata."""
import csv, os
import pandas as pd

GP = "/home/shuaiw/borg/revision/kp_eces/gene_profile"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
PLSDB_META = "/groups/diamond/databases/plasmid/PLSDB/2023_11_23_v2/plsdb.tsv"

reps = pd.read_csv(f"{GP}/kp14_representatives.tsv", sep="\t")
rep_set = list(reps["representative"])

# --- PLSDB metadata (acc -> description, species) ---
meta = {}
with open(PLSDB_META, encoding="utf-8", errors="replace") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    acc_c = "NUCCORE_ACC"; desc_c = "NUCCORE_Description"; sp_c = "TAXONOMY_species"
    for r in rd:
        meta[r[acc_c]] = (r.get(desc_c, ""), r.get(sp_c, ""))

# --- best PLSDB hit per rep from anicalc output ---
best = {}
ani_f = f"{GP}/kp14_plsdb_ani.tsv"
if os.path.exists(ani_f):
    a = pd.read_csv(ani_f, sep="\t")
    # anicalc columns: qname tname num_alns pid qcov tcov
    for q, sub in a.groupby("qname"):
        sub = sub.copy()
        sub["afs"] = sub[["qcov", "tcov"]].max(axis=1)
        sub = sub.sort_values(["pid", "afs"], ascending=False)
        top = sub.iloc[0]
        acc = str(top["tname"]).split()[0]
        desc, sp = meta.get(acc, ("", ""))
        known = bool(top["pid"] >= 95 and top["afs"] >= 85)
        best[q] = dict(best_ref=acc, db="PLSDB", ref_description=desc or acc, host=sp,
                       ANI=round(float(top["pid"]), 2), qcov=round(float(top["qcov"]), 1),
                       tcov=round(float(top["tcov"]), 1), aln_frac_smaller=round(float(top["afs"]), 1),
                       known_ANI=known)

# --- update kp14_ani_known.tsv: keep existing rows for reps still present, recompute missing ones ---
akf = f"{GP}/kp14_ani_known.tsv"
existing = {}
if os.path.exists(akf):
    for r in csv.DictReader(open(akf), delimiter="\t"):
        existing[r["representative"]] = r
cols = ["representative", "type", "best_ref", "db", "ref_description", "ANI", "qcov", "tcov",
        "aln_frac_smaller", "known_ANI"]
rep_type = dict(zip(reps["representative"], reps["rep_type"]))
rows = []
for rp in rep_set:
    if rp in existing:
        rows.append({c: existing[rp].get(c, "") for c in cols})
    elif rp in best:
        b = best[rp]
        rows.append({"representative": rp, "type": rep_type.get(rp, "plasmid"),
                     "best_ref": b["best_ref"], "db": b["db"], "ref_description": b["ref_description"],
                     "ANI": b["ANI"], "qcov": b["qcov"], "tcov": b["tcov"],
                     "aln_frac_smaller": b["aln_frac_smaller"], "known_ANI": b["known_ANI"]})
    else:  # no PLSDB hit at all
        rows.append({"representative": rp, "type": rep_type.get(rp, "plasmid"), "best_ref": "",
                     "db": "PLSDB", "ref_description": "", "ANI": 0, "qcov": 0, "tcov": 0,
                     "aln_frac_smaller": 0, "known_ANI": False})
ak = pd.DataFrame(rows, columns=cols)
ak.to_csv(akf, sep="\t", index=False)
ak.to_csv(f"{FIG}/kp14_ani_known.tsv", sep="\t", index=False)

# --- master table ---
ann = pd.read_csv(f"{GP}/kp14_annotation.tsv", sep="\t")
mob = dict(zip(ann["representative"], ann.get("predicted_mobility", pd.Series(dtype=str))))
knownmap = {r["representative"]: r for _, r in ak.iterrows()}
host_lookup = {rp: (meta.get(str(knownmap[rp]["best_ref"]).split()[0], ("", ""))[1]
                    if str(knownmap[rp]["best_ref"]) else "") for rp in rep_set}
mrows = []
for _, r in reps.iterrows():
    rp = r["representative"]; k = knownmap.get(rp, {})
    is_known = str(k.get("known_ANI", "")).lower() in ("true", "1")
    mrows.append({
        "cluster": r["cluster"], "cluster_members_n": r["n_members"], "representative": rp,
        "length_bp": r["rep_length"], "mobility": mob.get(rp, ""), "type": r["rep_type"],
        "circular": r["rep_circular"],
        "known_in_refseq": is_known,
        "refseq_genome": k.get("best_ref", "") if is_known else "",
        "refseq_genome_host": host_lookup.get(rp, "").replace("_", " ") if is_known else "",
    })
m = pd.DataFrame(mrows).sort_values("length_bp", ascending=False)
m.to_csv(f"{FIG}/kp14_master_table.csv", index=False)
print("known reps:", int(ak["known_ANI"].astype(str).str.lower().isin(["true", "1"]).sum()), "/", len(ak))
print(m.to_string(index=False))
