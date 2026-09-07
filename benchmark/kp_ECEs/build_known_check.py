#!/usr/bin/env python3
"""Known-ECE verification: for each of the 14 Kp cluster representatives, pair it with its best
type-appropriate database match (plasmid -> PLSDB, virus -> RefSeq viral), report identity/coverage +
each side's length and GC, and extract the representative and reference FASTAs for manual re-mapping.
"""
import os, csv, subprocess, glob
import pandas as pd

GP = "/home/shuaiw/borg/revision/kp_eces/gene_profile"
OUT = "/home/shuaiw/borg/revision/kp_eces/known_check"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
PLSDB_FNA = "/groups/diamond/databases/plasmid/PLSDB/2023_11_23_v2/plsdb.fna"
VIRAL_DB = "/shared/db/refseq_viral/viral"
REPS_FNA = f"{GP}/kp14_reps.fna"
os.makedirs(f"{OUT}/reps", exist_ok=True)
os.makedirs(f"{OUT}/refs", exist_ok=True)


def sh(cmd):
    return subprocess.run(cmd, shell=True, text=True, capture_output=True)


def seqkit_stats(fna):
    """return {name: (length, gc)} via seqkit fx2tab -n -l -g"""
    r = sh(f"conda run -n tldr seqkit fx2tab -n -l -g {fna}")
    d = {}
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) >= 3:
            d[p[0].split()[0]] = (int(float(p[1])), round(float(p[2]), 2))
    return d


reps = pd.read_csv(f"{GP}/kp14_representatives.tsv", sep="\t")
repmeta = {r["representative"]: (r["cluster"], r["rep_type"], int(r["rep_length"]))
           for _, r in reps.iterrows()}
rep_gc = seqkit_stats(REPS_FNA)

# PLSDB hits per rep
plsdb = {r["representative"]: r for r in csv.DictReader(open(f"{GP}/kp14_ani_known.tsv"), delimiter="\t")}

# best viral hit per virus rep (max aln_frac_smaller then pid)
viral = {}
for r in csv.DictReader(open(f"{GP}/kp14_viral_ani.tsv"), delimiter="\t"):
    q = r["qname"]; afs = max(float(r["qcov"]), float(r["tcov"])); pid = float(r["pid"])
    key = (afs, pid)
    if q not in viral or key > viral[q][0]:
        viral[q] = (key, r)

# assemble rows (type-appropriate DB; phage-plasmid gets an extra cross-DB PLSDB row)
def known(ani, afs):
    return bool(float(ani) >= 95 and float(afs) >= 85)

rows = []          # dicts
plsdb_accs = set()
viral_accs = set()
for _, rp in reps.iterrows():
    rep = rp["representative"]; typ = rp["rep_type"]; clu = rp["cluster"]
    if typ == "virus" and rep in viral:
        _, v = viral[rep]
        afs = round(max(float(v["qcov"]), float(v["tcov"])), 1)
        rows.append(dict(cluster=clu, representative=rep, rep_type=typ, search_db="RefSeq viral",
                         ANI=round(float(v["pid"]), 2), qcov=round(float(v["qcov"]), 1),
                         tcov=round(float(v["tcov"]), 1), aln_frac_smaller=afs,
                         known=known(v["pid"], afs), db_ece_id=v["tname"].split()[0],
                         ref_description=""))
        viral_accs.add(v["tname"].split()[0])
        # phage-plasmid cross-DB PLSDB row
        if rep in plsdb and rep == "infant_25_51_C":
            p = plsdb[rep]
            rows.append(dict(cluster=clu, representative=rep, rep_type=typ,
                             search_db="PLSDB (cross-type)", ANI=float(p["ANI"]),
                             qcov=float(p["qcov"]), tcov=float(p["tcov"]),
                             aln_frac_smaller=float(p["aln_frac_smaller"]),
                             known=known(p["ANI"], p["aln_frac_smaller"]),
                             db_ece_id=p["best_ref"], ref_description=p["ref_description"]))
            plsdb_accs.add(p["best_ref"])
    else:
        p = plsdb[rep]
        rows.append(dict(cluster=clu, representative=rep, rep_type=typ, search_db="PLSDB",
                         ANI=float(p["ANI"]), qcov=float(p["qcov"]), tcov=float(p["tcov"]),
                         aln_frac_smaller=float(p["aln_frac_smaller"]),
                         known=known(p["ANI"], p["aln_frac_smaller"]),
                         db_ece_id=p["best_ref"], ref_description=p["ref_description"]))
        plsdb_accs.add(p["best_ref"])

# --- extract reference FASTAs ---
# PLSDB: one seqkit grep pass
pat = "|".join(sorted(a for a in plsdb_accs if a))
sh(f'conda run -n tldr seqkit grep -n -r -p "{pat}" {PLSDB_FNA} -o {OUT}/refs/_plsdb_refs.fna')
# viral: blastdbcmd per accession
for acc in sorted(viral_accs):
    if acc:
        sh(f'conda run -n mod blastdbcmd -db {VIRAL_DB} -entry {acc} -out {OUT}/refs/_viral_{acc}.fna')

# build acc -> single-record fasta (rename per rep) and gather ref stats
def extract_one(src, acc, dest):
    sh(f'conda run -n tldr seqkit grep -n -r -p "{acc}" {src} -o {dest}')

ref_stats = {}
for _, rp in reps.iterrows():
    rep = rp["representative"]
for row in rows:
    rep = row["representative"]; acc = row["db_ece_id"]
    dest = f"{OUT}/refs/{rep}__{acc}.fna"
    if row["search_db"].startswith("RefSeq viral"):
        src = f"{OUT}/refs/_viral_{acc}.fna"
    else:
        src = f"{OUT}/refs/_plsdb_refs.fna"
    extract_one(src, acc, dest)
    st = seqkit_stats(dest)
    if st:
        L, G = list(st.values())[0]
        ref_stats[(rep, acc)] = (L, G)

# copy rep fastas individually
for rep in repmeta:
    sh(f'conda run -n tldr seqkit grep -p "{rep}" {REPS_FNA} -o {OUT}/reps/{rep}.fna')

# finalize table
for row in rows:
    rep = row["representative"]
    clu, typ, L = repmeta[rep]
    row["rep_length"] = L
    row["rep_GC"] = rep_gc.get(rep, (0, 0))[1]
    rl, rg = ref_stats.get((rep, row["db_ece_id"]), ("", ""))
    row["db_ece_length"] = rl
    row["db_ece_GC"] = rg

cols = ["cluster", "representative", "rep_type", "rep_length", "rep_GC", "search_db", "ANI",
        "qcov", "tcov", "aln_frac_smaller", "known", "db_ece_id", "db_ece_length", "db_ece_GC",
        "ref_description"]
df = pd.DataFrame(rows)[cols].sort_values(["rep_length", "search_db"], ascending=[False, True])
df.to_csv(f"{OUT}/kp14_known_check.tsv", sep="\t", index=False)
df.to_csv(f"{FIG}/kp14_known_check.tsv", sep="\t", index=False)
print(f"reps: {df['representative'].nunique()}; rows: {len(df)}; known: {df['known'].sum()}")
print(df.to_string(index=False))
