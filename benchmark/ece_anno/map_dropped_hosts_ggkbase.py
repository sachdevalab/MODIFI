#!/usr/bin/env python3
"""Map soil_s3_2 dropped (GTDB-unclassified) host contigs -> ggKbase scaffold names, from the
minimap2 PAF (host contigs vs ggKbase SR-VP_..._90cm metaMDBG assembly). Reports per host contig:
best ggKbase scaffold + all scaffolds covering >=20% of the host, with identity/coverage + bin."""
import csv, re, os, subprocess
from collections import defaultdict

O = "/home/shuaiw/borg/revision/ggkbase_map"
PAF = f"{O}/minimap2.paf"
GGK = "/home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI.contigs.fa"
GTDB = {  # from earlier lookup
 'soil_s3_2_581_L':'Unclassified Bacteria','soil_s3_2_534_L':'Unclassified Bacteria',
 'soil_s3_2_1215_L':'Unclassified Bacteria','soil_s3_2_5179_L':'Unclassified Bacteria',
 'soil_s3_2_3231_L':'Unclassified Bacteria','soil_s3_2_7923_L':'Unclassified Bacteria',
 'soil_s3_2_730_L':'Unclassified Bacteria','soil_s3_2_1142_L':'Unclassified Bacteria',
 'soil_s3_2_10279_L':'Unclassified Bacteria','soil_s3_2_6009_L':'Unclassified Bacteria',
 'soil_s3_2_4167_L':'Unclassified','soil_s3_2_1188_L':'Unclassified Bacteria',
 'soil_s3_2_14992_L':'Unclassified Bacteria'}

# parse PAF -> per (query,target) best alignment aggregate
rows = []
for ln in open(PAF):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 12: continue
    q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
    t, tlen, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
    matches, alnlen = int(f[9]), int(f[10])
    rows.append(dict(q=q, qlen=qlen, qcov=(qe-qs)/qlen, t=t, tlen=tlen,
                     tcov=(te-ts)/tlen, ident=matches/alnlen if alnlen else 0,
                     matches=matches, qspan=qe-qs))

# ggKbase header -> bin/id, only for targets we hit
targets = set(r["t"] for r in rows)
bin_of, id_of = {}, {}
for ln in subprocess.run(["grep", "^>", GGK], capture_output=True, text=True).stdout.splitlines():
    name = ln[1:].split()[0]
    if name not in targets: continue
    m = re.search(r'bin="([^"]*)"', ln); bin_of[name] = m.group(1) if m else ""
    m = re.search(r'id=(\S+)', ln); id_of[name] = m.group(1) if m else ""

# aggregate per query x target (a host may hit a scaffold in multiple blocks)
agg = defaultdict(lambda: dict(matches=0, qspan=0, qlen=0, tlen=0, tcov=0, ident_num=0, ident_den=0))
for r in rows:
    a = agg[(r["q"], r["t"])]
    a["matches"] += r["matches"]; a["qspan"] += r["qspan"]; a["qlen"] = r["qlen"]
    a["tlen"] = r["tlen"]; a["tcov"] = max(a["tcov"], r["tcov"])
    a["ident_num"] += r["matches"]; a["ident_den"] += r["matches"]/r["ident"] if r["ident"] else 0

by_q = defaultdict(list)
for (q, t), a in agg.items():
    qcov = a["qspan"]/a["qlen"] if a["qlen"] else 0
    ident = a["ident_num"]/a["ident_den"] if a["ident_den"] else 0
    by_q[q].append(dict(t=t, qcov=qcov, tcov=a["tcov"], ident=ident, matches=a["matches"],
                        tlen=a["tlen"], bin=bin_of.get(t, ""), id=id_of.get(t, "")))

out = []
for q in GTDB:
    hits = sorted(by_q.get(q, []), key=lambda h: -h["matches"])
    if not hits:
        out.append(dict(run2_contig=q, gtdb=GTDB[q], ggkbase_scaffold="NO_HIT", ggkbase_bin="",
                        ggkbase_id="", identity="", query_cov="", target_cov="", n_scaffolds=0)); continue
    best = hits[0]
    covering = [h for h in hits if h["qcov"] >= 0.20]
    total_qcov = min(1.0, sum(h["qcov"] for h in covering))
    out.append(dict(run2_contig=q, gtdb=GTDB[q],
                    ggkbase_scaffold=best["t"], ggkbase_bin=best["bin"], ggkbase_id=best["id"],
                    identity=round(best["ident"], 4), query_cov=round(best["qcov"], 3),
                    target_cov=round(best["tcov"], 3), n_scaffolds=len(covering),
                    total_query_cov=round(total_qcov, 3),
                    other_scaffolds=";".join(h["t"] for h in covering[1:])))

cols = ["run2_contig","gtdb","ggkbase_scaffold","ggkbase_bin","ggkbase_id","identity",
        "query_cov","target_cov","n_scaffolds","total_query_cov","other_scaffolds"]
with open(f"{O}/soil_s3_2_dropped_host_ggkbase_map.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
    for r in out: w.writerow(r)

print(f"{'run2_contig':20} {'best ggKbase scaffold':52} {'idn':>5} {'qcov':>5} {'#sc':>3} {'totcov':>6}")
for r in out:
    print(f"{r['run2_contig']:20} {r['ggkbase_scaffold']:52} {r.get('identity',''):>5} {r.get('query_cov',''):>5} {r.get('n_scaffolds',''):>3} {r.get('total_query_cov',''):>6}")
print(f"\nwrote {O}/soil_s3_2_dropped_host_ggkbase_map.csv")
