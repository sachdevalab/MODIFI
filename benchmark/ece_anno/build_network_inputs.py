#!/usr/bin/env python3
"""Step 1 of the network rebuild: assemble revised per-sample inputs for plot_linkage_data.

For each sample writes under borg/revision/network/per_sample/<sample>/:
  - host_summary.csv : all host rows for that sample's revised-passing ECEs, taken from the
      best-available linkage source (meth_all > meth_missing > run2 geNomad > newpass). Same schema
      as MODIFI host_summary; the network applies its own final_score>0.8 & specificity<0.001 filter.
  - all_mge.tsv : seq_name<TAB>type for that sample's revised-passing ECEs (so MGE_type is correct
      for new-caller ECEs; only plasmid/virus, no 'novel').
Also writes borg/revision/network/passing_contigs_by_sample.tsv (sample<TAB>contig) for fasta extraction.
"""
import csv, glob, os
from collections import defaultdict

A = "/home/shuaiw/borg/revision/ece_anno"
L = "/home/shuaiw/borg/revision/ece_linkage"
RUN2 = "/home/shuaiw/borg/paper/run2"
OUT = "/home/shuaiw/borg/revision/network"
os.makedirs(f"{OUT}/per_sample", exist_ok=True)

# passing set: (sample,MGE) -> type ; also methods
passing = {}; methods = {}
for r in csv.DictReader(open(f"{A}/expanded/filterpass_revised_final.csv")):
    passing[(r["sample"], r["MGE"])] = r["MGE_type"]
    methods[(r["sample"], r["MGE"])] = r.get("methods", "")

# gather host_summary rows per source, keyed by sample -> MGE -> list(row dicts)
# priority: 0 meth_all/ocean chunks, 1 meth_missing, 2 run2 geNomad, 3 newpass meth
sources = []  # (priority, glob_or_paths, sample_from_path_fn)
def add_dir(priority, pattern, samplefn):
    sources.append((priority, pattern, samplefn))

meth_all = glob.glob(f"{L}/*/meth_all/host_summary.csv")
ocean_chunks = glob.glob(f"{L}/ocean_1/meth_all_c*/host_summary.csv")
meth_missing = glob.glob(f"{L}/*/meth_missing/host_summary.csv")
run2_geno = glob.glob(f"{RUN2}/*/*_methylation4/host_summary.csv")   # methy_v=4 (manuscript default)
newpass = glob.glob(f"{L}/*/meth/host_summary.csv")

def samp_ml(p): return p.split("/ece_linkage/")[1].split("/")[0]
def samp_run2(p): return p.split("/run2/")[1].split("/")[0]

# per sample -> per MGE -> {priority: [rows]}
store = defaultdict(lambda: defaultdict(dict))
HDR = None
def ingest(paths, priority, samplefn, ocean_force=None):
    global HDR
    for p in paths:
        s = ocean_force or samplefn(p)
        try:
            rows = list(csv.DictReader(open(p)))
        except Exception:
            continue
        if rows and HDR is None:
            HDR = list(rows[0].keys())
        for r in rows:
            key = (s, r["MGE"])
            if key not in passing:
                continue
            store[s][r["MGE"]].setdefault(priority, []).append(r)

ingest(meth_all, 0, samp_ml)
ingest(ocean_chunks, 0, samp_ml, ocean_force="ocean_1")
ingest(meth_missing, 1, samp_ml)
ingest(run2_geno, 2, samp_run2)
ingest(newpass, 3, samp_ml)

# write per-sample outputs
samples = sorted(set(s for (s, m) in passing))
n_rows_total = 0; n_mge_written = 0; n_mge_norow = 0
contigs_out = open(f"{OUT}/passing_contigs_by_sample.tsv", "w")
for s in samples:
    d = f"{OUT}/per_sample/{s}"; os.makedirs(d, exist_ok=True)
    hs_rows = []
    mge_types = []
    for (ss, mge), typ in passing.items():
        if ss != s:
            continue
        mge_types.append((mge, typ))
        contigs_out.write(f"{s}\t{mge}\n")
        byprio = store.get(s, {}).get(mge, {})
        if not byprio:
            n_mge_norow += 1
            continue
        prio = min(byprio.keys())
        hs_rows.extend(byprio[prio])
        n_mge_written += 1
    # host_summary.csv
    cols = HDR
    with open(f"{d}/host_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for r in hs_rows:
            w.writerow(r); n_rows_total += 1
    # all_mge.tsv
    with open(f"{d}/all_mge.tsv", "w") as fh:
        fh.write("seq_name\ttype\tlength\tmethods\tGTDB\n")
        for mge, typ in mge_types:
            fh.write(f"{mge}\t{typ}\t\t{methods[(s,mge)]}\t\n")
contigs_out.close()
print(f"samples: {len(samples)}")
print(f"passing MGEs with >=1 host row written: {n_mge_written}")
print(f"passing MGEs with NO host row (any source): {n_mge_norow}")
print(f"total host_summary rows written: {n_rows_total}")
print(f"contig list: {OUT}/passing_contigs_by_sample.tsv ({sum(1 for _ in open(f'{OUT}/passing_contigs_by_sample.tsv'))} contigs)")
