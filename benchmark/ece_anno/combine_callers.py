#!/usr/bin/env python3
"""Combine ECE calls from geNomad + VirSorter2 + VIBRANT + VirSorter1 + PlasX across the
64 metagenomes, restricted to >=5x-depth, >=1kb contigs. No in-house circular.

For each caller, prefer the freshly-run revision-tree output
(/home/shuaiw/borg/revision/ece_callers/<sample>/...) and fall back to the pre-existing
run2 output (/home/shuaiw/borg/paper/run2/<sample>/...). Writes a combined, deduped ECE
table with method provenance. Read-only w.r.t. caller outputs."""
import os, glob, re
import pandas as pd
from Bio import SeqIO

RUN2 = "/home/shuaiw/borg/paper/run2"
REV = "/home/shuaiw/borg/revision/ece_callers"
OUT = "/home/shuaiw/borg/revision/ece_anno/expanded"
PLASX_MIN = 0.5
# VirSorter1: keep all free-phage categories (1/2/3), exclude prophages -- handled by section
# tracking in vs1_calls() (prophages are integrated, i.e. not extrachromosomal).
os.makedirs(OUT, exist_ok=True)

samples = sorted(pd.read_csv(
    "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/genome_data_all_samples.csv")["sample"].unique())
env_map = dict(zip(*[pd.read_csv(
    "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/genome_data_all_samples.csv")[c]
    for c in ("sample", "environment")]))


def depth_map(s):
    for mm in ("methylation4", "methylation3"):
        p = f"{RUN2}/{s}/{s}_{mm}/mean_depth.csv"
        if os.path.exists(p):
            d = pd.read_csv(p)
            return dict(zip(d.contig, d.depth)), dict(zip(d.contig, d.length))
    return {}, {}


def genomad_calls(s):
    gp = f"{s}.hifiasm.p_ctg.rename"; out = {}
    for kind in ("plasmid", "virus"):
        p = f"{RUN2}/{s}/Genomad/{gp}_summary/{gp}_{kind}_summary.tsv"
        if not os.path.exists(p):
            continue
        for _, r in pd.read_csv(p, sep="\t").iterrows():
            name = str(r["seq_name"])
            if "|provirus" in name or name == "seq_name":
                continue
            out[name] = kind
    return out


def vs2_calls(s):
    p = f"{REV}/{s}/vs2_ok/final-viral-score.tsv"          # fresh run
    if not os.path.exists(p):
        p = f"{RUN2}/{s}/virsorter2/final-viral-score.tsv"  # existing
    out = set()
    if os.path.exists(p):
        for _, r in pd.read_csv(p, sep="\t").iterrows():
            f = str(r["seqname"]).split("||")
            if len(f) > 1 and "partial" in f[1]:   # skip provirus regions (contig||..partial)
                continue
            out.add(f[0])
    return out


def vibrant_calls(s):
    fs = glob.glob(f"{REV}/{s}/vibrant_out/**/*.phages_combined.fna", recursive=True) \
        or glob.glob(f"{RUN2}/{s}/vibrant/**/*.phages_combined.fna", recursive=True)
    out = set()
    if fs:
        fs.sort(key=len, reverse=True)
        for rec in SeqIO.parse(fs[0], "fasta"):
            fld = rec.id.split("_")
            if len(fld) >= 2 and fld[-2] == "fragment":
                continue
            out.add(rec.id)
    return out


def vs1_calls(s):
    # VIRSorter_global-phage-signal.csv has 6 sections by "## N -" header: sections 1-3 are
    # "Complete phage contigs" (FREE PHAGE), sections 4-6 are "Prophages" (INTEGRATED). The
    # Category column (col 5) resets to 1/2/3 within EACH section, so free-vs-prophage lives only
    # in the section header. Keep only free-phage contigs (all categories 1/2/3); exclude prophages.
    p = f"{REV}/{s}/vs1_out/VIRSorter_global-phage-signal.csv"
    out = set()
    if os.path.exists(p):
        in_phage = False
        for line in open(p):
            if line.startswith("##"):
                if "Complete phage" in line:
                    in_phage = True
                elif "Prophage" in line:
                    in_phage = False
                # the repeated "## Contig_id,..." header leaves in_phage unchanged
                continue
            if line.startswith("Contig_id") or not in_phage:
                continue
            c = line.split(",")
            if len(c) < 5:
                continue
            out.add(re.sub(r"^VIRSorter_", "", c[0]))
    return out


def plasx_calls(s):
    for p in (glob.glob(f"{REV}/{s}/{s}-plasx-scores.txt")
              + glob.glob(f"{RUN2}/{s}/*plasmid_nr-scores*")
              + glob.glob(f"{RUN2}/{s}/plasx/*scores*")):
        d = pd.read_csv(p, sep="\t", header=None, names=["c", "sc"])
        d["sc"] = pd.to_numeric(d["sc"], errors="coerce")
        return set(d.loc[d.sc >= PLASX_MIN, "c"])
    return set()


rows = []
cover = {k: 0 for k in ("genomad", "vs2", "vibrant", "vs1", "plasx")}
for s in samples:
    dm, lm = depth_map(s)
    if not dm:
        continue
    gm = genomad_calls(s); vs2 = vs2_calls(s); vib = vibrant_calls(s)
    vs1 = vs1_calls(s); px = plasx_calls(s)
    for k, v in (("genomad", gm), ("vs2", vs2), ("vibrant", vib), ("vs1", vs1), ("plasx", px)):
        cover[k] += bool(v)
    per = {}
    for c, k in gm.items():
        per.setdefault(c, set()).add("genomad_" + k)
    for c in vs2:
        per.setdefault(c, set()).add("virsorter2")
    for c in vib:
        per.setdefault(c, set()).add("vibrant")
    for c in vs1:
        per.setdefault(c, set()).add("virsorter1")
    for c in px:
        per.setdefault(c, set()).add("plasx")
    for c, m in per.items():
        dp = dm.get(c); ln = lm.get(c)
        if dp is None or dp < 5 or (ln or 0) < 1000:
            continue
        viral = any(x in m for x in ("virsorter2", "vibrant", "virsorter1", "genomad_virus"))
        plasmidic = any(x in m for x in ("plasx", "genomad_plasmid"))
        typ = ("virus" if "genomad_virus" in m else "plasmid" if "genomad_plasmid" in m
               else "virus" if viral else "plasmid")
        rows.append({"sample": s, "environment": env_map.get(s, "NA"), "MGE": c,
                     "MGE_type": typ, "mge_len": ln, "MGE_cov": dp,
                     "methods": ",".join(sorted(m)), "type_conflict": int(viral and plasmidic)})
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/combined_metagenome_eces.csv", index=False)
print("caller coverage (samples w/ any call):", cover, "of", len(samples))
print(f"\ncombined ECEs (>=5x,>=1kb): {len(df)}  "
      f"(plasmid {(df.MGE_type=='plasmid').sum()}, virus {(df.MGE_type=='virus').sum()})")
allm = df.methods.str.get_dummies(sep=",")
print("\nECEs carrying each method:"); print(allm.sum().sort_values(ascending=False).to_string())
has_gm = df.methods.str.contains("genomad")
print(f"\ngeNomad-called: {has_gm.sum()} | NON-geNomad (NEW): {(~has_gm).sum()}")
print(f"  new by type: plasmid {((~has_gm)&(df.MGE_type=='plasmid')).sum()}, "
      f"virus {((~has_gm)&(df.MGE_type=='virus')).sum()}")
print(f"\nwrote {OUT}/combined_metagenome_eces.csv")
