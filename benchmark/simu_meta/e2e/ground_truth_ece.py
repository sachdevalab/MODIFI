#!/usr/bin/env python3
"""ground_truth_ece.py — select ECE contigs by mapping to known curated ECEs (no geNomad).

Each input isolate's curated ECE contigs (from toy.manifest.csv `ece_contigs`) are pulled
from its input genome; skani maps every ASSEMBLED contig against these known ECE
sequences. An assembled contig is called an ECE when it maps to a curated ECE (ANI>=95
and the contig is mostly that ECE, query aligned-fraction >= 50%). Its type (plasmid/virus)
is inherited from the matched curated ECE (via the pool MGE table).

Writes toy.mge.tsv (seq_name  mge_type  length) — the MODIFI --mge_file — and
ece_truth.tsv (assembled ECE contig -> matched curated ECE, origin isolate, type).

Usage: ground_truth_ece.py <C4_toy_dir>
"""
import os, sys, glob, subprocess, tempfile, re
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/simu_meta")
import build_community as bc

SKANI = "/shared/software/bin/skani"
MIN_ANI, MIN_QAF = 95.0, 50.0


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].strip().split()[0], []
        else:
            seq.append(line.strip())
    if name: yield name, "".join(seq)


def main():
    D = sys.argv[1]
    man = pd.read_csv(os.path.join(D, "toy.manifest.csv"))
    _, _, mge = bc.load_pool()
    ece_type = mge.drop_duplicates("mge_contig").set_index("mge_contig")["mge_type"].to_dict()
    contigs = os.path.join(D, "toy.contigs.fa")
    fai = pd.read_csv(contigs + ".fai", sep="\t", header=None)
    clen = dict(zip(fai[0], fai[1]))

    with tempfile.TemporaryDirectory() as tmp:
        # 1. curated ECE reference sequences, one file per ECE contig
        refdir = os.path.join(tmp, "eces"); os.makedirs(refdir)
        ref_info = {}   # ref_file -> (mge_contig, origin_sample, type)
        for _, r in man.iterrows():
            if not isinstance(r.get("ece_contigs"), str) or not r["ece_contigs"]:
                continue
            genome = dict(read_fasta(r["genome"]))
            for ec in r["ece_contigs"].split(";"):
                if ec in genome:
                    p = os.path.join(refdir, re.sub(r"[^A-Za-z0-9._-]", "_", ec) + ".fa")
                    with open(p, "w") as o: o.write(f">{ec}\n{genome[ec]}\n")
                    ref_info[p] = (ec, r["sample"], ece_type.get(ec, "plasmid"))
        # 2. assembled contigs, one file per contig (query)
        qdir = os.path.join(tmp, "asm"); os.makedirs(qdir)
        qmap = {}
        for n, s in read_fasta(contigs):
            p = os.path.join(qdir, re.sub(r"[^A-Za-z0-9._-]", "_", n) + ".fa")
            with open(p, "w") as o: o.write(f">{n}\n{s}\n")
            qmap[p] = n
        ql = os.path.join(tmp, "ql"); rl = os.path.join(tmp, "rl")
        open(ql, "w").write("\n".join(qmap) + "\n")
        open(rl, "w").write("\n".join(ref_info) + "\n")
        res = os.path.join(tmp, "sk.tsv")
        subprocess.run([SKANI, "dist", "--ql", ql, "--rl", rl, "-o", res, "-t", "16",
                        "--min-af", "5", "-s", "85"], check=True)
        sk = pd.read_csv(res, sep="\t")

    sk["contig"] = sk["Query_file"].map(qmap)
    sk["ref_mge"] = sk["Ref_file"].map(lambda f: ref_info[f][0])
    sk["origin"] = sk["Ref_file"].map(lambda f: ref_info[f][1])
    sk["type"] = sk["Ref_file"].map(lambda f: ref_info[f][2])
    hit = sk[(sk["ANI"] >= MIN_ANI) & (sk["Align_fraction_query"] >= MIN_QAF)]
    hit = hit.sort_values("ANI", ascending=False).drop_duplicates("contig")

    mge_rows, truth = [], []
    for _, r in hit.iterrows():
        mge_rows.append((r["contig"], r["type"], clen.get(r["contig"], 0)))
        truth.append(dict(ece_contig=r["contig"], matched_curated_ece=r["ref_mge"],
                          origin_sample=r["origin"], mge_type=r["type"],
                          ani=round(r["ANI"], 2), qaf=round(r["Align_fraction_query"], 1)))
    mdf = pd.DataFrame(mge_rows, columns=["seq_name", "mge_type", "length"]).drop_duplicates("seq_name")
    mdf.to_csv(os.path.join(D, "toy.mge.tsv"), sep="\t", index=False)
    pd.DataFrame(truth).to_csv(os.path.join(D, "ece_truth.tsv"), sep="\t", index=False)
    print(f"[gt_ece] {len(mdf)} assembled ECE contigs "
          f"({(mdf.mge_type=='plasmid').sum()} plasmid, {(mdf.mge_type=='virus').sum()} virus) "
          f"from {len(ref_info)} curated ECEs -> toy.mge.tsv")


if __name__ == "__main__":
    main()
