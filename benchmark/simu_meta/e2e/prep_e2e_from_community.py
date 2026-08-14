#!/usr/bin/env python3
"""prep_e2e_from_community.py — set up an end-to-end (de-novo) working dir from an existing
reference-based community (bg_80/150/300, one replicate).

Reuses the community's merged CCS BAM as the read source; run_e2e.sh then assembles it
de novo and links ECEs contig-level (same pipeline validated on the toys). Creates:
  <OUT>/toy.bam            symlink -> the community merged BAM (reads, kinetics preserved)
  <OUT>/input_genomes/     symlinks -> every isolate assembly (skani ground-truth origin)
  <OUT>/toy.manifest.csv   sample, species, role, target_dp, genome, ece_contigs, n_ece

Usage: prep_e2e_from_community.py <label>   e.g. bg_80  (community must be under C1/<label>/)
"""
import os, sys
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/simu_meta")
import build_community as bc

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir"


def main():
    label = sys.argv[1]
    SRC = f"{ROOT}/C1/{label}"
    OUT = f"{ROOT}/C4_{label}"
    gdir = f"{OUT}/input_genomes"
    os.makedirs(gdir, exist_ok=True)

    man = pd.read_csv(f"{SRC}/{label}.manifest.csv")
    _, _, mge = bc.load_pool()
    ece_by = mge.groupby("prefix")["mge_contig"].apply(lambda s: ";".join(sorted(set(s))))

    # reads symlink
    src_bam = f"{SRC}/{label}.bam"
    dst_bam = f"{OUT}/toy.bam"
    if not os.path.exists(dst_bam):
        os.symlink(src_bam, dst_bam)

    rows = []
    for _, r in man.iterrows():
        s = r["sample"]
        g = r["genome"]
        dst = f"{gdir}/{s}.fa"
        if not os.path.exists(dst) and isinstance(g, str) and os.path.exists(g):
            os.symlink(g, dst)
        eces = ece_by.get(s, "") if str(r.get("role")) == "donor" else ""
        rows.append(dict(sample=s, species=r.get("species"), role=r.get("role"),
                         target_dp=r.get("target_dp"), genome=dst, ece_contigs=eces,
                         n_ece=(len(eces.split(";")) if eces else 0)))
    out_man = pd.DataFrame(rows)
    out_man.to_csv(f"{OUT}/toy.manifest.csv", index=False)
    n_donor_ece = int((out_man["n_ece"] > 0).sum())
    print(f"[prep_e2e] {OUT}: {len(out_man)} isolates ({(out_man.role=='donor').sum()} donor / "
          f"{(out_man.role=='background').sum()} background); {n_donor_ece} with ECEs, "
          f"{int(out_man.n_ece.sum())} curated ECEs total; reads -> {src_bam}")


if __name__ == "__main__":
    main()
