#!/usr/bin/env python3
"""build_e2e_toy2.py — complexity-driven toy for the completeness gradient.

MODIFI drops contigs <5x, so we CANNOT starve genomes below 5x to force incompleteness.
Instead we raise COMMUNITY COMPLEXITY at low-but-usable depth: con-specific strains and
con-generic species tangle hifiasm_meta's overlap graph, so many genomes assemble
INCOMPLETELY even at 5-12x. That yields a real CheckM2 completeness spread while every
genome stays >=5x (MODIFI-usable).

Composition (all donors carry >=1 curated ECE):
  - multi-strain species (>=3 strains): 3 con-specific strains each  -> strain tangle
  - con-generic species (same genus): a few species per shared genus -> shared-region tangle
  - isolated distinct-genus species                                  -> clean, ~complete
Depth ladder low: 5-12x (all >= MODIFI's floor).
Emits toy.bam / toy.fq.gz / input_genomes/ / toy.manifest.csv (same layout as toy1).
"""
import os, sys, subprocess
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/simu_meta")
import build_community as bc

OUT = "/home/shuaiw/borg/paper/simu_meta_dir/C4_toy2"
SEED = 42
THREADS = 16
# low, all >=5x; assigned round-robin so complexity (not depth) drives the spread
DEPTHS = [5, 6, 7, 8, 6, 5, 8, 10, 12, 7]


def main():
    os.makedirs(OUT, exist_ok=True)
    gdir = os.path.join(OUT, "input_genomes"); os.makedirs(gdir, exist_ok=True)
    donors, bg, mge = bc.load_pool()
    ece_iso = set(mge["prefix"])
    d = donors[donors["Sample"].isin(ece_iso)].dropna(subset=["Average_DP", "Species"]).copy()
    d = d[(d["Species"].str.strip() != "") & (d["Average_DP"] >= 15)]  # need native>=~15 to hit 5-12x
    d["genus"] = d["Species"].str.split(" ").str[0]

    picked = []
    # 1. multi-strain species: 3 con-specific strains each (top 3 such species)
    st = d.groupby("Species")["Strain"].nunique()
    multi = st[st >= 3].index.tolist()[:3]
    for sp in multi:
        g = d[d["Species"] == sp].sort_values("Sample").groupby("Strain", as_index=False).first()
        picked.append(g.head(3).assign(role_grp="strain:" + sp))
    used_sp = set(multi)
    # 2. con-generic species: genera with >=2 UNused species -> 2 species each (top 2 genera)
    remain = d[~d["Species"].isin(used_sp)]
    ggen = remain.groupby("genus")["Species"].nunique()
    congen = ggen[ggen >= 2].index.tolist()[:2]
    for gn in congen:
        sps = remain[remain["genus"] == gn].sort_values("Sample").groupby("Species", as_index=False).first().head(2)
        picked.append(sps.assign(role_grp="genus:" + gn))
        used_sp |= set(sps["Species"])
    # 3. isolated distinct-genus species to fill toward ~25 (clean assemblers)
    used_gen = set(pd.concat(picked)["genus"])
    iso = (d[~d["genus"].isin(used_gen)].sort_values("Average_DP", ascending=False)
           .groupby("Species", as_index=False).first())
    iso = iso.groupby("genus", as_index=False).first().head(12)
    picked.append(iso.assign(role_grp="isolated"))

    pick = pd.concat(picked, ignore_index=True).drop_duplicates("Sample").reset_index(drop=True)
    pick["target_dp"] = [DEPTHS[i % len(DEPTHS)] for i in range(len(pick))]
    print(f"[toy2] {len(pick)} isolates: "
          f"{sum(pick.role_grp.str.startswith('strain'))} strain-tangle, "
          f"{sum(pick.role_grp.str.startswith('genus'))} genus-tangle, "
          f"{sum(pick.role_grp=='isolated')} isolated")

    prep_dir = os.path.join(OUT, "prepped_bams"); os.makedirs(prep_dir, exist_ok=True)
    bams, rows = [], []
    for _, r in pick.iterrows():
        s = r["Sample"]
        mb, downs = bc.prep_isolate_bam(s, r["Average_DP"], r["target_dp"],
                                        os.path.join(prep_dir, f"{s}.bam"), SEED, threads=THREADS)
        bams.append(mb)
        dst = os.path.join(gdir, f"{s}.fa")
        if not os.path.exists(dst): subprocess.run(["cp", r["genome"], dst], check=True)
        eces = sorted(mge[mge.prefix == s]["mge_contig"].unique())
        rows.append(dict(sample=s, species=r["Species"], genus=r["genus"], role_grp=r["role_grp"],
                         native_dp=r["Average_DP"], target_dp=r["target_dp"], downsampled=downs,
                         n_ece=len(eces), ece_contigs=";".join(eces), genome=dst))

    merged = os.path.join(OUT, "toy.bam")
    print(f"[toy2] merging {len(bams)} BAMs")
    subprocess.run(["samtools", "merge", "-f", "-c", "-p", "-@", str(THREADS), merged] + bams, check=True)
    subprocess.run(["samtools", "index", "-@", str(THREADS), merged], check=True)
    fq = os.path.join(OUT, "toy.fq.gz")
    with open(fq, "wb") as fo:
        p1 = subprocess.Popen(["samtools", "fastq", "-@", str(THREADS), merged], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip", "-c"], stdin=p1.stdout, stdout=fo)
        p1.stdout.close(); p2.communicate()
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "toy.manifest.csv"), index=False)
    for b in bams:
        if os.path.dirname(os.path.abspath(b)) == os.path.abspath(prep_dir):
            try: os.remove(b)
            except OSError: pass
    print(f"[toy2] done: {len(rows)} genomes, {fq}")


if __name__ == "__main__":
    main()
