#!/usr/bin/env python3
"""build_e2e_toy.py — build the toy end-to-end community (Part H/C4 first test).

10 donor isolates, distinct species, each with >=1 curated ECE, reads mixed at a DEPTH
LADDER so de-novo assembly completeness spans a real range (low depth -> fragmented ->
low completeness). Emits:
  toy.bam        merged CCS reads (kinetics preserved) -> MODIFI
  toy.fq.gz      samtools fastq of toy.bam            -> hifiasm_meta
  input_genomes/ the 10 known source assemblies       -> skani ground truth
  toy.manifest.csv   isolate, species, target/native depth, ECE contigs
Reuses build_community helpers (prep_isolate_bam, load_pool, ccs_bam_path).
"""
import os, sys, subprocess
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/simu_meta")
import build_community as bc

OUT = "/home/shuaiw/borg/paper/simu_meta_dir/C4_toy"
DEPTH_LADDER = [4, 8, 12, 20, 30, 40, 50, 60, 80, 100]   # ~ per-isolate target depths
SEED = 42
THREADS = 16


def main():
    os.makedirs(OUT, exist_ok=True)
    gdir = os.path.join(OUT, "input_genomes"); os.makedirs(gdir, exist_ok=True)
    donors, bg, mge = bc.load_pool()
    ece_iso = set(mge["prefix"])
    d = donors[donors["Sample"].isin(ece_iso)].dropna(subset=["Average_DP", "Species"]).copy()
    d = d[d["Species"].str.strip() != ""]
    d = d[d["Average_DP"] >= 50].sort_values("Average_DP", ascending=False)
    pick = d.groupby("Species", as_index=False).first().head(len(DEPTH_LADDER)).reset_index(drop=True)
    # assign ascending depth to isolates (ladder) so completeness varies across genomes
    pick["target_dp"] = DEPTH_LADDER[:len(pick)]
    print(f"[toy] {len(pick)} isolates; depth ladder {DEPTH_LADDER}")

    prep_dir = os.path.join(OUT, "prepped_bams"); os.makedirs(prep_dir, exist_ok=True)
    bams, rows = [], []
    for _, r in pick.iterrows():
        s = r["Sample"]
        out_bam = os.path.join(prep_dir, f"{s}.bam")
        merge_bam, downs = bc.prep_isolate_bam(s, r["Average_DP"], r["target_dp"],
                                               out_bam, SEED, threads=THREADS)
        bams.append(merge_bam)
        # copy source genome for skani ground truth
        g = r["genome"]; dst = os.path.join(gdir, f"{s}.fa")
        if not os.path.exists(dst):
            subprocess.run(["cp", g, dst], check=True)
        eces = sorted(mge[mge.prefix == s]["mge_contig"].unique())
        rows.append(dict(sample=s, species=r["Species"], native_dp=r["Average_DP"],
                         target_dp=r["target_dp"], downsampled=downs,
                         n_ece=len(eces), ece_contigs=";".join(eces), genome=dst))
        print(f"  {s:14} {r['Species'][:28]:28} target={r['target_dp']:>3}x native={r['Average_DP']:.0f} ECEs={len(eces)}")

    merged = os.path.join(OUT, "toy.bam")
    print(f"[toy] merging {len(bams)} BAMs -> {merged}")
    subprocess.run(["samtools", "merge", "-f", "-c", "-p", "-@", str(THREADS), merged] + bams, check=True)
    subprocess.run(["samtools", "index", "-@", str(THREADS), merged], check=True)

    fq = os.path.join(OUT, "toy.fq.gz")
    print(f"[toy] samtools fastq -> {fq}")
    with open(fq, "wb") as fo:
        p1 = subprocess.Popen(["samtools", "fastq", "-@", str(THREADS), merged], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(["gzip", "-c"], stdin=p1.stdout, stdout=fo)
        p1.stdout.close(); p2.communicate()
        if p2.returncode != 0:
            raise RuntimeError("samtools fastq | gzip failed")

    man = pd.DataFrame(rows)
    man.to_csv(os.path.join(OUT, "toy.manifest.csv"), index=False)
    # clean prepped subsampled copies (keep merged)
    for b in bams:
        if os.path.dirname(os.path.abspath(b)) == os.path.abspath(prep_dir):
            try: os.remove(b)
            except OSError: pass
    print(f"[toy] done: {merged}, {fq}, {len(man)} genomes in {gdir}")


if __name__ == "__main__":
    main()
