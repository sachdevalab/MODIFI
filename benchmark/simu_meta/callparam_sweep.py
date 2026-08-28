#!/usr/bin/env python3
"""callparam_sweep.py — assess --min_score / --min_cov by re-running ONLY the downstream steps
(split from the retained aligned BAM, then motif/profile/merge/host), reusing the retained
per-contig gffs. No re-alignment, no IPD recompute. Writes to a fresh isolated out_dir and
never overwrites an existing host_summary.csv.

Retained inputs (per orphan_300 rep, under modifi/<label>/): gffs/{ctg}.gff (per-locus IPD
calls with score in col 6 and a coverage= attribute), align_bam/aligned.bam. `split` with
--aligned_bam regenerates contigs/ + bams/ + mean_depth.csv from that aligned BAM (no re-align);
motif/profile read the staged gffs + regenerated contigs; host runs the linkage.

  --min_score V : symlink the retained gffs, run main.py --min_score V (motif/profile filter
                  gff loci by score column 6).
  --min_cov  V  : min_cov is baked into the gff at the compare stage (loci with coverage < min_cov
                  were dropped; the run used the minimum, 1). We emulate a STRICTER floor by
                  dropping gff loci whose coverage= attribute < V, then run at default min_score.
                  (Looser than 1 is unreachable.)
"""
import argparse
import glob
import os
import subprocess
import sys

MODIFI = "/home/shuaiw/MODIFI/main.py"
PY = "/home/shuaiw/miniconda3/envs/modifi/bin/python"


def _cov(line):
    for kv in line.rstrip("\n").split("\t")[-1].split(";"):
        if kv.startswith("coverage="):
            try:
                return int(kv.split("=")[1])
            except ValueError:
                return None
    return None


def stage_gffs(src_gffs, dst_gffs, min_cov):
    """Symlink retained per-contig gffs into dst; for min_cov>1 write coverage-filtered copies."""
    os.makedirs(dst_gffs, exist_ok=True)
    for g in glob.glob(os.path.join(src_gffs, "*.gff")):
        if g.endswith(".reprocess.gff"):
            continue
        dst = os.path.join(dst_gffs, os.path.basename(g))
        if os.path.exists(dst):
            continue
        if min_cov <= 1:
            os.symlink(g, dst)
        else:
            with open(g) as fi, open(dst, "w") as fo:
                for line in fi:
                    if line.startswith("#"):
                        fo.write(line)
                        continue
                    c = _cov(line)
                    if c is None or c >= min_cov:
                        fo.write(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="retained modifi/<label> dir (cached gffs + align_bam)")
    ap.add_argument("--ref", required=True, help="community ref .fa (with .fai)")
    ap.add_argument("--mge_file", required=True)
    ap.add_argument("--out_dir", required=True, help="fresh isolated work dir for this (param,value)")
    ap.add_argument("--min_score", type=int, default=30)
    ap.add_argument("--min_cov", type=int, default=1)
    ap.add_argument("--threads", type=int, default=64)
    a = ap.parse_args()

    out_hs = os.path.join(a.out_dir, "host_summary.csv")
    if os.path.exists(out_hs):
        print(f"[skip] {out_hs} already exists, not overwriting.")
        return

    aligned = os.path.join(a.src, "align_bam", "aligned.bam")
    if not os.path.exists(aligned):
        sys.exit(f"[error] missing aligned bam: {aligned}")

    os.makedirs(a.out_dir, exist_ok=True)
    stage_gffs(os.path.join(a.src, "gffs"), os.path.join(a.out_dir, "gffs"), a.min_cov)

    cmd = [PY, MODIFI, "--aligned_bam", aligned, "-r", a.ref, "-o", a.out_dir,
           "--mge_file", a.mge_file,
           "--run_steps", "split", "motif", "profile", "merge", "host",
           "--min_score", str(a.min_score), "--no-clean", "--threads", str(a.threads)]
    print(f"[run] min_score={a.min_score} min_cov={a.min_cov} -> {a.out_dir}")
    print("      " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    print(f"[done] {out_hs}" if os.path.exists(out_hs) else f"[warn] no host_summary at {out_hs}")


if __name__ == "__main__":
    main()
