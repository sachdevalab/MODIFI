#!/usr/bin/env python3
"""label_contig_origin.py — assign each assembled contig its true source isolate.

skani-dist each assembled contig against the known input isolate genomes; the best hit
(highest ANI, gated on aligned fraction) is the contig's origin. Chimeric contigs (a
strong second hit to a DIFFERENT isolate) are flagged. This is the de-novo ground truth
the concat-reference C1 gets for free from SRA-prefix contig names.

Usage: label_contig_origin.py <contigs.fa> <input_genomes_dir> <out.tsv> [threads]
Output TSV: contig  origin_sample  ani  align_frac  chimeric  second_sample  second_ani
"""
import sys, os, glob, subprocess, tempfile
import pandas as pd

SKANI = "/shared/software/bin/skani"
MIN_ANI = 95.0      # same-species/strain hit
MIN_AF = 15.0       # aligned fraction (%) of the (query) contig — skani reports ref/query AF


def split_contigs(fa, tmpdir):
    """Write each contig to its own fasta (skani query needs per-genome files for AF)."""
    paths = []
    import re
    name = None; seq = []
    def flush():
        if name:
            p = os.path.join(tmpdir, re.sub(r"[^A-Za-z0-9._-]", "_", name) + ".fa")
            with open(p, "w") as o:
                o.write(f">{name}\n" + "".join(seq) + "\n")
            paths.append((name, p))
    with open(fa) as f:
        for line in f:
            if line.startswith(">"):
                flush(); name = line[1:].strip().split()[0]; seq = []
            else:
                seq.append(line.strip())
        flush()
    return paths


def main():
    contigs, gdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
    threads = sys.argv[4] if len(sys.argv) > 4 else "16"
    refs = sorted(glob.glob(os.path.join(gdir, "*.fa")) + glob.glob(os.path.join(gdir, "*.fasta")))
    sample_of = {r: os.path.splitext(os.path.basename(r))[0] for r in refs}
    print(f"[label] {len(refs)} reference genomes; splitting contigs ...")

    with tempfile.TemporaryDirectory() as tmp:
        qpaths = split_contigs(contigs, tmp)
        qfiles = os.path.join(tmp, "qlist.txt"); rfiles = os.path.join(tmp, "rlist.txt")
        with open(qfiles, "w") as o: o.write("\n".join(p for _, p in qpaths) + "\n")
        with open(rfiles, "w") as o: o.write("\n".join(refs) + "\n")
        res = os.path.join(tmp, "skani.tsv")
        # skani dist: all query contigs vs all ref genomes
        subprocess.run([SKANI, "dist", "--ql", qfiles, "--rl", rfiles, "-o", res,
                        "-t", str(threads), "--min-af", "5", "-s", "90"], check=True)
        sk = pd.read_csv(res, sep="\t")
    # columns: Ref_file, Query_file, ANI, Align_fraction_ref, Align_fraction_query, Ref_name, Query_name
    qname = {p: n for n, p in qpaths}
    sk["contig"] = sk["Query_file"].map(qname)
    sk["origin_sample"] = sk["Ref_file"].map(sample_of)
    sk["af"] = sk[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1)
    sk = sk.sort_values(["contig", "ANI"], ascending=[True, False])

    rows = []
    for contig, g in sk.groupby("contig"):
        g = g[g["af"] >= MIN_AF]
        if g.empty:
            rows.append(dict(contig=contig, origin_sample="UNMAPPED", ani=0.0, align_frac=0.0,
                             chimeric=False, second_sample="", second_ani=0.0)); continue
        g = g.sort_values("ANI", ascending=False)
        best = g.iloc[0]
        others = g[g["origin_sample"] != best["origin_sample"]]
        sec = others.iloc[0] if len(others) else None
        chimeric = bool(sec is not None and sec["ANI"] >= MIN_ANI and sec["af"] >= 30)
        rows.append(dict(contig=contig, origin_sample=best["origin_sample"],
                         ani=round(best["ANI"], 2), align_frac=round(best["af"], 1),
                         chimeric=chimeric,
                         second_sample=(sec["origin_sample"] if sec is not None else ""),
                         second_ani=(round(sec["ANI"], 2) if sec is not None else 0.0)))
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    print(f"[label] {len(df)} contigs labeled ({(df.origin_sample=='UNMAPPED').sum()} unmapped, "
          f"{df.chimeric.sum()} chimeric) -> {out}")


if __name__ == "__main__":
    main()
