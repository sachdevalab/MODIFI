#!/usr/bin/env python3
"""genomad_to_mge.py — convert geNomad plasmid+virus summaries to a MODIFI --mge_file.

Output TSV (3-col, the format build_community writes and MODIFI reads):
    seq_name    mge_type    length
Usage: genomad_to_mge.py <genomad_summary_dir> <out_mge.tsv>
  where <genomad_summary_dir> holds *_plasmid_summary.tsv and *_virus_summary.tsv
"""
import sys, os, glob, re
import pandas as pd


def read(path, mge_type):
    if not path or not os.path.exists(path):
        return []
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        name = str(r["seq_name"])
        if re.search(r"\|provirus", name) or name == "seq_name":
            continue
        out.append((name, mge_type, int(r["length"])))
    return out


def main():
    sdir, out = sys.argv[1], sys.argv[2]
    pl = glob.glob(os.path.join(sdir, "*_plasmid_summary.tsv"))
    vi = glob.glob(os.path.join(sdir, "*_virus_summary.tsv"))
    rows = []
    for p in pl:
        rows += read(p, "plasmid")
    for v in vi:
        rows += read(v, "virus")
    df = pd.DataFrame(rows, columns=["seq_name", "mge_type", "length"]).drop_duplicates("seq_name")
    df.to_csv(out, sep="\t", index=False)
    print(f"[genomad_to_mge] {len(df)} ECEs ({(df.mge_type=='plasmid').sum()} plasmid, "
          f"{(df.mge_type=='virus').sum()} virus) -> {out}")


if __name__ == "__main__":
    main()
