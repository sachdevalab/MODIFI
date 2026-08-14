#!/usr/bin/env python3
"""frag_ref.py — fragment ladder_58's host genomes into synthetic contigs (Part H synthetic).

Randomly chop every HOST contig into ~FRAG-bp pieces (ECE contigs kept WHOLE — they are the
query, not the MAG). Methylation is then called ONCE on this fragmented reference; afterwards
MAG completeness (retain a subset of a host's fragments) and contamination (add foreign
fragments to a host's bin) are applied purely via --bin_file at linkage time. So this heavy
step runs once; the completeness x contamination grid is cheap linkage-only re-runs.

Emits <OUT>/ladder58_frag.ref.fa and <OUT>/fragment_map.tsv (frag, orig, isolate, len).
Usage: frag_ref.py [frag_bp]   (default 500000)
"""
import os, sys
import pandas as pd

SRC = "/home/shuaiw/borg/paper/simu_meta_dir/C1/ladder_58_rep2"
OUT = "/home/shuaiw/borg/paper/simu_meta_dir/C1/ladder58_frag"
FRAG = int(sys.argv[1]) if len(sys.argv) > 1 else 500_000


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
    os.makedirs(OUT, exist_ok=True)
    ece = set(pd.read_csv(f"{SRC}/ladder_58_rep2.mge_list.tsv", sep="\t")["seq_name"])
    out_ref = f"{OUT}/ladder58_frag.ref.fa"
    rows, n_ece, n_frag = [], 0, 0
    with open(out_ref, "w") as o:
        for name, seq in read_fasta(f"{SRC}/ladder_58_rep2.ref.fa"):
            if name in ece:                                  # ECE: keep whole (query)
                o.write(f">{name}\n{seq}\n"); n_ece += 1
                continue
            iso = name.split("_")[0]                         # host contig -> fragments
            L = len(seq)
            for i, s in enumerate(range(0, L, FRAG)):
                fn = f"{name}_f{i}"
                frag = seq[s:s + FRAG]
                o.write(f">{fn}\n{frag}\n")
                rows.append(dict(frag=fn, orig=name, isolate=iso, length=len(frag)))
                n_frag += 1
    pd.DataFrame(rows).to_csv(f"{OUT}/fragment_map.tsv", sep="\t", index=False)
    print(f"[frag_ref] {n_ece} ECEs kept whole; {n_frag} host fragments "
          f"({rows.__len__() and int(pd.DataFrame(rows).groupby('isolate').size().mean())} avg/isolate) "
          f"@ {FRAG//1000}kb -> {out_ref}")


if __name__ == "__main__":
    main()
