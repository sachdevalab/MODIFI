#!/usr/bin/env python
"""
Per-gene MTase summary TSV for the Klebsiella pneumoniae 248_1 cluster.

One row per predicted MTase gene (MT/IIG) in the 15 MAGs, with:
  contig, gene, gene_annotation (REBASE homolog), HMM, methylation_type,
  motif (REBASE recognition motif), motif_detected (is that motif detected as
  modified in the contig), matched_detected_motif (evidence).

Reuses the parsing/matching logic in build_mtase_matrix.py.
"""

import os
import pandas as pd

from build_mtase_matrix import (
    MAGS, OUT_DIR, ACTIVITY_FRACTION,
    load_mtases, load_detected_motifs, motifs_equivalent, clean_systype,
)


def main():
    rows = []
    for mag in MAGS:
        detected = load_detected_motifs(mag)
        for _, r in load_mtases(mag).iterrows():
            motif = r["Homolog motif"]
            if motif:
                hits = sorted({f"{dm}({frac:.2f})" for dm, frac in detected
                               if frac >= ACTIVITY_FRACTION and motifs_equivalent(motif, dm)})
                detected_flag = "yes" if hits else "no"
                matched = ";".join(hits)
            else:
                detected_flag = "NA"      # no REBASE motif to assess
                matched = ""
            rows.append({
                "contig": mag,
                "gene": r["Gene"],
                "gene_annotation": r["REBASE homolog"] or "NA",
                "HMM": r["HMM"],
                "methylation_type": r["Predicted methylation"] or "NA",
                "motif": motif or "NA",
                "motif_detected": detected_flag,
                "matched_detected_motif": matched,
                "system_type": clean_systype(r["System Type"]),
            })

    df = pd.DataFrame(rows).sort_values(["contig", "HMM", "gene"])
    out = os.path.join(OUT_DIR, "Kp248_1_mtase_summary.tsv")
    df.to_csv(out, sep="\t", index=False)
    print(f"{len(df)} MTase genes across {df['contig'].nunique()} contigs\n")
    print(df.to_string(index=False))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
