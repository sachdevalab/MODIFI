#!/usr/bin/env python3
"""
Find contigs across soil metagenomes that share >3 motifs with the Klingon genome.

Reference motifs: /home/shuaiw/borg/paper/natasha/soil_115/motifs/FINAL_Chon_Klingon.motifs.csv
Search path:      /home/shuaiw/borg/paper/gg_run3/soil_*/soil_*_methylation4/motifs/*.motifs.csv
"""

import glob
import os
import pandas as pd

KLINGON_MOTIFS_FILES = {
    "Chon_soil115": "/home/shuaiw/borg/paper/natasha/soil_115/motifs/FINAL_Chon_Klingon.motifs.csv",
    "Ghos_soil100": "/home/shuaiw/borg/paper/natasha/soil_100/motifs/FINAL_Ghos_Klingon.motifs.csv",
    "Chon_soil90":  "/home/shuaiw/borg/paper/natasha/soil_90/motifs/FINAL_Chon_Klingon.motifs.csv",
}
METHYLATION_GLOB = "/home/shuaiw/borg/paper/gg_run3/soil_*/soil_*_methylation4"
MIN_SHARED_MOTIFS = 2
OUTPUT_FILE = "contigs_sharing_gt2_klingon_motifs.tsv"

# Load all Klingon motif sets
klingon_motifs = {}
for name, path in KLINGON_MOTIFS_FILES.items():
    motifs = set(pd.read_csv(path)["motifString"].dropna().unique())
    klingon_motifs[name] = motifs
    print(f"{name}: {len(motifs)} unique motifs: {sorted(motifs)}")
print()

results = []

methylation_dirs = sorted(glob.glob(METHYLATION_GLOB))
print(f"Found {len(methylation_dirs)} methylation4 directories\n")

for meth_dir in methylation_dirs:
    sample = os.path.basename(meth_dir).replace("_methylation4", "")
    motifs_dir = os.path.join(meth_dir, "motifs")

    if not os.path.isdir(motifs_dir):
        print(f"  [SKIP] No motifs/ subfolder in {meth_dir}")
        continue

    contig_files = glob.glob(os.path.join(motifs_dir, "*.motifs.csv"))
    print(f"{sample}: {len(contig_files)} contigs")

    for cfile in contig_files:
        contig_name = os.path.basename(cfile).replace(".motifs.csv", "")
        try:
            df = pd.read_csv(cfile)
        except Exception as e:
            print(f"  [WARN] Could not read {cfile}: {e}")
            continue

        if "motifString" not in df.columns:
            continue

        contig_motifs = set(df["motifString"].dropna().unique())

        row = {"sample": sample, "contig": contig_name, "contig_total_motifs": len(contig_motifs)}
        max_shared = 0
        for name, ref_motifs in klingon_motifs.items():
            shared = contig_motifs & ref_motifs
            row[f"n_shared_{name}"] = len(shared)
            row[f"shared_motifs_{name}"] = ";".join(sorted(shared)) if shared else ""
            max_shared = max(max_shared, len(shared))

        if max_shared > MIN_SHARED_MOTIFS:
            results.append(row)

print(f"\nContigs with >{MIN_SHARED_MOTIFS} shared motifs with Klingon: {len(results)}\n")

if results:
    sort_col = f"n_shared_{list(klingon_motifs.keys())[0]}"
    out_df = pd.DataFrame(results).sort_values(sort_col, ascending=False)
    out_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(out_df.to_string(index=False))
    print(f"\nResults saved to {OUTPUT_FILE}")
else:
    print("No contigs found matching the criterion.")
