## given a gff file, only keep the lines with score >= 30


import os
import pandas as pd
import sys


# raw_gff = "/home/shuaiw/borg/paper/borg_data/filter_motifs/BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.reprocess.gff"
raw_gff = "/home/shuaiw/borg/paper/borg_data/batch_export/soil_80/gffs/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.gff"
# raw_gff = ""
filtered_gff = raw_gff.replace(".gff", "_score30.gff")


def filter_gff(raw_gff, filtered_gff, score_cutoff=30):
    with open(raw_gff, 'r') as infile, open(filtered_gff, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            score = float(parts[5])
            if score >= score_cutoff:
                outfile.write(line)
    print(f"Filtered GFF file saved to {filtered_gff}")


filter_gff(raw_gff, filtered_gff, score_cutoff=30)