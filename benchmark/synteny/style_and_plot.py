#!/usr/bin/env python3
"""Style the LoVis4u annotation tables for the E. faecalis R-M inversion synteny plot.

LoVis4u's first pass (mmseqs clustering) writes feature/locus annotation tables and
colours genes gray. This script edits those tables to:
  - give each gene a short label (REase/MTase/S1/XerC/S2; flanks unlabelled),
  - colour every gene by its homology (mmseqs) group so identical proteins share a
    colour and the specificity-subunit swap (the inversion) is visible as a colour
    rearrangement plus crossing homology links,
  - give each locus a clean panel label (I_1 DOL 14 (G1) etc.).

Then re-runs LoVis4u with -faf/-laf and -sgc-off (keep our colours) + -hl (homology links).

Run inside the lovis4u conda env, from /home/shuaiw/borg/revision/synteny.
"""

import csv
import subprocess
from pathlib import Path

DATA = Path("/home/shuaiw/borg/revision/synteny")
SRC_CSV = DATA / "efaecalis_rm_inversion_genes.csv"
IN_DIR = DATA / "lovis_input"
CLUST_OUT = DATA / "lovis_out"          # first-pass output (tables + mmseqs clusters)
FEAT_IN = CLUST_OUT / "feature_annotation_table.tsv"
LOCUS_IN = CLUST_OUT / "locus_annotation_table.tsv"
FEAT_OUT = DATA / "styled_feature_annotation.tsv"
LOCUS_OUT = DATA / "styled_locus_annotation.tsv"
FINAL_OUT = DATA / "lovis_out_final"

# fixed colours for the conserved core roles; distinct colours cycled for the
# divergent specificity-subunit (S) groups so the inversion reads out by colour.
ROLE_COLOUR = {
    "REase": "#4C72B0",   # blue
    "MTase": "#DD8452",   # orange
    "XerC":  "#55A868",   # green
    "flank": "#D9D9D9",   # light gray
}
S_PALETTE = ["#C44E52", "#8172B3", "#4FB0C6", "#CCB974"]  # red, purple, teal, gold

LABEL = {"REase": "REase", "MTase": "MTase", "S1": "S1", "S2": "S2", "XerC": "XerC",
         "flank": ""}

LOCUS_DESC = {
    "I_1_DOL14_G1": "I_1 DOL 14 (G1)",
    "I_1_DOL35_G2": "I_1 DOL 35 (G2)",
    "I_2": "I_2",
}


def load_roles():
    roles = {}
    with open(SRC_CSV) as f:
        for r in csv.DictReader(f):
            roles[r["gene"]] = r["role"]
    return roles


def main():
    roles = load_roles()

    rows = list(csv.DictReader(open(FEAT_IN), delimiter="\t"))

    # map each mmseqs group -> a colour. Core roles get fixed colours; S groups cycle.
    group_roles = {}
    for r in rows:
        group_roles.setdefault(r["group"], set()).add(roles.get(r["feature_id"], "flank"))
    group_colour = {}
    s_i = 0
    for g, rset in group_roles.items():
        if "REase" in rset:
            group_colour[g] = ROLE_COLOUR["REase"]
        elif "MTase" in rset:
            group_colour[g] = ROLE_COLOUR["MTase"]
        elif "XerC" in rset:
            group_colour[g] = ROLE_COLOUR["XerC"]
        elif rset <= {"flank"}:
            group_colour[g] = ROLE_COLOUR["flank"]
        else:  # specificity-subunit group
            group_colour[g] = S_PALETTE[s_i % len(S_PALETTE)]
            s_i += 1

    # rewrite feature table
    fieldnames = rows[0].keys()
    for r in rows:
        role = roles.get(r["feature_id"], "flank")
        r["name"] = LABEL.get(role, role)
        r["fill_colour"] = group_colour[r["group"]]
        r["stroke_colour"] = "#000000" if role != "flank" else "#B0B0B0"
        r["show_label"] = "1" if role != "flank" else "0"
    with open(FEAT_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {FEAT_OUT}")

    # rewrite locus table (clean panel labels)
    lrows = list(csv.DictReader(open(LOCUS_IN), delimiter="\t"))
    for r in lrows:
        r["description"] = LOCUS_DESC.get(r["sequence_id"], r["sequence_id"])
    with open(LOCUS_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=lrows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(lrows)
    print(f"Wrote {LOCUS_OUT}")

    # re-run lovis4u with our styled tables
    cmd = [
        "lovis4u", "-gb", str(IN_DIR), "-o", str(FINAL_OUT),
        "-faf", str(FEAT_OUT), "-laf", str(LOCUS_OUT),
        "-sgc-off",           # keep our fill/stroke colours
        "-hl",                # homology-link track
        "-cl-off",            # keep our locus order G1 -> G2 -> I_2
        "-mmsi", "0.9", "-mc", "0.5",
        "-llp", "left", "-lls", "description",   # clean panel labels (I_1 DOL 14 (G1) ...)
        # no -safl: labels driven by per-feature show_label (core=1, flanks=0)
    ]
    print("Running:", " ".join(cmd))
    r = subprocess.run(cmd, cwd=str(DATA), capture_output=True, text=True)
    print(r.stdout[-2000:])
    if r.returncode != 0:
        print("STDERR:", r.stderr[-2000:])
    print(f"Final plot: {FINAL_OUT / 'lovis4u.pdf'}")


if __name__ == "__main__":
    main()
