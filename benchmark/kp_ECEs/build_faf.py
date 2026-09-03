#!/usr/bin/env python3
"""Build the LoVis4u feature-annotation (-faf) table that marks resistance genes on the
infant_15_35_C cluster synteny plot.

Inputs:
  - feature_annotation_table.tsv  (from a first LoVis4u pass; gives the exact feature_id +
    locus_id + coordinates that LoVis4u assigned to each CDS)
  - synteny_amrfinder.tsv         (AMRFinderPlus nucleotide output; contig coords of AMR/STRESS genes)
Output:
  - synteny_faf.tsv   short-form LoVis4u faf: feature_id  name  fill_colour  show_label
  - sourcedata CSV    per-member resistance inventory

Colour: antibiotic AMR = #B83900 ; metal/biocide (STRESS) = #238B8B.
Labels: every AMR gene; for STRESS, one label per distinct (locus, gene symbol) to avoid clutter
(all STRESS genes are still coloured).
"""
import sys, csv
from collections import defaultdict

AMR_COL = "#B83900"      # antibiotic resistance (brick red)
STRESS_COL = "#542788"   # metal/biocide resistance (deep violet, distinct from pastel synteny palette)

feat_tbl, amr_tbl, out_faf, out_src = sys.argv[1:5]


def parse_coords(s):
    # LoVis4u coordinates look like "start:end:strand" (1-based)
    parts = str(s).split(":")
    a, b = int(parts[0]), int(parts[1])
    return min(a, b), max(a, b)


# --- features per locus ---
feats = defaultdict(list)   # locus_id -> [(fid, start, end)]
with open(feat_tbl) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        fid = row.get("feature_id") or row.get("")
        loc = row["locus_id"]
        try:
            s, e = parse_coords(row["coordinates"])
        except Exception:
            continue
        feats[loc].append((fid, s, e))


def best_feature(locus, s, e):
    best, bov = None, 0
    for fid, fs, fe in feats.get(locus, []):
        ov = min(e, fe) - max(s, fs)
        if ov > bov:
            bov, best = ov, fid
    return best


# --- resistance hits ---
rows = []          # (feature_id, name, colour, etype, locus)
seen_feat = {}     # feature_id -> chosen (prefer AMR over STRESS)
per_member = defaultdict(lambda: {"AMR": [], "STRESS": []})
unmapped = []
with open(amr_tbl) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        et = row["Element type"]
        if et not in ("AMR", "STRESS"):
            continue
        loc = row["Contig id"]
        s, e = int(row["Start"]), int(row["Stop"])
        s, e = min(s, e), max(s, e)
        sym = row["Gene symbol"]
        fid = best_feature(loc, s, e)
        per_member[loc]["AMR" if et == "AMR" else "STRESS"].append(sym)
        if fid is None:
            unmapped.append((loc, s, e, sym, et))
            continue
        rows.append((fid, sym, et, loc))

# resolve one row per feature_id (prefer AMR)
by_fid = {}
for fid, sym, et, loc in rows:
    if fid not in by_fid or (et == "AMR" and by_fid[fid][1] != "AMR"):
        by_fid[fid] = (sym, et, loc)

# decide labels: label every AMR gene; for STRESS label only ONE occurrence of each distinct gene
# symbol (globally) to avoid clutter -- all STRESS genes are still coloured.
stress_labeled = set()
faf = []
for fid, (sym, et, loc) in by_fid.items():
    if et == "AMR":
        colour, show = AMR_COL, 1
    else:
        colour = STRESS_COL
        show = 1 if sym not in stress_labeled else 0
        stress_labeled.add(sym)
    faf.append((fid, sym, colour, show))

with open(out_faf, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["feature_id", "name", "fill_colour", "show_label", "group_type"])
    for fid, sym, colour, show in faf:
        w.writerow([fid, sym, colour, show, "labeled"])

# source data: per-member inventory
with open(out_src, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["member", "n_AMR", "n_STRESS", "AMR_genes", "STRESS_genes"])
    for m in sorted(per_member):
        a = per_member[m]["AMR"]; s = per_member[m]["STRESS"]
        w.writerow([m, len(a), len(s), ";".join(sorted(set(a))), ";".join(sorted(set(s)))])

n_amr = sum(1 for _, _, c, _ in faf if c == AMR_COL)
n_str = sum(1 for _, _, c, _ in faf if c == STRESS_COL)
n_lab = sum(sh for *_, sh in faf)
print(f"faf rows: {len(faf)}  (AMR CDS {n_amr}, STRESS CDS {n_str}); labelled {n_lab}")
print(f"unmapped resistance hits (no overlapping CDS): {len(unmapped)}")
for u in unmapped[:10]:
    print("  UNMAPPED", *u)
