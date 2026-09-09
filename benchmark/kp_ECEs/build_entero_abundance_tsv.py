#!/usr/bin/env python3
"""Raw source-data TSV behind fig_entero_abundance.pdf (panel b: Enterobacteriaceae abundance per
infant gut). Reproduces the figure's processing from entero_abundance_sourcedata.csv:
clean GTDB suffixes -> group per (sample, species, genus) -> species_lab (kept if it reaches >=3%
within Enterobacteriaceae in any gut, else 'Other Enterobacteriaceae') -> infant samples only.
Each row is one independent data point (one species in one infant gut) with its raw values."""
import re, csv
from collections import defaultdict

FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
SRC = f"{FIG}/entero_abundance_sourcedata.csv"
OUT = f"{FIG}/fig_entero_abundance_sourcedata.tsv"


def clean(s):
    return re.sub(r"_[A-Z]+ ", " ", re.sub(r"_[A-Z]+$", "", s))


rows = list(csv.DictReader(open(SRC)))
# group by (sample, disp species, genus), sum numeric raw columns
num = ["rel_abundance_entero", "rel_abundance_among_genomes", "genome_mean_depth",
       "mass_fraction_wholemeta", "n_contigs", "total_length_bp"]
agg = defaultdict(lambda: {k: 0.0 for k in num} | {"n_genomes": 0})
for r in rows:
    key = (r["sample"], clean(r["species"]), clean(r["genus"]))
    for k in num:
        agg[key][k] += float(r[k])
    agg[key]["n_genomes"] += 1

# species_lab: kept if max within-Entero fraction across ALL samples >= 3% (matches the R figure)
maxrel = defaultdict(float)
for (s, sp, g), v in agg.items():
    maxrel[sp] = max(maxrel[sp], v["rel_abundance_entero"])
keep = {sp for sp, m in maxrel.items() if m >= 0.03}

out = []
for (s, sp, g), v in agg.items():
    if not s.startswith("infant_"):        # figure shows infant guts only
        continue
    out.append(dict(sample=s, genus=g, species=sp,
                    species_lab=sp if sp in keep else "Other Enterobacteriaceae",
                    rel_abundance_among_genomes=round(v["rel_abundance_among_genomes"], 6),
                    rel_abundance_within_entero=round(v["rel_abundance_entero"], 6),
                    genome_mean_depth=round(v["genome_mean_depth"], 4),
                    mass_fraction_wholemeta=round(v["mass_fraction_wholemeta"], 6),
                    n_genomes=int(v["n_genomes"]), total_length_bp=int(v["total_length_bp"])))

out.sort(key=lambda r: (int(r["sample"].split("_")[1]), -r["rel_abundance_among_genomes"]))
cols = ["sample", "genus", "species", "species_lab", "rel_abundance_among_genomes",
        "rel_abundance_within_entero", "genome_mean_depth", "mass_fraction_wholemeta",
        "n_genomes", "total_length_bp"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(out)
print(f"wrote {OUT}: {len(out)} rows (infant sample x species points), "
      f"{len({r['sample'] for r in out})} infant guts, "
      f"{len({r['species'] for r in out})} species; kept-in-legend: {len(keep)}")
