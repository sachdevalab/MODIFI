#!/usr/bin/env python3
"""
Step 4 - Enterobacteriaceae relative-abundance table for the infant guts (reviewer point 5.1).

For each infant gut sample, join per-contig coverage (abundance_<sample>.tsv) to per-contig GTDB-Tk
taxonomy (run2/<sample>/GTDB/gtdbtk.bac120.summary.tsv) and per-contig length. Restrict to
Enterobacteriaceae contigs, aggregate to species as length-weighted mean coverage (genome-size
neutral), and express each species as a fraction of the sample's Enterobacteriaceae community.

Also records the Enterobacteriaceae fraction of the whole classified community per sample (context:
are these guts Enterobacteriaceae-dominated?), and flags whether Klebsiella is present / dominant.

Outputs a tidy long-form source-data CSV consumed by plot_abundance_entero.R.
"""
import os
import re
import glob
import pandas as pd
import numpy as np

R2 = "/groups/banfield/projects/multienv/methylation/data/borg/paper/run2"
AB = "/groups/banfield/projects/multienv/methylation/binning/test"
LEN_TSV = "/home/shuaiw/borg/revision/kp_eces/infant_contig_lengths.tsv"
OUT_DIR = "/home/shuaiw/borg/revision/kp_eces"
FIG_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"


def parse_lineage(classification):
    """Return dict of rank->name from a GTDB lineage string."""
    out = {}
    for tok in str(classification).split(";"):
        tok = tok.strip()
        m = re.match(r"^([dpcofgs])__(.*)$", tok)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def main():
    lengths = pd.read_csv(LEN_TSV, sep="\t", header=None, names=["contig", "length"])
    len_map = dict(zip(lengths["contig"], lengths["length"]))

    samples = sorted(
        [os.path.basename(os.path.dirname(os.path.dirname(p)))
         for p in glob.glob(f"{R2}/infant_*/GTDB/gtdbtk.bac120.summary.tsv")],
        key=lambda s: int(s.split("_")[1]))

    rows = []            # per sample x species Enterobacteriaceae abundance
    ctx_rows = []        # per sample context (Enterobacteriaceae community fraction)
    for s in samples:
        gtdb = f"{R2}/{s}/GTDB/gtdbtk.bac120.summary.tsv"
        abund = f"{AB}/{s}/abundance_{s}.tsv"
        if not (os.path.exists(gtdb) and os.path.exists(abund)):
            continue
        tax = pd.read_csv(gtdb, sep="\t")[["user_genome", "classification"]]
        tax = tax.rename(columns={"user_genome": "contig"})
        cov = pd.read_csv(abund, sep="\t")
        cov.columns = ["contig", "coverage"]
        cov["coverage"] = pd.to_numeric(cov["coverage"], errors="coerce").fillna(0.0)
        cov["length"] = cov["contig"].map(len_map)
        cov = cov.dropna(subset=["length"])
        cov["length"] = cov["length"].astype(float)
        # true whole-metagenome denominator = ALL assembled contigs (classified or not),
        # so the Enterobacteriaceae fraction is not inflated by dropping unclassified biomass.
        total_meta_mass = float((cov["coverage"] * cov["length"]).sum())
        df = tax.merge(cov, on="contig", how="inner")
        df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce").fillna(0.0)
        # taxonomy ranks
        lin = df["classification"].map(parse_lineage)
        df["family"] = lin.map(lambda d: d.get("f", ""))
        df["genus"] = lin.map(lambda d: d.get("g", ""))
        df["species"] = lin.map(lambda d: d.get("s", ""))
        # weighted "mass" = coverage * length (aligned bases ~ organism abundance)
        df["mass"] = df["coverage"] * df["length"]

        classified = df[df["species"] != ""]

        # --- cell-abundance metric: length-weighted mean coverage (depth) per genome ---
        # Depth is genome-size neutral (copies of a genome), unlike coverage*length (aligned bases,
        # which favours larger genomes). Each GTDB-classified species collapses its contigs to ONE
        # length-weighted mean depth (a reconstructed genome). Unclassified biomass is mostly short
        # contig fragments and is NOT a set of genomes, so it is NOT used as a denominator (counting
        # each fragment as a genome would massively inflate community depth). We therefore express
        # relative abundance among reconstructed genomes; total assembled base fraction of the whole
        # metagenome is also recorded (mass_fraction_*) for the biomass context.
        sp_depth = {}       # species -> length-weighted mean depth (one reconstructed genome)
        sp_genus, sp_fam = {}, {}
        for sp_name, g in classified.groupby("species"):
            tl = g["length"].sum()
            sp_depth[sp_name] = (g["mass"].sum() / tl) if tl else 0.0
            sp_genus[sp_name] = g["genus"].iloc[0]
            sp_fam[sp_name] = g["family"].iloc[0]
        total_genome_depth = sum(sp_depth.values())                 # over reconstructed genomes
        ent_species = [x for x in sp_depth if sp_fam[x] == "Enterobacteriaceae"]
        ent_depth = sum(sp_depth[x] for x in ent_species)

        # whole-metagenome base fractions (aligned bases; handles fragments correctly)
        ent_mass = classified[classified["family"] == "Enterobacteriaceae"]["mass"].sum()

        ctx_rows.append({
            "sample": s,
            # cell-abundance (depth) share of Enterobacteriaceae among reconstructed genomes
            "entero_depthfrac_of_genomes": round(ent_depth / total_genome_depth, 4) if total_genome_depth else 0.0,
            # whole-metagenome base fraction (standard relative abundance, includes unclassified in denom)
            "entero_massfrac_of_metagenome": round(ent_mass / total_meta_mass, 4) if total_meta_mass else 0.0,
            "n_entero_species": len(ent_species),
            "klebsiella_present": any("Klebsiella" in str(sp_genus[x]) for x in ent_species),
        })
        if ent_depth <= 0:
            continue
        # per Enterobacteriaceae species: relative abundance BY DEPTH (cell abundance).
        #  - rel_abundance_entero: fraction of the sample's Enterobacteriaceae depth
        #  - rel_abundance_among_genomes: fraction of ALL reconstructed genomes' depth (genome-size
        #    neutral; low value = the taxon is a minor member of the resolved community)
        ent = classified[classified["family"] == "Enterobacteriaceae"]
        for sp, g in ent.groupby("species"):
            depth = sp_depth[sp]
            sp_mass = g["mass"].sum()
            rows.append({
                "sample": s,
                "species": re.sub(r"^s__", "", sp),
                "genus": g["genus"].iloc[0],
                "rel_abundance_entero": round(depth / ent_depth, 5),
                "rel_abundance_among_genomes": round(depth / total_genome_depth, 5) if total_genome_depth else 0.0,
                "genome_mean_depth": round(depth, 3),
                "mass_fraction_wholemeta": round(sp_mass / total_meta_mass, 5) if total_meta_mass else 0.0,
                "n_contigs": len(g),
                "total_length_bp": int(g["length"].sum()),
            })

    ab_df = pd.DataFrame(rows).sort_values(["sample", "rel_abundance_entero"],
                                           ascending=[True, False])
    ctx_df = pd.DataFrame(ctx_rows)
    os.makedirs(FIG_DIR, exist_ok=True)
    src = os.path.join(FIG_DIR, "entero_abundance_sourcedata.csv")
    ab_df.to_csv(src, index=False)
    ab_df.to_csv(os.path.join(OUT_DIR, "entero_abundance.tsv"), sep="\t", index=False)
    ctx_df.to_csv(os.path.join(FIG_DIR, "entero_community_context_sourcedata.csv"), index=False)

    print(f"[abundance] samples with Enterobacteriaceae: {ab_df['sample'].nunique()}")
    print(f"[abundance] wrote {src} ({len(ab_df)} sample x species rows)")
    # quick summary: is Klebsiella the dominant Enterobacteriaceae where present?
    kleb = ab_df[ab_df["genus"] == "Klebsiella"]
    print(f"[abundance] guts with Klebsiella: {kleb['sample'].nunique()}")
    # dominant Enterobacteriaceae species per gut
    dom = ab_df.loc[ab_df.groupby("sample")["rel_abundance_entero"].idxmax()]
    kleb_dom = dom[dom["genus"] == "Klebsiella"]["sample"].nunique()
    print(f"[abundance] guts where Klebsiella is the DOMINANT Enterobacteriaceae: "
          f"{kleb_dom}/{dom['sample'].nunique()}")
    print("\n[abundance] per-gut Enterobacteriaceae context:")
    print(ctx_df.to_string(index=False))


if __name__ == "__main__":
    main()
