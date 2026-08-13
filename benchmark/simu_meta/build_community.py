#!/usr/bin/env python3
"""
build_community.py — generalized in-silico mock-community builder for the MODIFI
reviewer-response benchmark (dataset C1 and beyond).

Generalizes benchmark/linkage/merge_isolate.py:
  - donor selection (pure + MGE-bearing + Motif_Num>1), species grouping, optional
    multi-strain via dRep Cdb.csv secondary_cluster
  - optional background hosts (pure isolates that contribute NO scored ECEs) to build
    community complexity toward the full ~153-species pool
  - per-isolate downsampling to a target depth (donors >=30x, background ~10x) BEFORE
    merge, preserving PacBio kinetics tags (fi/fp/ri/rp)
  - per-isolate read-group tagging (harmless here; the cross-strain mis-mapping ground
    truth in Part E depends on it)
  - concatenate references, merge BAMs, write per-community MGE list, and emit a SLURM
    run script using the CURRENT main.py CLI (-o / -b / -r / --read_type / --mge_file)

Everything a mock community carries is REAL isolate reads (real IPD kinetics); ground
truth = SRA-accession prefix on each contig (e.g. ERR10042285_1_L -> ERR10042285).

Usage examples:
  # C1 donor-only nested ladder {10,25,40,66}, one representative strain per species
  python build_community.py ladder --sizes 10,25,40,66

  # A single deep + background community (~80 genomes)
  python build_community.py community --n-species 40 --strains-per-species 1 \
        --n-background 40 --label mix_bg80
"""

import os
import glob
import argparse
import subprocess
import textwrap
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

# ---------------------------------------------------------------------------- paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODIFI_MAIN = "/home/shuaiw/MODIFI/main.py"
# Use the `modifi` conda env explicitly (bare `python` resolves to a pyenv shim that
# lacks MODIFI's deps). The env also provides pbmm2/pbindex/samtools/pbmotifmaker.
MODIFI_ENV_BIN = "/home/shuaiw/miniconda3/envs/modifi/bin"
MODIFI_PY = os.path.join(MODIFI_ENV_BIN, "python")

ISOLATION_SUMMARY = "/home/shuaiw/borg/paper/isolation/GTDB_tree/anno/isolation_sample_summary.tsv"
BAM_DIR = "/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/"
MGE_CSV = "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv"
CDB_CSV = "/home/shuaiw/borg/paper/specificity/iso_99_out/data_tables/Cdb.csv"

OUT_ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
RUNS_DIR = os.path.join(SCRIPT_DIR, "cmds")   # generated SLURM scripts (git-ignored)

DEPTH_DONOR = 30.0
DEPTH_BG = 10.0
SEED = 42


# ------------------------------------------------------------------- data loading
def species_of(lineage):
    if not isinstance(lineage, str):
        return None
    return next((p[3:] for p in lineage.split(";") if p.startswith("s__")), None)


def has_ccs_bam(sample):
    return bool(glob.glob(os.path.join(BAM_DIR, f"{sample}*ccs.bam")))


def ccs_bam_path(sample):
    hits = glob.glob(os.path.join(BAM_DIR, f"{sample}*ccs.bam"))
    if not hits:
        raise FileNotFoundError(f"No ccs.bam for sample {sample}")
    return hits[0]


def load_pool():
    """Return (donors, background, mge_df).

    Restricted to the manuscript's high-quality set: pure + Average_DP >= 10x
    (= 1,420 isolates / 154 species / 7 phyla). donors are unchanged by this filter
    (561/58 — all already >=10x); it removes low-coverage background candidates.

    donors     : HQ-pure + MGE + Motif_Num>1 + has bam + in MGE table, with Species/Strain
    background : HQ-pure isolates that are NOT donors but have a genome + bam (complexity)
    """
    df = pd.read_csv(ISOLATION_SUMMARY, sep="\t")
    pure = df[(df["Pure_anno"] == "pure") & (df["Average_DP"] >= 10)].copy()  # manuscript HQ set (1,420)
    pure["Species"] = pure["Lineage"].apply(species_of)

    mge_df = pd.read_csv(MGE_CSV)
    mge_samples = set(mge_df["prefix"])

    is_donor = (
        (pure["MGE_bool"] == 1)
        & (pure["Motif_Num"] > 1)
        & pure["Sample"].isin(mge_samples)
        & pure["Sample"].apply(has_ccs_bam)
        & pure["Species"].notna()
    )
    donors = pure[is_donor].copy()

    # strain label from dRep Cdb.csv (genome name prefix -> secondary_cluster)
    cdb = pd.read_csv(CDB_CSV)
    cdb["Sample"] = cdb["genome"].str.extract(r"^([A-Za-z]+\d+)")
    strain_map = cdb.drop_duplicates("Sample").set_index("Sample")["secondary_cluster"]
    donors["Strain"] = donors["Sample"].map(strain_map).fillna(donors["Sample"])

    background = pure[~pure.index.isin(donors.index)].copy()
    background = background[
        background["Sample"].apply(has_ccs_bam) & background["genome"].apply(
            lambda g: isinstance(g, str) and os.path.exists(g)
        )
    ]

    print(f"[pool] pure={len(pure)}  donors={len(donors)} "
          f"({donors['Species'].nunique()} species)  background={len(background)}")
    return donors, background, mge_df


# --------------------------------------------------------------- bam preparation
def prep_isolate_bam(sample, native_dp, target_dp, out_bam, seed, threads):
    """Return (bam_path_to_merge, downsampled) for an isolate CCS BAM at ~target_dp.

    IMPORTANT: we do NOT rewrite the @RG. The original PacBio @RG carries
    `DS:READTYPE=CCS;Ipd:CodecV1=...;PulseWidth:CodecV1=...`, which pbmm2 requires to
    recognise the reads as alignable CCS/HiFi. `samtools addreplacerg` strips it and
    silently zeroes alignment coverage for some isolates — so we preserve the header.
    `samtools view -s` keeps the full header (and all kinetics tags). If native depth is
    below target we use the original BAM as-is (no copy).

    Read origin for later (Part E) is recoverable from the preserved original @RG IDs.
    """
    src = ccs_bam_path(sample)
    frac = 1.0 if (not native_dp or native_dp <= 0) else float(target_dp) / float(native_dp)

    if frac >= 1.0:
        return src, False  # not downsampled: merge the original BAM directly

    # samtools -s takes SEED.FRACTION as one float, e.g. 42.714 = seed 42, keep 71.4%
    s_arg = f"{seed}.{int(round(frac * 1000)):03d}"
    subprocess.run(
        ["samtools", "view", "-@", str(threads), "-b", "-s", s_arg, "-o", out_bam, src],
        check=True,
    )
    return out_bam, True


# ----------------------------------------------------------------- build one community
def build_one(label, donor_rows, background_rows, mge_df, seed, threads,
              keep_prepped=False):
    """donor_rows / background_rows: DataFrames with Sample, genome, Average_DP, Species."""
    out_dir = os.path.join(OUT_ROOT, label)
    prep_dir = os.path.join(out_dir, "prepped_bams")
    os.makedirs(prep_dir, exist_ok=True)

    merged_ref = os.path.join(out_dir, f"{label}.ref.fa")
    merged_bam = os.path.join(out_dir, f"{label}.bam")

    # collect isolate specs (ordered: donors then background), then subsample in parallel
    specs = ([(row, "donor", DEPTH_DONOR) for _, row in donor_rows.iterrows()] +
             [(row, "background", DEPTH_BG) for _, row in background_rows.iterrows()])
    for row, _, _ in specs:
        genome = row["genome"]
        if not (isinstance(genome, str) and os.path.exists(genome)):
            raise FileNotFoundError(f"Reference not found for {row['Sample']}: {genome}")

    def prep_one(spec):
        row, role, target_dp = spec
        sample = row["Sample"]
        out_bam = os.path.join(prep_dir, f"{sample}.{role}.bam")
        merge_bam, downsampled = prep_isolate_bam(
            sample, row.get("Average_DP"), target_dp, out_bam, seed, threads=2)
        return {
            "sample": sample, "role": role, "species": row.get("Species"),
            "strain": row.get("Strain", sample), "native_dp": row.get("Average_DP"),
            "target_dp": target_dp, "downsampled": downsampled, "genome": row["genome"],
            "merge_bam": merge_bam,
        }

    # parallel subsampling (samtools is I/O-bound; run many at once)
    n_workers = min(len(specs), max(1, threads // 2))
    print(f"[{label}] subsampling {len(specs)} isolates with {n_workers} workers ...")
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        results = list(ex.map(prep_one, specs))   # preserves input order
    manifest = results
    ref_files = [r["genome"] for r in results]
    bam_files = [r["merge_bam"] for r in results]

    # ---- merge references (contigs are SRA-prefixed => globally unique)
    print(f"[{label}] concatenating {len(ref_files)} references -> {merged_ref}")
    with open(merged_ref, "wb") as fout:
        for ref in ref_files:
            with open(ref, "rb") as fin:
                fout.write(fin.read())
    subprocess.run(["samtools", "faidx", merged_ref], check=True)

    # ---- MGE list: donors only, deduped on seq_name
    donor_samples = set(donor_rows["Sample"])
    mge_subset = (
        mge_df[mge_df["prefix"].isin(donor_samples)][["mge_contig", "mge_type", "mge_length"]]
        .drop_duplicates("mge_contig")
        .rename(columns={"mge_contig": "seq_name", "mge_length": "length"})
    )
    mge_out = os.path.join(out_dir, f"{label}.mge_list.tsv")
    mge_subset.to_csv(mge_out, index=False, sep="\t")
    print(f"[{label}] MGE list: {len(mge_subset)} ECEs -> {mge_out}")

    # ---- merge BAMs (unaligned CCS; MODIFI aligns internally with pbmm2)
    print(f"[{label}] merging {len(bam_files)} BAMs -> {merged_bam}")
    subprocess.run(
        ["samtools", "merge", "-f", "-c", "-p", "-@", str(threads), merged_bam] + bam_files,
        check=True,
    )
    subprocess.run(["samtools", "index", "-@", str(threads), merged_bam], check=True)

    # ---- manifest (ground-truth membership)
    man_df = pd.DataFrame(manifest)
    man_out = os.path.join(out_dir, f"{label}.manifest.csv")
    man_df.to_csv(man_out, index=False)
    n_down = int(man_df["downsampled"].sum()) if len(man_df) else 0
    print(f"[{label}] manifest -> {man_out}  "
          f"(genomes={len(man_df)}, downsampled={n_down}, "
          f"native-depth-kept={len(man_df) - n_down})")

    if not keep_prepped:
        for b in bam_files:
            # only delete our subsampled copies, NEVER the original isolate BAMs
            if os.path.dirname(os.path.abspath(b)) == os.path.abspath(prep_dir):
                try:
                    os.remove(b)
                except OSError:
                    pass

    emit_run_script(label, merged_bam, merged_ref, mge_out, out_dir)
    return out_dir


# MODIFI run resources (decoupled from build threads). 32 cores schedules on a partial
# node far sooner than a full 64-core request.
RUN_THREADS = 32
RUN_MEM_GB = 120


def emit_run_script(label, merged_bam, merged_ref, mge_out, out_dir,
                    run_threads=RUN_THREADS, mem_gb=RUN_MEM_GB):
    os.makedirs(RUNS_DIR, exist_ok=True)
    work_dir = os.path.join(out_dir, "modifi")
    script_path = os.path.join(RUNS_DIR, f"modifi_{label}.sh")
    content = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=sm_{label}
        #SBATCH --partition=standard
        #SBATCH --cpus-per-task={run_threads}
        #SBATCH --mem={mem_gb}G
        set -euo pipefail
        export PATH="{MODIFI_ENV_BIN}:$PATH"
        mkdir -p "{work_dir}"

        {MODIFI_PY} {MODIFI_MAIN} \\
            -o "{work_dir}/{label}" \\
            -b "{merged_bam}" \\
            -r "{merged_ref}" \\
            --read_type hifi \\
            --mge_file "{mge_out}" \\
            --threads {run_threads}
    """)
    with open(script_path, "w") as f:
        f.write(content)
    os.chmod(script_path, 0o755)
    print(f"[{label}] run script -> {script_path}")


# ----------------------------------------------------------------------- ladder mode
def build_ladder(sizes, seed, threads, keep_prepped, only=None, tag=""):
    suffix = f"_{tag}" if tag else ""
    donors, _, mge_df = load_pool()
    reps = donors.groupby("Species").first().reset_index()  # 1 representative / species
    max_size = max(sizes)
    assert len(reps) >= max_size, (
        f"Only {len(reps)} ECE-donor species available; requested max size {max_size}. "
        f"Reduce --sizes (donor-species ceiling is {len(reps)})."
    )
    # sample the full base ONCE so every --only subset is nested & reproducible
    base = reps.sample(n=max_size, random_state=seed).reset_index(drop=True)

    communities = {s: base.iloc[:s].copy() for s in sizes}
    for i in range(len(sizes) - 1):
        assert set(communities[sizes[i]]["Species"]).issubset(
            set(communities[sizes[i + 1]]["Species"])
        ), "nesting invariant violated"

    to_build = only if only else sizes
    empty_bg = donors.iloc[0:0]
    for s in to_build:
        build_one(f"ladder_{s}{suffix}", communities[s], empty_bg, mge_df, seed, threads,
                  keep_prepped)


# --------------------------------------------------------------- single community mode
def build_community(n_species, strains_per_species, n_background, label, seed,
                    threads, keep_prepped, tag=""):
    if tag:
        label = f"{label}_{tag}"
    donors, background, mge_df = load_pool()

    # pick species, then up to K distinct strains per species
    species_pool = donors["Species"].dropna().unique()
    rng = pd.Series(species_pool).sample(
        n=min(n_species, len(species_pool)), random_state=seed
    ).tolist()

    picked = []
    for sp in rng:
        grp = donors[donors["Species"] == sp]
        strains = grp.groupby("Strain").first().reset_index()
        k = min(strains_per_species, len(strains))
        picked.append(strains.sample(n=k, random_state=seed))
    donor_rows = pd.concat(picked, ignore_index=True)

    # background: restrict to isolates with NO detected ECE (MGE_bool==0) so no unscored
    # ECE contig can be spuriously predicted as a host; prefer species not among the
    # donors to maximise community richness, then fill from the rest.
    if n_background > 0:
        bg_pool = background[background["MGE_bool"] == 0].copy()
        donor_species = set(donor_rows["Species"].dropna())
        new_sp = bg_pool[~bg_pool["Species"].isin(donor_species)]
        rest = bg_pool[bg_pool["Species"].isin(donor_species)]
        # one isolate per new species first (richness), then extras if needed
        new_first = new_sp.groupby("Species", group_keys=False).apply(
            lambda g: g.sample(1, random_state=seed))
        ordered = pd.concat([new_first.sample(frac=1, random_state=seed),
                             new_sp.drop(new_first.index), rest], ignore_index=False)
        bg_rows = ordered.drop_duplicates("Sample").head(n_background)
    else:
        bg_rows = background.iloc[0:0]

    print(f"[{label}] donors={len(donor_rows)} "
          f"({donor_rows['Species'].nunique()} species) background={len(bg_rows)} "
          f"({bg_rows['Species'].nunique() if len(bg_rows) else 0} bg species) "
          f"=> {len(donor_rows) + len(bg_rows)} genomes")
    build_one(label, donor_rows, bg_rows, mge_df, seed, threads, keep_prepped)


# ---------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    pl = sub.add_parser("ladder", help="donor-only nested species ladder")
    pl.add_argument("--sizes", default="10,25,40,58",
                    help="comma-separated nested sizes (default 10,25,40,58; "
                         "donor-species ceiling ~58)")
    pl.add_argument("--only", default=None,
                    help="build only these sizes (subset of --sizes) while keeping the "
                         "full --sizes base fixed for nesting; e.g. --only 10")

    pc = sub.add_parser("community", help="single deep+background community")
    pc.add_argument("--n-species", type=int, required=True)
    pc.add_argument("--strains-per-species", type=int, default=1)
    pc.add_argument("--n-background", type=int, default=0)
    pc.add_argument("--label", required=True)

    for p in (pl, pc):
        p.add_argument("--seed", type=int, default=SEED)
        p.add_argument("--threads", type=int, default=20)
        p.add_argument("--tag", default="", help="replicate/label suffix, e.g. rep2 -> ladder_10_rep2")
        p.add_argument("--keep-prepped", action="store_true",
                       help="keep per-isolate downsampled BAMs (default: delete after merge)")

    args = ap.parse_args()
    os.makedirs(OUT_ROOT, exist_ok=True)

    if args.mode == "ladder":
        sizes = sorted(int(x) for x in args.sizes.split(","))
        only = sorted(int(x) for x in args.only.split(",")) if args.only else None
        build_ladder(sizes, args.seed, args.threads, args.keep_prepped, only=only, tag=args.tag)
    else:
        build_community(args.n_species, args.strains_per_species, args.n_background,
                        args.label, args.seed, args.threads, args.keep_prepped, tag=args.tag)


if __name__ == "__main__":
    main()
