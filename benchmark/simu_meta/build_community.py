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
  # C1 donor-only nested ladder {25,40,58}, one representative strain per species
  # (ladder_10 dropped: 10 genomes are too few for a reliable unmodified IPD baseline)
  python build_community.py ladder --sizes 25,40,58

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


def prep_orphan_bam(sample, genome, ece_contigs, native_dp, target_dp,
                    out_prefix, seed, threads):
    """Build an ORPHAN: the ECE present with real methylation, host genome ABSENT.

    Returns (ece_ref_fa, ece_reads_bam, downsampled). We keep ONLY the ECE contig(s)
    and ONLY the reads that map to them, so the community sees the ECE but not its host:
      1. pbmm2-align the isolate CCS BAM to ITS OWN assembly (kinetics preserved).
      2. collect read NAMES that align over the ECE contig(s).
      3. subset the ORIGINAL unaligned CCS BAM by those names (`samtools view -N`) — this
         keeps the @RG DS:READTYPE=CCS + fi/fp/ri/rp tags intact (NOT addreplacerg/reset,
         which strip them); downsample to target_dp with the same frac as a planted donor.
      4. extract ONLY the ECE contig(s) from the assembly as the orphan reference.
    """
    src = ccs_bam_path(sample)
    aln = f"{out_prefix}.selfaln.bam"
    names = f"{out_prefix}.ece_reads.txt"
    ece_bam0 = f"{out_prefix}.ece_reads.bam"
    ece_fa = f"{out_prefix}.ece.fa"
    pbmm2 = os.path.join(MODIFI_ENV_BIN, "pbmm2")

    # 1. align to own assembly. Run in the (out-of-repo) prep dir so `pbmm2 --sort`'s
    # samtools temp files land there, not in the git working tree.
    workdir = os.path.dirname(os.path.abspath(out_prefix))
    subprocess.run([pbmm2, "align", "--preset", "CCS", "-j", str(threads),
                    genome, src, aln, "--sort"], check=True, cwd=workdir)
    subprocess.run(["samtools", "index", "-@", str(threads), aln], check=True)
    # 2. read names over the ECE contig(s)
    with open(names, "w") as fh:
        view = subprocess.run(["samtools", "view", "-@", str(threads), aln] + list(ece_contigs),
                              check=True, capture_output=True, text=True)
        seen = set()
        for line in view.stdout.splitlines():
            rn = line.split("\t", 1)[0]
            if rn not in seen:
                seen.add(rn)
                fh.write(rn + "\n")
    # 3. subset original unaligned CCS BAM by those names (tags preserved)
    subprocess.run(["samtools", "view", "-@", str(threads), "-N", names,
                    "-b", "-o", ece_bam0, src], check=True)
    # downsample the ECE reads with the same frac a planted donor would use
    frac = 1.0 if (not native_dp or native_dp <= 0) else float(target_dp) / float(native_dp)
    if frac < 1.0:
        ece_bam = f"{out_prefix}.ece_reads.ds.bam"
        s_arg = f"{seed}.{int(round(frac * 1000)):03d}"
        subprocess.run(["samtools", "view", "-@", str(threads), "-b", "-s", s_arg,
                        "-o", ece_bam, ece_bam0], check=True)
        downsampled = True
        os.remove(ece_bam0)
    else:
        ece_bam = ece_bam0
        downsampled = False
    # 4. ECE-only reference (chromosome omitted)
    with open(ece_fa, "w") as fout:
        subprocess.run(["samtools", "faidx", genome] + list(ece_contigs),
                       check=True, stdout=fout)
    for tmp in (aln, aln + ".bai", names):
        try:
            os.remove(tmp)
        except OSError:
            pass
    return ece_fa, ece_bam, downsampled


# ----------------------------------------------------------------- build one community
def build_one(label, donor_rows, background_rows, mge_df, seed, threads,
              keep_prepped=False, orphan_rows=None):
    """donor_rows / background_rows / orphan_rows: DataFrames with Sample, genome,
    Average_DP, Species. orphan_rows contribute ONLY their ECE (host removed) — the
    false-positive negative class (Part B)."""
    out_dir = os.path.join(OUT_ROOT, label)
    prep_dir = os.path.join(out_dir, "prepped_bams")
    os.makedirs(prep_dir, exist_ok=True)

    merged_ref = os.path.join(out_dir, f"{label}.ref.fa")
    merged_bam = os.path.join(out_dir, f"{label}.bam")

    # collect isolate specs (ordered: donors, background, then orphans), subsample in parallel.
    # orphans (host removed, ECE only) are optional and default empty.
    if orphan_rows is None:
        orphan_rows = donor_rows.iloc[0:0]
    ece_by_sample = mge_df.groupby("prefix")["mge_contig"].apply(list)
    specs = ([(row, "donor", DEPTH_DONOR) for _, row in donor_rows.iterrows()] +
             [(row, "background", DEPTH_BG) for _, row in background_rows.iterrows()] +
             [(row, "orphan", DEPTH_DONOR) for _, row in orphan_rows.iterrows()])
    for row, _, _ in specs:
        genome = row["genome"]
        if not (isinstance(genome, str) and os.path.exists(genome)):
            raise FileNotFoundError(f"Reference not found for {row['Sample']}: {genome}")

    def prep_one(spec):
        row, role, target_dp = spec
        sample = row["Sample"]
        if role == "orphan":
            eces = ece_by_sample.get(sample, [])
            ece_fa, merge_bam, downsampled = prep_orphan_bam(
                sample, row["genome"], eces, row.get("Average_DP"), target_dp,
                os.path.join(prep_dir, f"{sample}.orphan"), seed, threads=2)
            genome = ece_fa                       # ECE-only reference (host omitted)
        else:
            out_bam = os.path.join(prep_dir, f"{sample}.{role}.bam")
            merge_bam, downsampled = prep_isolate_bam(
                sample, row.get("Average_DP"), target_dp, out_bam, seed, threads=2)
            genome = row["genome"]
        return {
            "sample": sample, "role": role, "species": row.get("Species"),
            "strain": row.get("Strain", sample), "native_dp": row.get("Average_DP"),
            "target_dp": target_dp, "downsampled": downsampled, "genome": genome,
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

    # ---- MGE list: planted donors + orphans (MODIFI must score orphan ECEs too), deduped
    donor_samples = set(donor_rows["Sample"])
    orphan_samples = set(orphan_rows["Sample"])
    scored_samples = donor_samples | orphan_samples
    mge_subset = (
        mge_df[mge_df["prefix"].isin(scored_samples)][["mge_contig", "mge_type", "mge_length"]]
        .drop_duplicates("mge_contig")
        .rename(columns={"mge_contig": "seq_name", "mge_length": "length"})
    )
    mge_out = os.path.join(out_dir, f"{label}.mge_list.tsv")
    mge_subset.to_csv(mge_out, index=False, sep="\t")
    # ---- orphan ECE list (true host absent) -> evaluator's --orphans negative class
    orphan_eces = sorted(set(
        mge_df[mge_df["prefix"].isin(orphan_samples)]["mge_contig"])) if orphan_samples else []
    if orphan_eces:
        orphan_out = os.path.join(out_dir, f"{label}.orphans.txt")
        with open(orphan_out, "w") as fh:
            fh.write("\n".join(orphan_eces) + "\n")
        print(f"[{label}] MGE list: {len(mge_subset)} ECEs "
              f"({len(mge_subset)-len(orphan_eces)} planted + {len(orphan_eces)} orphan) -> {mge_out}")
    else:
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
    # 1 representative isolate / species, drawn per-seed (not deterministic .first()):
    # across replicates the strain chosen for each species varies, so ladder_58 (all
    # species present) gets genuine compositional variance rather than a fixed set.
    # Reads are additionally re-subsampled per seed inside prep_isolate_bam.
    reps = (donors.sample(frac=1, random_state=seed)         # shuffle rows per-seed
                  .groupby("Species", as_index=False).first()  # 1 random isolate / species
                  .reset_index(drop=True))
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
                    threads, keep_prepped, tag="", orphan_frac=0.0, species=None):
    if tag:
        label = f"{label}_{tag}"
    donors, background, mge_df = load_pool()

    # pick species, then up to K distinct strains per species. --species restricts the
    # donor-species pool to a fixed set (e.g. a single focal species for a strain panel).
    species_pool = donors["Species"].dropna().unique()
    if species:
        species_pool = [s for s in species_pool if s in set(species)]
        if not species_pool:
            raise SystemExit(f"none of --species {species} found among donor species")
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

    # ---- orphan split (Part B): remove the host of a seed-varying orphan_frac of donor
    # SPECIES -> those ECEs become host-absent orphans (false-positive negatives).
    orphan_rows = donor_rows.iloc[0:0]
    orphan_species = set()
    if orphan_frac > 0:
        sp_shuf = pd.Series(sorted(donor_rows["Species"].unique())).sample(
            frac=1, random_state=seed).tolist()
        n_orphan = int(round(len(sp_shuf) * orphan_frac))
        orphan_species = set(sp_shuf[:n_orphan])
        orphan_rows = donor_rows[donor_rows["Species"].isin(orphan_species)].copy()
        donor_rows = donor_rows[~donor_rows["Species"].isin(orphan_species)].copy()

    # background: restrict to ECE-free isolates (MGE_bool==0) so no unscored ECE contig
    # can be spuriously predicted as a host. Select ROUND-ROBIN across species (one strain
    # per species per round, new/non-donor species ordered first) to MAXIMISE species
    # diversity and avoid a few strain-rich species (e.g. H. influenzae, B. pertussis)
    # dominating the background.
    if n_background > 0:
        bg_pool = background[background["MGE_bool"] == 0].copy()
        # orphan species must stay ABSENT -> never let a con-specific host slip in via bg
        if orphan_species:
            bg_pool = bg_pool[~bg_pool["Species"].isin(orphan_species)]
        donor_species = set(donor_rows["Species"].dropna())
        # per-species strain order shuffled by seed (so replicates draw different strains);
        # species ordered new/non-donor first, shuffled within each block by seed.
        groups = {sp: g.sample(frac=1, random_state=seed) for sp, g in bg_pool.groupby("Species")}
        new_sp = pd.Series([s for s in groups if s not in donor_species]).sample(
            frac=1, random_state=seed).tolist()
        ovl_sp = pd.Series([s for s in groups if s in donor_species]).sample(
            frac=1, random_state=seed).tolist()
        sp_order = new_sp + ovl_sp
        selected = []
        r = 0
        while len(selected) < n_background:
            added = 0
            for sp in sp_order:
                if len(selected) >= n_background:
                    break
                g = groups[sp]
                if r < len(g):
                    selected.append(g.iloc[r])
                    added += 1
            if added == 0:
                break  # pool exhausted
            r += 1
        bg_rows = pd.DataFrame(selected)
    else:
        bg_rows = background.iloc[0:0]

    print(f"[{label}] planted-donors={len(donor_rows)} "
          f"({donor_rows['Species'].nunique()} species) "
          f"orphans={len(orphan_rows)} ({len(orphan_species)} species) "
          f"background={len(bg_rows)} "
          f"({bg_rows['Species'].nunique() if len(bg_rows) else 0} bg species) "
          f"=> {len(donor_rows) + len(bg_rows)} host genomes + {len(orphan_rows)} orphan ECEs")
    build_one(label, donor_rows, bg_rows, mge_df, seed, threads, keep_prepped,
              orphan_rows=orphan_rows)


# ---------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    pl = sub.add_parser("ladder", help="donor-only nested species ladder")
    pl.add_argument("--sizes", default="25,40,58",
                    help="comma-separated nested sizes (default 25,40,58; ladder_10 "
                         "dropped — 10 genomes give too few contigs for a reliable "
                         "unmodified IPD baseline; "
                         "donor-species ceiling ~58)")
    pl.add_argument("--only", default=None,
                    help="build only these sizes (subset of --sizes) while keeping the "
                         "full --sizes base fixed for nesting; e.g. --only 10")

    pc = sub.add_parser("community", help="single deep+background community")
    pc.add_argument("--n-species", type=int, required=True)
    pc.add_argument("--strains-per-species", type=int, default=1)
    pc.add_argument("--n-background", type=int, default=0)
    pc.add_argument("--label", required=True)
    pc.add_argument("--orphan-frac", type=float, default=0.0,
                    help="fraction of donor SPECIES whose host is removed to create "
                         "host-absent orphan ECEs (Part B false-positive test); e.g. 0.5")
    pc.add_argument("--species", default=None,
                    help="restrict donor species to this exact name (e.g. 'Escherichia coli') "
                         "for a single-species strain panel; background stays diverse")

    for p in (pl, pc):
        p.add_argument("--seed", type=int, default=SEED)
        p.add_argument("--threads", type=int, default=20)
        p.add_argument("--tag", default="", help="replicate/label suffix, e.g. rep2 -> ladder_25_rep2")
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
                        args.label, args.seed, args.threads, args.keep_prepped, tag=args.tag,
                        orphan_frac=args.orphan_frac,
                        species=[args.species] if args.species else None)


if __name__ == "__main__":
    main()
