import os
import glob
import subprocess
import textwrap
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
runs_dir = os.path.join(SCRIPT_DIR, "runs_strain")
os.makedirs(runs_dir, exist_ok=True)

isolation_sample_summary = "/home/shuaiw/borg/paper/isolation//GTDB_tree/anno//isolation_sample_summary.tsv"
bam_dir = "/home/shuaiw/borg/paper/isolation//batch2_ccs_bam/"
out_root = "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/"
MGE_csv = "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv"
strain_info_csv = "/home/shuaiw/borg/paper/specificity/iso_99_out/data_tables/Cdb.csv"



df = pd.read_csv(isolation_sample_summary, sep="\t")
print (df)
# Keep only pure samples that have MGE
filtered = df[(df["Pure_anno"] == "pure") & (df["MGE_bool"] == 1) & (df["Motif_Num"] > 1)]
print (f"Filtered to {len(filtered)} pure samples with MGE")
# Extract species name from lineage (e.g. "d__Bacteria;...;s__Klebsiella grimontii" -> "Klebsiella grimontii")
filtered = filtered.copy()
filtered["Species"] = filtered["Lineage"].str.split(";").apply(
    lambda parts: next((p[3:] for p in parts if p.startswith("s__")), None)
)

# Load MGE table once
mge_df = pd.read_csv(MGE_csv)
mge_samples = set(mge_df["prefix"])

# Keep only samples that also appear in MGE_csv
filtered = filtered[filtered["Sample"].isin(mge_samples)]

# Keep only samples that have a ccs.bam file
def has_ccs_bam(sample):
    return bool(glob.glob(os.path.join(bam_dir, f"{sample}*ccs.bam")))

filtered = filtered[filtered["Sample"].apply(has_ccs_bam)]

# Load strain info and join to get secondary_cluster (strain name)
# Cdb.csv genome names encode sample ID as prefix: e.g. "ERR9843365_1_C.fa" -> "ERR9843365"
strain_df = pd.read_csv(strain_info_csv)
strain_df["Sample"] = strain_df["genome"].str.extract(r"^(ERR\d+)")
filtered = filtered.merge(strain_df[["Sample", "secondary_cluster"]], on="Sample", how="inner")
# For each species, pick 2 different strains (different secondary_cluster),
# one representative sample per strain
def pick_two_strains(group):
    strains = group.groupby("secondary_cluster").first().reset_index()
    if len(strains) < 2:
        return None
    return strains.sample(n=2, random_state=42)

strain_pairs = (
    filtered.groupby("Species")
    .apply(pick_two_strains)
    .dropna()
    .reset_index(drop=True)
)

# Species that have at least 2 strains
eligible_species = strain_pairs["Species"].unique()

strain_counts = filtered.groupby("Species")["secondary_cluster"].nunique().sort_values(ascending=False)
print("Strains per species:")
print(strain_counts.to_string())

print(f"Total samples in dataset: {len(df)}")
print(f"Samples after filtering (pure + MGE + ccs.bam + strain info): {len(filtered)}")
print(f"Total eligible species (>=2 strains): {len(eligible_species)}")

# One community using all eligible species (>=2 strains), 2 strains each
n_species = len(eligible_species)
size = n_species * 2
print(f"Building one community: {n_species} species x 2 strains = {size} genomes")
communities = {size: strain_pairs.copy()}

# Generate merged reference and BAM for each community
for size, community in communities.items():
    out_dir = os.path.join(out_root, f"mix_{size}")
    os.makedirs(out_dir, exist_ok=True)

    merged_ref = os.path.join(out_dir, f"mix_{size}.ref.fa")
    merged_bam = os.path.join(out_dir, f"mix_{size}.bam")

    # Write selected strains to CSV
    strain_out = os.path.join(out_dir, f"mix_{size}.strains.csv")
    community[["Species", "secondary_cluster", "genome", "Sample"]].to_csv(strain_out, index=False)
    print(f"Strain list ({len(community)} entries) -> {strain_out}")

    ref_files = []
    bam_files = []

    for _, row in community.iterrows():
        sample = row["Sample"]
        genome = row["genome"]

        # Genome reference
        if not os.path.exists(genome):
            raise FileNotFoundError(f"Reference not found: {genome}")
        ref_files.append(genome)

        # CCS BAM file
        hits = glob.glob(os.path.join(bam_dir, f"{sample}*ccs.bam"))
        if not hits:
            raise FileNotFoundError(f"No ccs.bam found for sample: {sample}")
        bam_files.append(hits[0])

    # Merge references
    print(f"Merging {size} references -> {merged_ref}")
    with open(merged_ref, "wb") as fout:
        for ref in ref_files:
            with open(ref, "rb") as fin:
                fout.write(fin.read())
    subprocess.run(["samtools", "faidx", merged_ref], check=True)

    # Output MGE list for this community
    samples = set(community["Sample"])
    mge_out = os.path.join(out_dir, f"mix_{size}.mge_list.csv")
    mge_subset = mge_df[mge_df["prefix"].isin(samples)][["mge_contig", "mge_type", "mge_length"]].rename(
        columns={"mge_contig": "seq_name", "mge_length": "length"}
    )
    mge_subset.to_csv(mge_out, index=False, sep="\t")
    print(f"MGE list ({len(mge_subset)} entries) -> {mge_out}")

    # Merge BAMs
    print(f"Merging {size} BAMs -> {merged_bam}")
    subprocess.run(
        ["samtools", "merge", "-f", "-@", "20", merged_bam] + bam_files,
        check=True
    )
    subprocess.run(["samtools", "index", "-@", "20", merged_bam], check=True)



    # Generate MODIFI run script
    work_dir = os.path.join(out_dir, "modifi")
    script_path = os.path.join(runs_dir, f"modifi_mix_{size}.sh")
    script_content = textwrap.dedent(f"""\
        #!/bin/bash
        #SBATCH --job-name=mix_{size}
        #SBATCH --partition=standard

        set -euo pipefail

        mkdir -p "{work_dir}"

        python /home/shuaiw/MODIFI/main.py \\
            --work_dir "{work_dir}/mix_{size}" \\
            --unaligned_bam "{merged_bam}" \\
            --whole_ref "{merged_ref}" \\
            --read_type hifi \\
            --mge_file "{mge_out}" \\
            --threads 64
    """)
    with open(script_path, "w") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    print(f"Wrote run script -> {script_path}")

    print(f"Done: mix_{size}")
    # break
