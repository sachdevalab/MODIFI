import os
import glob
import subprocess
import textwrap
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
runs_dir = os.path.join(SCRIPT_DIR, "runs")
os.makedirs(runs_dir, exist_ok=True)

isolation_sample_summary = "/home/shuaiw/borg/paper/isolation//GTDB_tree/anno//isolation_sample_summary.tsv"
bam_dir = "/home/shuaiw/borg/paper/isolation//batch2_ccs_bam/"
out_root = "/home/shuaiw/borg/paper/linkage/mixed_isolates/"
MGE_csv = "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv"

df = pd.read_csv(isolation_sample_summary, sep="\t")
print (df)
# Keep only pure samples that have MGE
filtered = df[(df["Pure_anno"] == "pure") & (df["MGE_bool"] == 1)]
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

# One representative per species
representatives = filtered.groupby("Species").first().reset_index()

print(f"Total samples in dataset: {len(df)}")
print(f"Samples after filtering (pure + MGE + ccs.bam): {len(filtered)}")
print(f"Total eligible species: {len(representatives)}")
assert len(representatives) >= 50, "Not enough species for 50-species community"

# Build nested subsets: 10 ⊆ 20 ⊆ 30 ⊆ 40 ⊆ 50
sizes = [10, 20, 30, 40, 50]
base = representatives.sample(n=50, random_state=42).reset_index(drop=True)
communities = {size: base.iloc[:size].copy() for size in sizes}

# Verify nesting
for i in range(len(sizes) - 1):
    s1, s2 = sizes[i], sizes[i + 1]
    assert set(communities[s1]["Species"]).issubset(set(communities[s2]["Species"]))

# Generate merged reference and BAM for each community
for size, community in communities.items():
    out_dir = os.path.join(out_root, f"mix_{size}")
    os.makedirs(out_dir, exist_ok=True)

    merged_ref = os.path.join(out_dir, f"mix_{size}.ref.fa")
    merged_bam = os.path.join(out_dir, f"mix_{size}.bam")

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
