"""
Read profile_contig_table.csv and each contig's motif file ({all_dir}/{sample_name}/{sample_name}_methylation4/motifs/{contig}.motifs.csv).
Extract motif_string, modification_fraction, num_motif_sites, num_modified_motifs and write a summary table.
"""
import os
import pandas as pd

# Hardcoded paths (same as profile_borg.py)
all_dir = "/home/shuaiw/borg/paper/gg_run3/"
seq_dir = "/home/shuaiw/borg/paper/borg_data/profile5/"
cluster = "profile"
contig_table_path = os.path.join(seq_dir, f"{cluster}_contig_table.csv")
out_path = os.path.join(seq_dir, f"{cluster}_motif_summary.csv")


def count_motif_composition(motif_seq):
    """GC content of motif sequence (only A,T,C,G). Same logic as motif_density_vs_gc.py."""
    motif_seq = str(motif_seq).upper()
    g_count = motif_seq.count('G')
    c_count = motif_seq.count('C')
    a_count = motif_seq.count('A')
    t_count = motif_seq.count('T')
    length = g_count + c_count + a_count + t_count
    gc_percent = ((g_count + c_count) / length * 100) if length > 0 else 0
    return gc_percent

def genomad_paths_gg_run2(all_dir_gg2, prefix):
    """
    Resolve GenoMad output paths for gg_run2. Supports two layouts:
    (1) .../Genomad/{stem}.contigs_summary/{stem}.contigs_plasmid_summary.tsv  (e.g. soil_100)
    (2) .../Genomad/{stem}_summary/{stem}_plasmid_summary.tsv                  (e.g. soil_1)
    Returns (plasmid_tsv_path, virus_tsv_path); either can be None if not found.
    """
    genomad_dir = os.path.join(all_dir_gg2, prefix, "Genomad")
    if not os.path.isdir(genomad_dir):
        return None, None
    stem = None
    suffix = None  # ".contigs_summary" or "_summary"
    for name in sorted(os.listdir(genomad_dir)):
        d = os.path.join(genomad_dir, name)
        if not os.path.isdir(d):
            continue
        if name.endswith(".contigs_summary"):
            stem = name[:-len(".contigs_summary")]
            suffix = ".contigs_summary"
            break
        if name.endswith("_summary"):
            stem = name[:-len("_summary")]
            suffix = "_summary"
            break
    if stem is None:
        return None, None
    summary_dir = os.path.join(genomad_dir, stem + suffix)
    if suffix == ".contigs_summary":
        plasmid_path = os.path.join(summary_dir, stem + ".contigs_plasmid_summary.tsv")
        virus_path = os.path.join(summary_dir, stem + ".contigs_virus_summary.tsv")
    else:
        plasmid_path = os.path.join(summary_dir, stem + "_plasmid_summary.tsv")
        virus_path = os.path.join(summary_dir, stem + "_virus_summary.tsv")
    return plasmid_path, virus_path


# Motif file path: same as My_contig in sample_object (profile_borg.py uses this via given_species_drep_fast)
def motif_file_path(sample_name, contig, base_dir, use_methylation4=True):
    sub = f"{sample_name}_methylation4" if use_methylation4 else f"{sample_name}_methylation3"
    return os.path.join(base_dir, sample_name, sub, "motifs", f"{contig}.motifs.csv")

def main():
    contig_table = pd.read_csv(contig_table_path)
    contig_table["contig"] = contig_table["contig"].astype(str)

    rows = []
    for _, row in contig_table.iterrows():
        contig = row["contig"]
        sample_name = row["sample_name"]
        borg_ref = row.get("borg_ref", "")
        borg_type = row.get("borg_type", "")
        borg_indicator = row.get("borg_indicator", "")

        path = motif_file_path(sample_name, contig, all_dir, use_methylation4=True)
        if not os.path.exists(path):
            path = motif_file_path(sample_name, contig, all_dir, use_methylation4=False)
        if not os.path.exists(path):
            continue

        motif_df = pd.read_csv(path)
        # Expected columns: motifString, centerPos, fraction, nDetected; nGenome or similar for total sites
        for _, m in motif_df.iterrows():
            motif_str = str(m.get("motifString", ""))
            center_pos = m.get("centerPos", "")
            motif_string = f"{motif_str}_{center_pos}" if motif_str else ""
            motif_gc_percent = count_motif_composition(motif_str) if motif_str else None
            fraction = m.get("fraction", "")
            n_detected = m.get("nDetected", m.get("num_modified_motifs", ""))
            n_genome = m.get("nGenome", m.get("num_motif_sites", m.get("nTotal", "")))

            rows.append({
                "contig": contig,
                "sample_name": sample_name,
                "borg_ref": borg_ref,
                "borg_type": borg_type,
                "borg_indicator": borg_indicator,
                "motif_string": motif_string,
                "motif_gc_percent": motif_gc_percent,
                "modification_fraction": fraction,
                "num_motif_sites": n_genome,
                "num_modified_motifs": n_detected,
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print(f"Loaded motifs for {contig_table['contig'].nunique()} contigs; saved {len(df)} rows -> {out_path}")


def plot_profile_motif_distributions(
    csv_path=None,
    out_dir="/home/shuaiw/MODIFI/tmp/figures/borg_fig",
    min_modified_motifs=20,
    indicator_order=("Mp", "Mini_Chr", "Mp_Virus", "Borg", "mini_Borg"),
):
    """
    Read motif summary CSV, filter by num_modified_motifs, then plot:
    (1) modification fraction distribution, one subplot per borg_indicator;
    (2) motif GC-content distribution, one subplot per borg_indicator.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    if csv_path is None:
        csv_path = out_path
    df = pd.read_csv(csv_path)
    df = df[df["num_modified_motifs"] >= min_modified_motifs]
    present = df["borg_indicator"].dropna().unique()
    present = set(str(p).strip() for p in present if str(p).strip())
    indicators = [x for x in indicator_order if x in present]
    for ind in sorted(present):
        if ind not in indicator_order:
            indicators.append(ind)

    if len(indicators) == 0:
        print("No borg_indicator found in data; skipping plots.")
        return
    n = len(indicators)
    n_cols = min(3, n)
    n_rows = (n + n_cols - 1) // n_cols
    os.makedirs(out_dir, exist_ok=True)

    # Modification fraction distribution
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n == 1:
        axes = np.array([axes])
    axes = axes.ravel()
    for i, ind in enumerate(indicators):
        ax = axes[i]
        sub = df[df["borg_indicator"] == ind]
        sns.histplot(sub["modification_fraction"], bins=50, kde=True, ax=ax)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Modification Fraction")
        ax.set_ylabel("Count")
        ax.set_title(ind)
    for j in range(len(indicators), len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    out_fig = os.path.join(out_dir, "modification_fraction_distribution.pdf")
    plt.savefig(out_fig, bbox_inches="tight")
    plt.close()
    print(f"Saved plot -> {out_fig}")

    # GC-content distribution
    if "motif_gc_percent" not in df.columns:
        def _gc(s):
            if pd.isna(s) or not str(s).strip():
                return None
            parts = str(s).rsplit("_", 1)
            seq = parts[0] if parts else str(s)
            return count_motif_composition(seq)
        df["motif_gc_percent"] = df["motif_string"].apply(_gc)
    df_gc = df[df["motif_gc_percent"].notna()]
    if len(df_gc) > 0:
        fig2, axes2 = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
        if n == 1:
            axes2 = np.array([axes2])
        axes2 = axes2.ravel()
        for i, ind in enumerate(indicators):
            ax = axes2[i]
            sub = df_gc[df_gc["borg_indicator"] == ind]
            sns.histplot(sub["motif_gc_percent"], bins=30, kde=True, ax=ax)
            ax.set_xlim(0, 100)
            ax.set_xlabel("Motif GC content (%)")
            ax.set_ylabel("Count")
            ax.set_title(ind)
        for j in range(len(indicators), len(axes2)):
            axes2[j].set_visible(False)
        plt.tight_layout()
        out_fig_gc = os.path.join(out_dir, "motif_gc_distribution.pdf")
        plt.savefig(out_fig_gc, bbox_inches="tight")
        plt.close()
        print(f"Saved GC distribution plot -> {out_fig_gc}")

def run_gg_run2_genomad_motif_plot(
    all_dir_gg2="/home/shuaiw/methylation/data/borg/paper/gg_run2",
    min_depth=10,
    min_length=5000,
    out_fig_dir="/home/shuaiw/MODIFI/tmp/figures/borg_fig",
    out_csv=None,
):
    """
    Read GenoMad plasmid/virus results and depth from gg_run2 using sample_object.My_sample;
    keep contigs with depth > min_depth and length >= min_length, collect their motifs,
    store the motif DataFrame to CSV, and plot modification fraction distribution.
    """
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "isolation"))
    from sample_object import My_sample
    from find_borg import find_assembly

    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    # Same reference_fasta / fai as profile_borg: from find_assembly()
    fasta_dict = find_assembly()
    rows = []
    for prefix in sorted(os.listdir(all_dir_gg2)):
        sample_dir = os.path.join(all_dir_gg2, prefix)
        if not os.path.isdir(sample_dir):
            continue
        sample = My_sample(prefix, all_dir_gg2)
        reference = fasta_dict.get(prefix) if fasta_dict else None
        sample.reference_fasta = reference
        sample.fai = (sample.reference_fasta + ".fai") if sample.reference_fasta else None
        # GenoMad paths for gg_run2: .../Genomad/{stem}.contigs_summary/{stem}.contigs_plasmid_summary.tsv
        plasmid_tsv, virus_tsv = genomad_paths_gg_run2(all_dir_gg2, prefix)
        if plasmid_tsv is not None:
            sample.genomad_plasmid = plasmid_tsv
        if virus_tsv is not None:
            sample.genomad_virus = virus_tsv
        sample.collect_all_mges()
        ## read the depth file
        sample.read_depth()
        if not sample.depth_dict:
            continue
        if not sample.mge_genomad_dict:
            continue
        high_dp_mges = [
            ctg for ctg in sample.mge_genomad_dict
            if sample.depth_dict.get(ctg, 0) > min_depth
            and int(sample.mge_genomad_dict[ctg].length) >= min_length
        ]
        for contig in high_dp_mges:
            mge_type = sample.mge_genomad_dict[contig].type  # plasmid or virus
            path = motif_file_path(prefix, contig, all_dir_gg2, use_methylation4=True)
            if not os.path.exists(path):
                path = motif_file_path(prefix, contig, all_dir_gg2, use_methylation4=False)
            if not os.path.exists(path):
                continue
            motif_df = pd.read_csv(path)
            for _, m in motif_df.iterrows():
                motif_str = str(m.get("motifString", ""))
                center_pos = m.get("centerPos", "")
                motif_string = f"{motif_str}_{center_pos}" if motif_str else ""
                fraction = m.get("fraction", "")
                n_detected = m.get("nDetected", m.get("num_modified_motifs", ""))
                n_sites = m.get("nGenome", m.get("num_motif_sites", m.get("nTotal", "")))
                rows.append({
                    "contig": contig,
                    "sample_name": prefix,
                    "mge_type": mge_type,
                    "motif_string": motif_string,
                    "modification_fraction": fraction,
                    "num_motif_sites": n_sites,
                    "num_modified_motifs": n_detected,
                })
    if not rows:
        print("No motif rows from GenoMad plasmid/virus contigs; skipping plot.")
        return
    df = pd.DataFrame(rows)
    df = df[df["num_modified_motifs"] >= 20]
    if len(df) == 0:
        print("No rows after num_modified_motifs >= 20; skipping plot.")
        return
    os.makedirs(out_fig_dir, exist_ok=True)
    if out_csv is None:
        out_csv = os.path.join(out_fig_dir, "gg_run2_genomad_motif_summary.csv")
    df.to_csv(out_csv, index=False)
    print(f"Saved motif DataFrame -> {out_csv} ({len(df)} rows)")
    types = ["plasmid", "virus"]
    present = [t for t in types if (df["mge_type"] == t).any()]
    n = len(present)
    if n == 0:
        print("No plasmid/virus types in df; skipping plot.")
        return
    n_cols = min(2, n)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n == 1:
        axes = np.array([axes])
    axes = axes.ravel()
    for i, mge_type in enumerate(present):
        ax = axes[i]
        sub = df[df["mge_type"] == mge_type]
        sns.histplot(sub["modification_fraction"], bins=50, kde=True, ax=ax)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Modification Fraction")
        ax.set_ylabel("Count")
        ax.set_title(mge_type)
    for j in range(len(present), len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    out_fig = os.path.join(out_fig_dir, "gg_run2_genomad_modification_fraction_distribution.pdf")
    plt.savefig(out_fig, bbox_inches="tight")
    plt.close()
    print(f"Saved GenoMad modification fraction plot -> {out_fig}")


if __name__ == "__main__":
    # main()
    # plot_profile_motif_distributions()
    run_gg_run2_genomad_motif_plot()



# Uncomment to run GenoMad plasmid/virus motif plot for gg_run2:
# run_gg_run2_genomad_motif_plot()
