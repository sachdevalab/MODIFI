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


def _read_metadata_sample_env(meta_file):
    """Read prefix_table.tab: return dict sample -> environment (parts[1] -> parts[3])."""
    sample_env_dict = {}
    if not os.path.exists(meta_file):
        return sample_env_dict
    with open(meta_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            sample_env_dict[parts[1]] = parts[3]
    return sample_env_dict


def _collect_modification_fraction_rows(sample, mge, mge_type, all_dir, mge_len=None):
    """Read motif file for (sample, mge), return list of dicts with modification_fraction etc.
    If mge_len is provided, add it to each row."""
    path = motif_file_path(sample, mge, all_dir, use_methylation4=True)
    if not os.path.exists(path):
        path = motif_file_path(sample, mge, all_dir, use_methylation4=False)
    if not os.path.exists(path):
        return []
    mdf = pd.read_csv(path)
    frac_col = mdf.get("fraction", mdf.get("modification_fraction", None))
    if frac_col is None or len(frac_col) == 0:
        return []
    out = []
    for _, m in mdf.iterrows():
        motif_str = str(m.get("motifString", ""))
        center_pos = m.get("centerPos", "")
        motif_string = f"{motif_str}_{center_pos}" if motif_str else ""
        frac_val = m.get("fraction", m.get("modification_fraction", None))
        frac_num = pd.to_numeric(frac_val, errors="coerce")
        if pd.isna(frac_num):
            continue
        n_motifs = m.get("nGenome", m.get("num_motif_sites", m.get("nTotal", None)))
        n_motifs = pd.to_numeric(n_motifs, errors="coerce")
        num_modified = m.get("nDetected", m.get("num_modified_motifs", None))
        num_modified = pd.to_numeric(num_modified, errors="coerce")
        row = {
            "sample": sample,
            "MGE": mge,
            "MGE_type": mge_type,
            "motif_string": motif_string,
            "n_motifs": n_motifs if pd.notna(n_motifs) else None,
            "num_modified_motifs": num_modified if pd.notna(num_modified) else None,
            "modification_fraction": float(frac_num),
        }
        if mge_len is not None:
            row["mge_len"] = mge_len
        out.append(row)
    return out


def run_linkable_mge_modification_fraction(
    linkage_csv=None,
    all_dir="/home/shuaiw/borg/paper/run2/",
    out_dir=None,
    meta_file=None,
    min_mge_length=5000,
):
    """
    For MGEs in mge_host_gc_cov.csv (linkable MGEs), only consider soil samples and
    MGEs with length > min_mge_length bp. Read motif files from all_dir (run2),
    compute modification fraction per (sample, MGE), write CSV and histograms.
    Also collect modification_fraction for plasmid/virus that cannot be linked to host
    in soils (same length filter); write unlinkable CSV and add two subplots.
    All outputs in out_dir.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "isolation"))
    from sample_object import My_sample

    if linkage_csv is None:
        linkage_csv = os.path.join(
            os.path.dirname(__file__),
            "../../tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv",
        )
    if out_dir is None:
        out_dir = "/home/shuaiw/MODIFI/tmp/figures/borg_fig"
    if meta_file is None:
        meta_file = "/home/shuaiw/MODIFI/assembly_pipe/prefix_table.tab"
    if not os.path.exists(linkage_csv):
        print(f"Linkage CSV not found: {linkage_csv}; skipping.")
        return
    ldf = pd.read_csv(linkage_csv)
    if "sample" not in ldf.columns or "MGE" not in ldf.columns or "MGE_type" not in ldf.columns:
        print("Linkage CSV must have columns sample, MGE, MGE_type; skipping.")
        return
    # Only consider soil samples
    if "environment" in ldf.columns:
        ldf = ldf[ldf["environment"] == "soil"]
    else:
        ldf = ldf[ldf["sample"].astype(str).str.lower().str.startswith("soil")]
    if len(ldf) == 0:
        print("No soil samples in linkage CSV; skipping.")
        return
    # Only consider MGEs with length > min_mge_length bp
    if "mge_len" in ldf.columns:
        ldf = ldf[ldf["mge_len"] > min_mge_length]
        if len(ldf) == 0:
            print(f"No soil linkable MGEs with length > {min_mge_length}; skipping.")
            return
    linkable_cols = ["sample", "MGE", "MGE_type"]
    if "mge_len" in ldf.columns:
        linkable_cols.append("mge_len")
    linkable = ldf[linkable_cols].drop_duplicates(subset=["sample", "MGE", "MGE_type"])
    linkable_set = set((row["sample"], row["MGE"]) for _, row in linkable.iterrows())

    # Linkable: collect modification fraction rows
    rows = []
    for _, row in linkable.iterrows():
        sample, mge, mge_type = row["sample"], row["MGE"], row["MGE_type"]
        mge_len = row.get("mge_len", None) if "mge_len" in row.index else None
        rows.extend(_collect_modification_fraction_rows(sample, mge, mge_type, all_dir, mge_len=mge_len))

    if not rows:
        print("No modification fraction rows for linkable MGEs; skipping CSV/plot.")
        return
    df = pd.DataFrame(rows)
    df = df[df["n_motifs"].notna() & (df["n_motifs"] >= 20)]
    if len(df) == 0:
        print("No rows after filtering n_motifs >= 20; skipping CSV/plot.")
        return
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, "linkable_mge_modification_fraction.csv")
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} linkable MGE modification fractions -> {out_csv}")

    # Unlinkable: soil plasmid/virus not in linkable set, length > min_mge_length
    sample_env_dict = _read_metadata_sample_env(meta_file)
    soil_samples = [s for s, e in sample_env_dict.items() if e == "soil"]
    unlinkable_list = []  # (sample, mge, mge_type, mge_len)
    for prefix in soil_samples:
        sample_dir = os.path.join(all_dir, prefix)
        if not os.path.isdir(sample_dir):
            continue
        sample_obj = My_sample(prefix, all_dir)
        mge_dict = sample_obj.read_MGE()
        if mge_dict is None:
            continue
        sample_obj.read_depth()  # populate length_dict for length filter
        length_dict = getattr(sample_obj, "length_dict", None) or {}
        for contig, mge_type in mge_dict.items():
            if mge_type not in ("plasmid", "virus"):
                continue
            if (prefix, contig) in linkable_set:
                continue
            ctg_len = length_dict.get(contig, 0)
            if ctg_len <= min_mge_length:
                continue
            unlinkable_list.append((prefix, contig, mge_type, ctg_len))

    unlinkable_rows = []
    for sample, mge, mge_type, mge_len in unlinkable_list:
        unlinkable_rows.extend(_collect_modification_fraction_rows(sample, mge, mge_type, all_dir, mge_len=mge_len))
    if unlinkable_rows:
        df_unlink = pd.DataFrame(unlinkable_rows)
        df_unlink = df_unlink[df_unlink["n_motifs"].notna() & (df_unlink["n_motifs"] >= 20)]
    else:
        df_unlink = pd.DataFrame(columns=["sample", "MGE", "MGE_type", "motif_string", "n_motifs", "num_modified_motifs", "modification_fraction", "mge_len"])
    out_csv_unlink = os.path.join(out_dir, "unlinkable_mge_modification_fraction.csv")
    df_unlink.to_csv(out_csv_unlink, index=False)
    print(f"Wrote {len(df_unlink)} unlinkable MGE modification fractions -> {out_csv_unlink}")

    # 4-panel figure: linkable plasmid, linkable virus, unlinkable plasmid, unlinkable virus
    fig, axes = plt.subplots(2, 2, figsize=(5 * 2, 4 * 2))
    axes = axes.ravel()
    panels = [
        ("plasmid (linkable)", df, "plasmid"),
        ("virus (linkable)", df, "virus"),
        ("plasmid (unlinkable)", df_unlink, "plasmid"),
        ("virus (unlinkable)", df_unlink, "virus"),
    ]
    for i, (title, data, mge_type) in enumerate(panels):
        ax = axes[i]
        sub = data[data["MGE_type"] == mge_type] if "MGE_type" in data.columns and len(data) > 0 else pd.DataFrame()
        if len(sub) == 0:
            ax.set_visible(False)
            continue
        sns.histplot(sub["modification_fraction"], bins=30, kde=True, ax=ax)
        ax.set_xlim(0, 1)
        ax.set_xlabel("Modification fraction")
        ax.set_ylabel("Count")
        ax.set_title(title)
    plt.tight_layout()
    out_pdf = os.path.join(out_dir, "linkable_mge_modification_fraction_distribution.pdf")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Saved modification fraction distributions -> {out_pdf}")


if __name__ == "__main__":
    # main()
    # plot_profile_motif_distributions()
    # run_gg_run2_genomad_motif_plot()
    # Linkable MGE modification fraction (for count_linkages.R 3-panel figure):
    run_linkable_mge_modification_fraction()
