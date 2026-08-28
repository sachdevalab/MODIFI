#!/usr/bin/env python3
"""
Per-sample host-motif SEQUENCE density (motifs/kb): each Methanoperedens strain vs the Borg.

Reviewer question: does the Borg carry the host's modification motifs at all, and at what
density in its sequence (motifs/kb) compared to the host?

Design:
  - HOST MOTIFS of a sample = motifStrings modified (fraction >= 0.4) on any Mp contig of the
    sample, collapsing reverse-complement duplicates so each recognition site is one motif.
  - Each Methanoperedens strain is counted INDEPENDENTLY: Mp contigs are grouped by their
    BORG_Ref lineage (a sample can carry several Mp strains), a strain's contigs pooled, and
    only strains with >= MIN_HOST_LEN of sequence are kept (so density is stable).
  - The Borg of each sample is the comparison element (mini-Borg is ignored here).
  - Motif sites are counted on both strands, IUPAC-aware (nt_search); palindromes counted once.

Output (to /home/shuaiw/MODIFI/tmp/rev_figs/borg/):
  - borg_hostmotif_density_scatter.pdf          : per-sample scatter, Mp strain vs Borg (motifs/kb)
  - borg_hostmotif_density.per_contig.csv        : cached per-(contig,motif) site counts
  - borg_hostmotif_density.sourcedata.csv        : per-(sample,motif,unit) densities behind the figure

Run:
  # first pass counts sequences (~10 min, disk-bound), then caches per-contig counts
  /home/shuaiw/miniconda3/envs/modifi/bin/python motif_density_hostmotifs_per_sample.py --threads 6
  # re-plot instantly from the cache (no recount)
  /home/shuaiw/miniconda3/envs/modifi/bin/python motif_density_hostmotifs_per_sample.py --from_counts
"""
import os, sys, csv, argparse, subprocess
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- paths ----
PROFILE_DIR = "/home/shuaiw/borg/paper/borg_data/profile4"
PROFILE_CSV = os.path.join(PROFILE_DIR, "profile_profile_df_filtered.csv")
SAMPLE_MAP  = os.path.join(PROFILE_DIR, "profile_profile_df_filtered.sample.name.csv")
FASTA_DIR   = "/home/shuaiw/borg/paper/curated_genome/unique"
OUT_DIR     = "/home/shuaiw/MODIFI/tmp/rev_figs/borg"

FRACTION_CUTOFF = 0.4       # motif is a host motif if modified at >= this fraction on an Mp contig
MIN_LEN         = 50000     # a contig needs >= this many bp to be a data point (density stability)

HOST_COLOR = "#ffbf72"
BORG_COLOR = "#3df3f7"

USABLE_SAMPLES = [
    "SR-VP_07_25_2022_A1_60cm",
    "SR-VP_07_25_2022_A1_80cm",
    "SR-VP_07_25_2022_A1_90cm",
    "SR-VP_07_25_2022_A1_100cm",
    "SR-VP_07_25_2022_A1_115cm",
    "SR-VP_9_9_2021_34_2B_1_4m",
]


def load_metadata():
    contig2sample = {}
    for row in csv.DictReader(open(SAMPLE_MAP)):
        contig2sample[row["contig"]] = row["sample"]
    contig2genome, contig2borg = {}, {}
    mp_motif_frac = defaultdict(float)            # (sample, motif_seq) -> max fraction on Mp
    for row in csv.DictReader(open(PROFILE_CSV)):
        c = row["contig"]
        contig2genome.setdefault(c, row["Genome"])
        contig2borg.setdefault(c, row["BORG_Ref"])
    for row in csv.DictReader(open(PROFILE_CSV)):
        c = row["contig"]
        if contig2genome.get(c) != "Mp":
            continue
        s = contig2sample.get(c)
        if s is None:
            continue
        motif_seq = row["motifString"].rsplit("_", 1)[0]
        try:
            frac = float(row["fraction"])
        except ValueError:
            frac = 0.0
        key = (s, motif_seq)
        mp_motif_frac[key] = max(mp_motif_frac[key], frac)
    return contig2sample, contig2genome, contig2borg, mp_motif_frac


def dedupe_revcomp(motifs):
    seen = {}
    for m in motifs:
        canon = min(m, str(Seq(m).reverse_complement()))
        seen.setdefault(canon, m)
    return sorted(set(seen.values()))


def build_contig_fasta_map():
    c2fa, c2len = {}, {}
    for fn in os.listdir(FASTA_DIR):
        if fn.endswith(".fa.fai"):
            fa = os.path.join(FASTA_DIR, fn[:-4])
            for line in open(os.path.join(FASTA_DIR, fn)):
                p = line.split("\t")
                c2fa[p[0]] = fa
                c2len[p[0]] = int(p[1])
    return c2fa, c2len


def extract_seqs(contigs, c2fa):
    by_fa = defaultdict(list)
    for c in contigs:
        fa = c2fa.get(c)
        if fa is None:
            sys.stderr.write(f"WARN: contig not in any fasta index: {c}\n")
            continue
        by_fa[fa].append(c)
    seqs = {}
    for fa, cs in by_fa.items():
        out = subprocess.run(["samtools", "faidx", fa, *cs],
                             capture_output=True, text=True, check=True).stdout
        name, buf = None, []
        for line in out.splitlines():
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def count_sites(seq_str, motif):
    """Double-stranded recognition-site count (IUPAC-aware); palindromes counted once."""
    seq_str = seq_str.upper()
    fwd = len(nt_search(seq_str, motif)) - 1
    rc = str(Seq(motif).reverse_complement())
    if rc == motif:
        return fwd
    return fwd + (len(nt_search(seq_str, rc)) - 1)


def count_sample(sample, host_motifs, contigs, c2fa, c2len):
    """Count each host motif on each contig (Mp + Borg) of the sample. Returns per-contig rows."""
    seqs = extract_seqs(contigs, c2fa)
    rows = []
    for c in contigs:
        if c not in seqs:
            continue
        for motif in host_motifs:
            rows.append((sample, c, motif, count_sites(seqs[c], motif), c2len.get(c, len(seqs[c]))))
    return rows


def _worker(args):
    return count_sample(*args)


def compute_counts(samples, threads):
    contig2sample, contig2genome, contig2borg, mp_motif_frac = load_metadata()
    c2fa, c2len = build_contig_fasta_map()
    sample_hostmotifs, tasks = {}, []
    for s in samples:
        motifs = [m for (ss, m), fr in mp_motif_frac.items() if ss == s and fr >= FRACTION_CUTOFF]
        sample_hostmotifs[s] = dedupe_revcomp(motifs)
        # Mp + Borg contigs only (mini-Borg ignored)
        contigs = [c for c, ss in contig2sample.items()
                   if ss == s and contig2genome.get(c) in ("Mp", "Borg")]
        tasks.append((s, sample_hostmotifs[s], contigs, c2fa, c2len))
        n_mp = sum(1 for c in contigs if contig2genome.get(c) == "Mp")
        n_bg = sum(1 for c in contigs if contig2genome.get(c) == "Borg")
        print(f"{s}: {len(sample_hostmotifs[s])} host motifs; Mp contigs={n_mp}, Borg contigs={n_bg}")
    if threads > 1:
        import multiprocessing as mp
        with mp.Pool(min(threads, len(tasks))) as pool:
            results = pool.map(_worker, tasks)
    else:
        results = [_worker(t) for t in tasks]
    rows = [r for rs in results for r in rs]
    df = pd.DataFrame(rows, columns=["sample", "contig", "motif", "sites", "length_bp"])
    df["genome"] = df["contig"].map(contig2genome)
    df["borg_ref"] = df["contig"].map(contig2borg)
    return df


def strain_label(borg_ref):
    """Short label for an Mp strain from its BORG_Ref lineage."""
    b = borg_ref.replace("<GTDB>", "").strip()
    if b.startswith("Methanoperedens sp"):
        return "Mp " + b.split("sp")[-1].strip()
    if b == "Methanoperedens":
        return "Mp (unresolved)"
    if b.startswith("Methanoperedens_44_19-type"):
        return "Mp 44_19-type"
    return b[:24]


def aggregate(df):
    """Per-CONTIG densities (no pooling): each Mp contig and each Borg contig is an
    independent genome unit, since a sample can carry several Mp strains and several Borg
    genomes. Contigs shorter than MIN_LEN are dropped so densities stay stable."""
    rows = []
    for _, r in df.iterrows():
        if r["length_bp"] < MIN_LEN or r["genome"] not in ("Mp", "Borg"):
            continue
        rows.append(dict(sample=r["sample"], genome=r["genome"], contig=r["contig"],
                         strain=strain_label(r["borg_ref"]) if r["genome"] == "Mp" else "Borg",
                         motif=r["motif"], sites=r["sites"], length_bp=r["length_bp"],
                         density_per_kb=r["sites"] / (r["length_bp"] / 1000.0)))
    return pd.DataFrame(rows)


def make_scatter(agg, out_pdf):
    """Per-sample scatter. Each point = one (Mp contig, Borg contig, host motif): x = that
    Mp contig's motif density, y = that Borg contig's density. Points are colored by the
    full (original) Mp host contig name; Borg contigs get distinct marker shapes when a
    sample has more than one Borg genome."""
    samples = [s for s in USABLE_SAMPLES if s in set(agg["sample"])]
    ncol = 3
    nrow = int(np.ceil(len(samples) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 4.6, nrow * 4.4), squeeze=False)
    cmap = plt.get_cmap("tab20")
    markers = ["o", "s", "^", "D", "v", "P"]
    short = lambda c: "_".join(c.split("_")[-2:])   # e.g. ..._METAMDBG_723659_L -> 723659_L
    for i, s in enumerate(samples):
        ax = axes[i // ncol][i % ncol]
        sub = agg[agg["sample"] == s]
        mp = sub[sub["genome"] == "Mp"]
        bg = sub[sub["genome"] == "Borg"]
        borg_contigs = sorted(bg["contig"].unique())
        mp_contigs = sorted(mp["contig"].unique())
        maxv = 0.0
        for bk, bc in enumerate(borg_contigs):
            bdens = bg[bg["contig"] == bc].set_index("motif")["density_per_kb"].to_dict()
            for k, mc in enumerate(mp_contigs):
                hc = mp[mp["contig"] == mc]
                xs, ys = [], []
                for _, r in hc.iterrows():
                    y = bdens.get(r["motif"])
                    if y is None:
                        continue
                    xs.append(r["density_per_kb"]); ys.append(y)
                if xs:
                    lab = short(mc) if len(borg_contigs) == 1 else f"{short(mc)}/B{bk+1}"
                    ax.scatter(xs, ys, s=22, alpha=0.75, color=cmap(k % 20),
                               marker=markers[bk % len(markers)],
                               edgecolor="black", linewidth=0.2, label=lab)
                    maxv = max(maxv, max(xs + ys))
        lim = maxv * 1.05 if maxv > 0 else 1
        ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        n_mp = mp["contig"].nunique(); n_bg = len(borg_contigs)
        ax.set_title(f"{s}  ({n_mp} Mp contigs, {n_bg} Borg)", fontsize=9)
        ax.set_xlabel("Mp host contig (motifs/kb)", fontsize=8)
        ax.set_ylabel("Borg contig (motifs/kb)", fontsize=8)
        ax.tick_params(labelsize=7)
        ncol_leg = 2 if n_mp > 6 else 1
        ax.legend(fontsize=6, frameon=False, loc="upper left", ncol=ncol_leg,
                  handletextpad=0.2, columnspacing=0.6, title="Mp contig", title_fontsize=6)
    for k in range(len(samples), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle(r"Host-motif sequence density: each $\it{Methanoperedens}$ host contig vs each Borg contig"
                 "\n(one point per host motif; color = Mp host contig; dashed line y = x; points on x-axis = absent in Borg)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved scatter figure: {out_pdf}")


def make_boxplot(agg, out_pdf):
    """Pool all samples: for each (sample, host motif) pair, host density = pooled Mp density
    (sum sites / sum length across the sample's Mp contigs) and Borg density = pooled Borg
    density. Box plot of the two groups (paired by sample+motif) with a Wilcoxon signed-rank
    p-value (Borg vs host)."""
    from scipy.stats import wilcoxon
    pairs = []
    for (sample, motif), g in agg.groupby(["sample", "motif"]):
        mp = g[g["genome"] == "Mp"]; bg = g[g["genome"] == "Borg"]
        if mp.empty or bg.empty:
            continue
        h = mp["sites"].sum() / (mp["length_bp"].sum() / 1000.0)
        b = bg["sites"].sum() / (bg["length_bp"].sum() / 1000.0)
        pairs.append((sample, motif, h, b))
    pdf = pd.DataFrame(pairs, columns=["sample", "motif", "host", "borg"])
    stat, p = wilcoxon(pdf["host"], pdf["borg"])           # paired, two-sided
    ratios = (pdf["borg"] / pdf["host"]).replace([np.inf], np.nan).dropna()
    med_ratio = float(np.median(ratios))

    fig, ax = plt.subplots(figsize=(4.2, 5.2))
    data = [pdf["host"].values, pdf["borg"].values]
    bp = ax.boxplot(data, tick_labels=["Mp host", "Borg"], widths=0.55, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch, col in zip(bp["boxes"], [HOST_COLOR, BORG_COLOR]):
        patch.set_facecolor(col); patch.set_alpha(0.85)
    # paired points + faint connecting lines
    rng = np.random.RandomState(0)
    x1 = 1 + (rng.rand(len(pdf)) - 0.5) * 0.12
    x2 = 2 + (rng.rand(len(pdf)) - 0.5) * 0.12
    for a, b, ya, yb in zip(x1, x2, pdf["host"], pdf["borg"]):
        ax.plot([a, b], [ya, yb], color="grey", lw=0.2, alpha=0.3, zorder=1)
    ax.scatter(x1, pdf["host"], s=10, color=HOST_COLOR, edgecolor="black", linewidth=0.2, zorder=2)
    ax.scatter(x2, pdf["borg"], s=10, color=BORG_COLOR, edgecolor="black", linewidth=0.2, zorder=2)
    ax.set_ylabel("host-motif sequence density (motifs/kb)", fontsize=9)
    ax.set_title("Host-motif density: Borg vs host\n"
                 f"(all samples pooled, n={len(pdf)} sample-motif pairs)", fontsize=10)
    ptxt = "p < 1e-300" if p == 0 else f"p = {p:.2e}"
    ymax = max(pdf["host"].max(), pdf["borg"].max())
    ax.plot([1, 2], [ymax * 1.04, ymax * 1.04], color="black", lw=0.8)
    ax.text(1.5, ymax * 1.06, f"Wilcoxon {ptxt}\nmedian Borg/host = {med_ratio:.2f}",
            ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, ymax * 1.18)
    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved box plot: {out_pdf}  (Wilcoxon {ptxt}, n={len(pdf)}, median Borg/host={med_ratio:.2f})")
    return pdf


def compute_pairs(agg):
    """One paired (host, Borg) density per (sample, host motif); host = pooled Mp density."""
    rows = []
    for (sample, motif), g in agg.groupby(["sample", "motif"]):
        mp = g[g["genome"] == "Mp"]; bg = g[g["genome"] == "Borg"]
        if mp.empty or bg.empty:
            continue
        h = mp["sites"].sum() / (mp["length_bp"].sum() / 1000.0)
        b = bg["sites"].sum() / (bg["length_bp"].sum() / 1000.0)
        rows.append((sample, motif, h, b))
    return pd.DataFrame(rows, columns=["sample", "motif", "host", "borg"])


def _panel_letter(ax, letter):
    ax.text(-0.16, 1.02, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="right")


def make_combined(agg, out_pdf):
    """One figure: (a-f) per-sample scatters + (g) pooled Borg-vs-host box plot.
    Bold panel letters, no titles, caption at the bottom."""
    from scipy.stats import wilcoxon
    samples = [s for s in USABLE_SAMPLES if s in set(agg["sample"])]
    cmap = plt.get_cmap("tab20")
    markers = ["o", "s", "^", "D", "v", "P"]
    short = lambda c: "_".join(c.split("_")[-2:])
    short_sample = lambda s: s.replace("SR-VP_07_25_2022_", "").replace("SR-VP_9_9_2021_", "")

    fig = plt.figure(figsize=(14, 15))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.15], hspace=0.42, wspace=0.28)
    letters = "abcdefghij"

    # (a-f) scatter panels
    for i, s in enumerate(samples):
        ax = fig.add_subplot(gs[i // 3, i % 3])
        sub = agg[agg["sample"] == s]
        mp = sub[sub["genome"] == "Mp"]; bg = sub[sub["genome"] == "Borg"]
        borg_contigs = sorted(bg["contig"].unique())
        mp_contigs = sorted(mp["contig"].unique())
        maxv = 0.0
        for bk, bc in enumerate(borg_contigs):
            bdens = bg[bg["contig"] == bc].set_index("motif")["density_per_kb"].to_dict()
            for k, mc in enumerate(mp_contigs):
                hc = mp[mp["contig"] == mc]
                xs, ys = [], []
                for _, r in hc.iterrows():
                    y = bdens.get(r["motif"])
                    if y is None:
                        continue
                    xs.append(r["density_per_kb"]); ys.append(y)
                if xs:
                    lab = short(mc) if len(borg_contigs) == 1 else f"{short(mc)}/B{bk+1}"
                    ax.scatter(xs, ys, s=20, alpha=0.75, color=cmap(k % 20),
                               marker=markers[bk % len(markers)],
                               edgecolor="black", linewidth=0.2, label=lab)
                    maxv = max(maxv, max(xs + ys))
        lim = maxv * 1.05 if maxv > 0 else 1
        ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("Mp host contig (motifs/kb)", fontsize=8)
        ax.set_ylabel("Borg contig (motifs/kb)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.text(0.97, 0.03, short_sample(s), transform=ax.transAxes,
                fontsize=8, va="bottom", ha="right", color="0.35")
        ncol_leg = 2 if len(mp_contigs) > 6 else 1
        ax.legend(fontsize=5.5, frameon=False, loc="upper left", ncol=ncol_leg,
                  handletextpad=0.2, columnspacing=0.6, title="Mp contig", title_fontsize=6)
        _panel_letter(ax, letters[i])

    # (g) box plot
    pdf = compute_pairs(agg)
    stat, p = wilcoxon(pdf["host"], pdf["borg"])
    ratios = (pdf["borg"] / pdf["host"]).replace([np.inf], np.nan).dropna()
    med_ratio = float(np.median(ratios))
    axb = fig.add_subplot(gs[2, 1])
    data = [pdf["host"].values, pdf["borg"].values]
    bp = axb.boxplot(data, tick_labels=["Mp host", "Borg"], widths=0.55, showfliers=False,
                     patch_artist=True, medianprops=dict(color="black"))
    for patch, col in zip(bp["boxes"], [HOST_COLOR, BORG_COLOR]):
        patch.set_facecolor(col); patch.set_alpha(0.85)
    rng = np.random.RandomState(0)
    x1 = 1 + (rng.rand(len(pdf)) - 0.5) * 0.12
    x2 = 2 + (rng.rand(len(pdf)) - 0.5) * 0.12
    for a, b, ya, yb in zip(x1, x2, pdf["host"], pdf["borg"]):
        axb.plot([a, b], [ya, yb], color="grey", lw=0.2, alpha=0.3, zorder=1)
    axb.scatter(x1, pdf["host"], s=9, color=HOST_COLOR, edgecolor="black", linewidth=0.2, zorder=2)
    axb.scatter(x2, pdf["borg"], s=9, color=BORG_COLOR, edgecolor="black", linewidth=0.2, zorder=2)
    axb.set_ylabel("host-motif density (motifs/kb)", fontsize=8)
    axb.tick_params(labelsize=8)
    ptxt = "p < 1e-300" if p == 0 else f"p = {p:.2e}"
    ymax = max(pdf["host"].max(), pdf["borg"].max())
    axb.plot([1, 2], [ymax * 1.04, ymax * 1.04], color="black", lw=0.8)
    axb.text(1.5, ymax * 1.06, f"Wilcoxon {ptxt}\nmedian Borg/host = {med_ratio:.2f}",
             ha="center", va="bottom", fontsize=8)
    axb.set_ylim(0, ymax * 1.20)
    _panel_letter(axb, letters[len(samples)])

    fig.subplots_adjust(left=0.06, right=0.98, top=0.98, bottom=0.04)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"Saved combined figure: {out_pdf}  (Wilcoxon {ptxt}, n={len(pdf)}, median Borg/host={med_ratio:.2f})")

    # Plain-text caption (printed to screen, NOT drawn on the figure).
    letter_map = "; ".join(f"{letters[i]}, {short_sample(s)}" for i, s in enumerate(samples))
    ab = f"{letters[0]}-{letters[len(samples)-1]}"
    gl = letters[len(samples)]
    caption = (
        "Host-motif sequence density in Methanoperedens hosts versus their Borgs. "
        "Host motifs are the modification motifs called on Methanoperedens (Mp) contigs of each sample "
        "(modification fraction >= 0.4), reverse-complement duplicates merged. Motif recognition sites "
        "were counted on both strands (IUPAC-aware) in each contig >= 50 kb; density = sites per kb. "
        f"({ab}) Per-sample scatter of each Mp host contig (colors) versus the sample's Borg contig; "
        "each point is one host motif; dashed line, y = x; points on the x-axis are motifs absent in the "
        f"Borg. Samples: {letter_map}. ({gl}) All samples pooled: host-motif density in the Mp host "
        f"versus the Borg, paired by sample and motif (n = {len(pdf)}); grey lines connect paired values; "
        f"two-sided Wilcoxon signed-rank {ptxt}; median Borg/host density ratio = {med_ratio:.2f}."
    )
    print("\n===== FIGURE CAPTION (Methanoperedens in italics) =====\n")
    print(caption)
    print()
    return pdf


def print_summary(agg):
    print("\n=== Summary: host motifs present in the Borg (per contig; each Mp contig vs each Borg contig) ===")
    for sample, sub in agg.groupby("sample"):
        mp = sub[sub["genome"] == "Mp"]
        bg = sub[sub["genome"] == "Borg"]
        ratios, absent, tot = [], 0, 0
        for bc in bg["contig"].unique():
            bdens = bg[bg["contig"] == bc].set_index("motif")["density_per_kb"].to_dict()
            for _, r in mp.iterrows():
                y = bdens.get(r["motif"])
                if y is None:
                    continue
                tot += 1
                if y == 0:
                    absent += 1
                if r["density_per_kb"] > 0:
                    ratios.append(y / r["density_per_kb"])
        med = np.median(ratios) if ratios else float("nan")
        print(f"{sample}: {mp['contig'].nunique()} Mp contigs x {bg['contig'].nunique()} Borg; "
              f"{absent}/{tot} contig-motif pairs absent in Borg; "
              f"median Borg/host density ratio = {med:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*", default=USABLE_SAMPLES)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out_dir", default=OUT_DIR)
    ap.add_argument("--from_counts", action="store_true",
                    help="reuse cached per-contig counts instead of recounting sequences")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    counts_csv = os.path.join(args.out_dir, "borg_hostmotif_density.per_contig.csv")

    if args.from_counts and os.path.exists(counts_csv):
        df = pd.read_csv(counts_csv)
        df = df[df["sample"].isin(args.samples)]
        print(f"Loaded cached per-contig counts: {counts_csv} ({len(df)} rows)")
    else:
        df = compute_counts(args.samples, args.threads)
        df.to_csv(counts_csv, index=False)
        print(f"Saved per-contig counts: {counts_csv} ({len(df)} rows)")

    agg = aggregate(df)
    sd = os.path.join(args.out_dir, "borg_hostmotif_density.sourcedata.csv")
    agg.to_csv(sd, index=False)
    print(f"Saved source data: {sd} ({len(agg)} rows)")

    paired = make_combined(agg, os.path.join(args.out_dir, "borg_hostmotif_density_combined.pdf"))
    paired.to_csv(os.path.join(args.out_dir, "borg_hostmotif_density_boxplot.sourcedata.csv"), index=False)
    print_summary(agg)


if __name__ == "__main__":
    main()
