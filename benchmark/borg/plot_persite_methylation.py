#!/usr/bin/env python
"""
Figure-6b-style per-site modification map for one motif across selected contigs of a sample.

For a given motif (here ACC, methylated base at position 1), every occurrence is drawn as a
triangle along each contig's genomic coordinate:
    triangle up   (^) = forward strand (+)
    triangle down (v) = reverse strand (-)
    filled            = modified (GFF modified_base call with score >= 30 at that base)
    hollow            = unmodified
One track (subplot) per contig; x-axis in kb.

Reuses the enumeration logic from benchmark/circos/get_circos_data.py:get_motif_sites and the
GFF parse from get_modified_ratio. Per-occurrence tag:
    forward: f"{contig}:{site+exact_pos}+"          (1-based methylated position)
    reverse: f"{contig}:{site+ (len-exact_pos+1)}-"
modified iff tag in modified_loci.

Run:
  /home/shuaiw/miniconda3/envs/methy3/bin/python benchmark/borg/plot_persite_methylation.py
"""
import os
import csv
import re
import random
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio.Seq import Seq
from Bio import SeqIO

# ---------------- config ----------------
SAMPLE = "SR-VP_07_25_2022_A1_100cm"
METHY_DIR = "/home/shuaiw/borg/paper/gg_run2/soil_100/soil_100_methylation4"
OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/borg"
SCORE_CUTOFF = 30

MOTIF = "ACC"      # motif string
EXACT_POS = 1      # 1-based position of the methylated base within the motif
MOTIF_ID = f"{MOTIF}_{EXACT_POS}"

WINDOW = 20_000    # show a single random WINDOW-bp window per contig (same width for all)
random.seed(42)    # reproducible window choice

MOTIF_COLOR = "#2c7fb8"   # per-motif color (fill = modified, hollow = unmodified)
TRACK_BG = "#eef3f8"      # light tinted background band per track
LEG_DARK = "#333333"      # neutral colour for the (status/strand) legend keys

def contig_path(name, sub, ext):
    return os.path.join(METHY_DIR, sub, f"{name}.{ext}")

# label, type, contig id; fa/gff resolved from the methylation dir unless given explicitly.
# The mini_chr (FINAL_272kb) has no per-contig fa/gff in this pipeline; fill fa/gff when found.
CONTIGS = [
    {"label": "Mp_724567_L",       "type": "Mp",        "name": f"{SAMPLE}_PACBIO-HIFI_METAMDBG_724567_L"},
    {"label": "Mp_724408_L",       "type": "Mp",        "name": f"{SAMPLE}_PACBIO-HIFI_METAMDBG_724408_L"},
    {"label": "Mp_724941_L",       "type": "Mp",        "name": f"{SAMPLE}_PACBIO-HIFI_METAMDBG_724941_L"},
    # The curated mini-chromosome was not split per-contig in gg_run2; it was processed in
    # gg_run3 (same sample, soil_100), so its fa/gff are taken from there.
    {"label": "Mini_Chr_272kb",    "type": "Mini_Chr",  "name": "FINAL_272kb_HMp-Plasmid-like_SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_39_43_complete",
        "fa":  "/home/shuaiw/borg/paper/gg_run3/soil_100/soil_100_methylation4/contigs/FINAL_272kb_HMp-Plasmid-like_SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_39_43_complete.fa",
        "gff": "/home/shuaiw/borg/paper/gg_run3/soil_100/soil_100_methylation4/gffs/FINAL_272kb_HMp-Plasmid-like_SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_39_43_complete.reprocess.gff"},
    {"label": "Iris_Borg",            "type": "Borg",      "name": f"{SAMPLE}_PACBIO-HIFI_METAMDBG_723659_L"},
    {"label": "Uranus_mini-Borg_69kb","type": "mini_Borg", "name": f"{SAMPLE}_PACBIO-HIFI_METAMDBG_717158_L"},
]

# ---------------- reused logic ----------------
def read_ref(fa):
    rec = next(SeqIO.parse(fa, "fasta"))
    return rec.id, str(rec.seq)

def read_modified_loci(gff, score_cutoff=SCORE_CUTOFF):
    loci = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            ref, pos, score, strand = p[0], int(p[3]), int(p[5]), p[6]
            if score < score_cutoff:
                continue
            loci[f"{ref}:{pos}{strand}"] = score
    return loci

def enumerate_sites(ref_id, seq, motif, exact_pos, modified_loci):
    """Return list of (pos1based, strand, modified_bool, score) for every occurrence."""
    motif_len = len(motif)
    rev_exact_pos = motif_len - exact_pos + 1
    revcomp = str(Seq(motif).reverse_complement())
    sites = []
    counts = dict(for_loci=0, for_mod=0, rev_loci=0, rev_mod=0)
    # forward: occurrences of the motif
    for m in re.finditer(f"(?={motif})", seq):   # overlapping matches
        site = m.start()
        pos = site + exact_pos                    # 1-based methylated base
        tag = f"{ref_id}:{pos}+"
        modified = tag in modified_loci
        sites.append((pos, "+", modified, modified_loci.get(tag, 0)))
        counts["for_loci"] += 1
        counts["for_mod"] += int(modified)
    # reverse: occurrences of the reverse complement
    for m in re.finditer(f"(?={revcomp})", seq):
        site = m.start()
        pos = site + rev_exact_pos
        tag = f"{ref_id}:{pos}-"
        modified = tag in modified_loci
        sites.append((pos, "-", modified, modified_loci.get(tag, 0)))
        counts["rev_loci"] += 1
        counts["rev_mod"] += int(modified)
    return sites, counts

# ---------------- gather ----------------
os.makedirs(OUT_DIR, exist_ok=True)
tracks = []
source_rows = []
for c in CONTIGS:
    fa = c.get("fa", contig_path(c["name"], "contigs", "fa"))
    gff = c.get("gff", contig_path(c["name"], "gffs", "reprocess.gff"))
    if not fa or not gff or not (os.path.exists(fa) and os.path.exists(gff)):
        print(f"[skip] {c['label']} ({c['name']}): fa/gff not found "
              f"(fa={fa}, gff={gff})")
        continue
    ref_id, seq = read_ref(fa)
    modified_loci = read_modified_loci(gff)
    sites, counts = enumerate_sites(ref_id, seq, MOTIF, EXACT_POS, modified_loci)
    length = len(seq)
    # pick a random WINDOW-bp window (whole contig if shorter than WINDOW)
    win_start = 0 if length <= WINDOW else random.randint(0, length - WINDOW)
    win_end = win_start + WINDOW
    wsites = [s for s in sites if win_start < s[0] <= win_end]
    wmod = sum(1 for _, _, m, _ in wsites if m)
    print(f"{c['label']:16s} len={length:>9,} bp | full for {counts['for_mod']}/{counts['for_loci']} "
          f"rev {counts['rev_mod']}/{counts['rev_loci']} | "
          f"window {win_start/1000:.0f}-{win_end/1000:.0f} kb: {wmod}/{len(wsites)} modified")
    tracks.append({**c, "ref_id": ref_id, "length": length, "sites": wsites,
                   "win_start": win_start, "win_end": win_end})
    for pos, strand, modified, score in wsites:
        source_rows.append([c["label"], c["type"], ref_id, win_start, win_end, pos, strand,
                            "modified" if modified else "unmodified", score])

# ---------------- source data ----------------
sd = os.path.join(OUT_DIR, f"persite_{MOTIF_ID}_{SAMPLE}.sourcedata.csv")
with open(sd, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["track", "type", "contig", "window_start", "window_end",
                "position_1based", "strand", "status", "gff_score"])
    w.writerows(source_rows)
print(f"Saved source data: {sd}  ({len(source_rows)} sites)")

# ---------------- plot (Fig-6b style) ----------------
# One baseline per contig; forward triangles above the line, reverse below the SAME line.
# Fill encodes status (solid = modified, hollow = unmodified); colour encodes the motif.
# Light tinted background band per track; legend split STATUS x STRAND.
n = len(tracks)
fig, axes = plt.subplots(n, 1, figsize=(12, 0.95 * n + 2.2), squeeze=False)
axes = axes[:, 0]
YF, YR = 0.34, -0.34   # forward just above / reverse just below the baseline

for ax, t in zip(axes, tracks):
    def xs(strand, mod):
        return [pos / 1000.0 for pos, s, m, sc in t["sites"] if s == strand and m == mod]
    ax.set_facecolor(TRACK_BG)
    ax.axhline(0, color="#9aa7b4", lw=0.8, zorder=1)          # genomic baseline
    # unmodified: hollow (outline only)
    ax.scatter(xs("+", False), [YF] * len(xs("+", False)), marker="^", s=26,
               facecolors="none", edgecolors=MOTIF_COLOR, linewidths=0.5, alpha=0.55, zorder=2)
    ax.scatter(xs("-", False), [YR] * len(xs("-", False)), marker="v", s=26,
               facecolors="none", edgecolors=MOTIF_COLOR, linewidths=0.5, alpha=0.55, zorder=2)
    # modified: solid fill
    ax.scatter(xs("+", True), [YF] * len(xs("+", True)), marker="^", s=30,
               facecolors=MOTIF_COLOR, edgecolors=MOTIF_COLOR, linewidths=0.3, alpha=0.98, zorder=3)
    ax.scatter(xs("-", True), [YR] * len(xs("-", True)), marker="v", s=30,
               facecolors=MOTIF_COLOR, edgecolors=MOTIF_COLOR, linewidths=0.3, alpha=0.98, zorder=3)
    nmod = sum(1 for _, _, m, _ in t["sites"] if m)
    ntot = len(t["sites"])
    ax.set_ylim(-1.0, 1.0)
    ax.set_xlim(t["win_start"] / 1000.0, t["win_end"] / 1000.0)
    ax.set_yticks([])
    for sp in ("left", "right", "top", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="x", length=3, labelsize=8)
    ax.set_ylabel(f"{t['label']}\n{t['length']/1000:.0f} kb",
                  rotation=0, ha="right", va="center", fontsize=8)
    ax.annotate(f"{nmod}/{ntot} modified", xy=(1, 1), xycoords="axes fraction",
                xytext=(-3, -3), textcoords="offset points", ha="right", va="top",
                fontsize=6, color="#666666")

axes[-1].set_xlabel("genomic coordinate (kb)", fontsize=9)

# motif header (colored square + name), Fig-6b style, top-left
fig.text(0.13, 0.975, "■", color=MOTIF_COLOR, fontsize=13, va="center")
fig.text(0.145, 0.975, f"{MOTIF} (methylated position {EXACT_POS})",
         fontsize=11, fontweight="bold", va="center")
fig.text(0.99, 0.975, f"{SAMPLE}  -  random {WINDOW//1000} kb window per contig",
         fontsize=8, color="#666666", ha="right", va="center")

# legend: STATUS (fill) x STRAND (triangle direction), neutral colour like the app
status_mod = Line2D([0], [0], marker="s", color="none", markerfacecolor=LEG_DARK,
                    markeredgecolor=LEG_DARK, markersize=9, label="Modified (methylated)")
status_unmod = Line2D([0], [0], marker="s", color="none", markerfacecolor="none",
                      markeredgecolor=LEG_DARK, markersize=9, label="Unmodified")
strand_fwd = Line2D([0], [0], marker="^", color="none", markerfacecolor=LEG_DARK,
                    markeredgecolor=LEG_DARK, markersize=9, label="+ (forward)")
strand_rev = Line2D([0], [0], marker="v", color="none", markerfacecolor=LEG_DARK,
                    markeredgecolor=LEG_DARK, markersize=9, label="- (reverse)")
leg1 = fig.legend(handles=[status_mod, status_unmod], title="STATUS", ncol=2,
                  loc="upper center", bbox_to_anchor=(0.42, 0.945), frameon=False,
                  fontsize=8, title_fontsize=8, handletextpad=0.3, columnspacing=1.2)
fig.add_artist(leg1)
fig.legend(handles=[strand_fwd, strand_rev], title="STRAND", ncol=2,
           loc="upper center", bbox_to_anchor=(0.72, 0.945), frameon=False,
           fontsize=8, title_fontsize=8, handletextpad=0.3, columnspacing=1.2)

fig.tight_layout(rect=[0.0, 0.0, 1, 0.87])
out_pdf = os.path.join(OUT_DIR, f"persite_{MOTIF_ID}_{SAMPLE}.pdf")
out_png = os.path.join(OUT_DIR, f"persite_{MOTIF_ID}_{SAMPLE}.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"Saved figure: {out_pdf}")
