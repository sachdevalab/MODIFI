#!/usr/bin/env python3
"""
MODIFI vs SOTA methylation callers — MOTIF-LEVEL comparison on 10 circular soil
genomes (test_100). Consolidated: per-contig GFF building, per-genome motif-set
aggregation, runtime summary, and figures (value/area-proportional Venns).

Motif discovery is done PER CONTIG (per genome) — bacterial methylation motifs are
genome-specific, so pooling genomes would dilute each MTase's signal. Per tool we
then take the UNION of motifs across the 10 genomes.

  - MODIFI-HiFi : its own native de-novo motifs (per-contig motifs/*.csv)
  - ipdSummary / fibertools / jasmine : motifMaker run per-contig on their calls
    (they do not discover motifs themselves)

Comparison is centred on MODIFI-HiFi vs each SOTA tool.

Usage:
  python motif_comparison.py gffs      # build per-contig GFFs for ipd/ft/jasmine
  python motif_comparison.py analyze   # aggregate motifs, Venns, table, runtime
"""
from __future__ import annotations
import csv, gzip, math, os, re, sys
from pathlib import Path

OUT = Path("/home/shuaiw/borg/paper/ipdsummary/compare_all")
FIG = Path("/home/shuaiw/MODIFI/tmp/rev_figs/compare_all")
MC = OUT / "motifs_compare"
PC = MC / "percontig"                     # per-contig GFFs + motifMaker outputs
PREF = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META"
CONTIGS = ["658", "327", "645", "479", "396", "514", "504", "517", "350", "316"]
CONTIG_NAMES = [f"{PREF}_{n}_C" for n in CONTIGS]
REF_DIR = OUT / "refs"                     # per-contig fastas <contig>.fasta (indexed)
MODIFI_NATIVE = OUT / "modifi_hifi/selected10/motifs"   # MODIFI per-contig native motifs

IPD_SUBREAD_GFF = Path("/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out/test_100.gff")
FT_PILEUP = OUT / "fibertools/selected10.m6a.pileup.bed"
JAS_BED = OUT / "jasmine/selected10.cpg.combined.bed.gz"

# motifMaker tools (not MODIFI): per-contig GFF dir + min-score
MM_TOOLS = {"ipdSummary": 30, "fibertools": 50, "jasmine": 50}
FT_FRAC, JAS_PCT, MINCOV = 0.50, 50.0, 4

T_IPD = Path("/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out/test_100.ipdSummary.time")
T_MODIFI_SUB = Path("/home/shuaiw/borg/paper/ipdsummary/soil_1/modifi.out/test_100/modifi.host.time")

# ---------------- fasta ----------------
def load_ref():
    seqs, name, buf = {}, None, []
    with open(REF_DIR / "selected10.fasta") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs

# ---------------- per-contig GFF building ----------------
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def _revcomp(s): return s.translate(_COMP)[::-1]

def _context41(seq, pos1, strand):
    i = pos1 - 1
    lo, hi = i - 20, i + 21
    ctx = "N" * max(0, -lo) + seq[max(0, lo):min(len(seq), hi)] + "N" * max(0, hi - len(seq))
    return _revcomp(ctx) if strand == "-" else ctx

def _gff_line(seq, src, pos1, score, strand, cov, context):
    return (f"{seq}\t{src}\tmodified_base\t{pos1}\t{pos1}\t{score}\t{strand}\t.\t"
            f"coverage={cov};context={context};IPDRatio=2.0\n")

def build_gffs():
    for t in MM_TOOLS:
        (PC / t).mkdir(parents=True, exist_ok=True)
    ref = load_ref()

    # ipdSummary: split the whole subread GFF by contig (it already carries context=)
    handles = {c: open(PC / "ipdSummary" / f"{c}.gff", "w") for c in CONTIG_NAMES}
    want = set(CONTIG_NAMES)
    with open(IPD_SUBREAD_GFF) as f:
        for line in f:
            if line[0] == "#":
                continue
            seq = line.split("\t", 1)[0]
            if seq in want:
                handles[seq].write(line)
    for h in handles.values():
        h.close()

    # fibertools 6mA: per-contig, strand from ref base, add context
    ftc = {c: 0 for c in CONTIG_NAMES}
    fh = {c: open(PC / "fibertools" / f"{c}.gff", "w") for c in CONTIG_NAMES}
    if FT_PILEUP.is_file():
        for line in open(FT_PILEUP):
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            seq = p[0]
            if seq not in ref:
                continue
            cov = int(p[3]); m6a = int(p[8])
            if cov < MINCOV or cov == 0 or m6a / cov < FT_FRAC:
                continue
            pos1 = int(p[1]) + 1
            base = ref[seq][pos1 - 1]
            strand = "+" if base == "A" else "-" if base == "T" else None
            if strand is None:
                continue
            fh[seq].write(_gff_line(seq, "fibertools", pos1, round(m6a / cov * 100, 1),
                                    strand, cov, _context41(ref[seq], pos1, strand))); ftc[seq] += 1
    for h in fh.values():
        h.close()

    # jasmine 5mC/CpG: per-contig, strand from ref base, add context
    jc = {c: 0 for c in CONTIG_NAMES}
    jh = {c: open(PC / "jasmine" / f"{c}.gff", "w") for c in CONTIG_NAMES}
    if JAS_BED.is_file():
        with gzip.open(JAS_BED, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[4] != "Total":
                    continue
                seq = p[0]
                if seq not in ref:
                    continue
                try:
                    meth = float(p[3]); cov = int(p[5])
                except ValueError:
                    continue
                if cov < MINCOV or meth < JAS_PCT:
                    continue
                pos1 = int(p[1]) + 1
                base = ref[seq][pos1 - 1]
                strand = "+" if base == "C" else "-" if base == "G" else None
                if strand is None:
                    continue
                jh[seq].write(_gff_line(seq, "jasmine", pos1, round(meth, 1), strand, cov,
                                        _context41(ref[seq], pos1, strand))); jc[seq] += 1
    for h in jh.values():
        h.close()

    for c in CONTIG_NAMES:
        n = c.split("META_")[1]
        ip = sum(1 for _ in open(PC / "ipdSummary" / f"{c}.gff"))
        print(f"  {n:8s} ipd={ip:>7}  fibertools={ftc[c]:>6}  jasmine={jc[c]:>5}")

# ---------------- motif-set aggregation (union across the 10 genomes) ----------------
# MODIFI default motif-retention criteria (main.py defaults) + RC collapse
from Bio.Seq import Seq
MIN_FRAC = 0.4
MIN_SITES = 30

def canon(m):
    """RC-collapse: a motif and its reverse complement map to one canonical
    string (palindromes map to themselves). Same Bio.Seq dependency as
    scripts/derep_motifs.py."""
    m = m.strip().upper()
    try:
        rc = str(Seq(m).reverse_complement())
    except Exception:
        rc = m
    return min(m, rc)

def _add(motifs, info, m, row):
    m = (m or "").strip().upper()
    if not m:
        return
    try:
        fr = float(row.get("fraction", 0) or 0)
        nd = int(float(row.get("nDetected", 0) or 0))
    except ValueError:
        return
    if fr < MIN_FRAC or nd < MIN_SITES:        # MODIFI default filter, all tools
        return
    cm = canon(m)                               # collapse reverse-complement pairs
    motifs.add(cm)
    if cm not in info or fr > float(info[cm].get("fraction", 0) or 0):
        info[cm] = row

def load_modifi_native():
    """Union of MODIFI-HiFi native motifs across per-contig motifs/*.csv."""
    motifs, info = set(), {}
    for p in sorted(MODIFI_NATIVE.glob("*.motifs.csv")):
        for row in csv.DictReader(open(p)):
            _add(motifs, info, row.get("motifString"), row)
    return motifs, info

def aggregate_motifmaker(tool):
    """Union of a tool's per-contig motifMaker motifs across the 10 genomes."""
    motifs, info = set(), {}
    for c in CONTIG_NAMES:
        p = PC / tool / f"{c}.motifs.csv"
        if p.is_file():
            for row in csv.DictReader(open(p)):
                _add(motifs, info, row.get("motifString"), row)
    return motifs, info

# ---------------- runtime ----------------
def parse_time(path):
    if not Path(path).is_file() or Path(path).stat().st_size == 0:
        return None
    t = Path(path).read_text(errors="replace")
    um = re.search(r"User time \(seconds\):\s*([\d.]+)", t)
    sm = re.search(r"System time \(seconds\):\s*([\d.]+)", t)
    mm = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", t)
    em = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", t)
    if not (um and sm and mm and em):
        return None
    parts = em.group(1).split(":")
    wall = (int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])) if len(parts) == 3 \
        else int(parts[0]) * 60 + float(parts[1])
    return dict(cpu_hr=(float(um.group(1)) + float(sm.group(1))) / 3600,
                wall_hr=wall / 3600, rss_gb=int(mm.group(1)) / 1024 / 1024)

def sum_times(paths):
    cpu = wall = rss = 0.0; n = 0
    for p in paths:
        r = parse_time(p)
        if r:
            cpu += r["cpu_hr"]; wall += r["wall_hr"]; rss = max(rss, r["rss_gb"]); n += 1
    return dict(cpu_hr=cpu, wall_hr=wall, rss_gb=rss, n=n)

# ---------------- UpSet plot ----------------
def draw_upset(fig, sets, names, colors):
    """UpSet plot: intersection-size bars (top) + membership dot-matrix (bottom)
    + per-tool totals (left). Unambiguous for 4+ sets, unlike a 4-ellipse Venn."""
    import numpy as np
    union = sorted(set().union(*sets))
    combos = {}
    for m in union:
        bits = tuple(i for i, s in enumerate(sets) if m in s)
        combos.setdefault(bits, 0)
        combos[bits] += 1
    items = sorted(combos.items(), key=lambda kv: (-kv[1], -len(kv[0])))
    xs = list(range(len(items)))
    sizes = [c for _, c in items]
    ntool = len(names)

    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4.5], height_ratios=[3, 1.6],
                          hspace=0.06, wspace=0.05)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_mat)

    # top: intersection sizes
    ax_bar.bar(xs, sizes, color="#444", width=0.6, zorder=3)
    for x, s in zip(xs, sizes):
        ax_bar.text(x, s, str(s), ha="center", va="bottom", fontsize=10)
    ax_bar.set_ylabel("motifs in intersection")
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.set_ylim(0, max(sizes) * 1.18)
    ax_bar.tick_params(labelbottom=False)
    ax_bar.set_xticks([])

    # bottom: membership matrix (row 0 at top)
    yrow = {i: ntool - 1 - i for i in range(ntool)}   # tool i -> y
    for x, (bits, _) in enumerate(items):
        ax_mat.axvspan(x - 0.5, x + 0.5, color="#f2f2f2" if x % 2 else "white", zorder=0)
        present = set(bits)
        for i in range(ntool):
            y = yrow[i]
            filled = i in present
            ax_mat.scatter(x, y, s=130, color=(colors[i] if filled else "#d9d9d9"), zorder=3)
        if len(present) > 1:
            ys = [yrow[i] for i in present]
            ax_mat.plot([x, x], [min(ys), max(ys)], color="#444", lw=2, zorder=2)
    ax_mat.set_xlim(-0.6, len(items) - 0.4)
    ax_mat.set_ylim(-0.6, ntool - 0.4)
    ax_mat.set_xticks([]); ax_mat.set_yticks([])
    for sp in ax_mat.spines.values():
        sp.set_visible(False)

    # left: per-tool totals (horizontal, pointing left) + tool names on far left
    totals = [len(s) for s in sets]
    ax_set.barh([yrow[i] for i in range(ntool)], totals, color=colors, height=0.5, zorder=3)
    for i in range(ntool):
        ax_set.text(totals[i], yrow[i], f"{totals[i]} ", va="center", ha="left",
                    fontsize=9, color="#333")
    ax_set.invert_xaxis()
    ax_set.set_xlabel("set size (motifs)")
    ax_set.spines[["top", "left", "right"]].set_visible(False)
    ax_set.set_yticks([yrow[i] for i in range(ntool)])
    ax_set.set_yticklabels([names[i] for i in range(ntool)])
    for i, t in enumerate(ax_set.get_yticklabels()):
        t.set_color(colors[i]); t.set_fontweight("bold")
    ax_set.tick_params(length=0)
    # ax_mat shares y with ax_set; suppress the duplicated (clipped) labels on it
    ax_mat.tick_params(labelleft=False, left=False)


# ---------------- analyze ----------------
OI = dict(blue="#0072B2", orange="#E69F00", green="#009E73", vermillion="#D55E00", purple="#CC79A7")

def analyze():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "figure.dpi": 130})
    FIG.mkdir(parents=True, exist_ok=True)

    sets, infos = {}, {}
    sets["MODIFI_hifi"], infos["MODIFI_hifi"] = load_modifi_native()
    for t in MM_TOOLS:
        sets[t], infos[t] = aggregate_motifmaker(t)
    order = ["MODIFI_hifi", "ipdSummary", "fibertools", "jasmine"]
    for nm in order:
        print(f"  {nm:15s} motifs: {len(sets[nm])}  {sorted(sets[nm])[:10]}")

    # motif x tool presence table
    all_motifs = sorted(set().union(*sets.values())) if sets else []
    with open(OUT / "motif_presence.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["motif"] + order + ["modType(any)"])
        for m in all_motifs:
            mt = ""
            for nm in order:
                if m in infos[nm]:
                    mt = infos[nm][m].get("modificationType", "") or mt
            w.writerow([m] + [("1" if m in sets[nm] else "0") for nm in order] + [mt])

    # one combined 4-set Venn (MODIFI-HiFi + all three SOTA tools)
    venn_sets = [sets["MODIFI_hifi"], sets["ipdSummary"], sets["fibertools"], sets["jasmine"]]
    venn_names = ["MODIFI (HiFi)", "ipdSummary", "fibertools (6mA)", "jasmine (5mC)"]
    venn_cols = [OI["blue"], OI["orange"], OI["green"], OI["purple"]]
    # figure source data (per-intersection) next to the figure, for reproducibility
    combo = {}
    for m in all_motifs:
        combo.setdefault(tuple(venn_names[i] for i in range(4) if m in venn_sets[i]), []).append(m)
    with open(FIG / "fig_motif_upset.intersections.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["tools_in_intersection", "n_motifs", "motifs"])
        for key, ms in sorted(combo.items(), key=lambda kv: -len(kv[1])):
            if key:
                w.writerow(["+".join(key), len(ms), ";".join(sorted(ms))])
    import shutil as _sh
    _sh.copyfile(OUT / "motif_presence.csv", FIG / "fig_motif_upset.presence.csv")
    # primary: UpSet plot (clear for 4 sets)
    fig = plt.figure(figsize=(11, 6))
    draw_upset(fig, venn_sets, venn_names, venn_cols)
    fig.suptitle("Methylation-motif concordance across tools (10 circular soil genomes; "
                 "motifs filtered frac≥0.4 & ≥30 sites, RC-collapsed)", fontsize=11)
    fig.savefig(FIG / "fig_motif_upset.png", bbox_inches="tight")
    fig.savefig(FIG / "fig_motif_upset.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG}/fig_motif_upset.png/.pdf "
          f"({len(set().union(*venn_sets))} distinct motifs across all tools)")

    # runtime figure (subread whole-set efficiency)
    ipd, mod = parse_time(T_IPD), parse_time(T_MODIFI_SUB)
    if ipd and mod:
        fig, axs = plt.subplots(1, 3, figsize=(10, 3.6))
        for ax, (k, lab) in zip(axs, [("wall_hr", "Wall-clock (h)"), ("cpu_hr", "CPU (h)"),
                                      ("rss_gb", "Peak RSS (GB)")]):
            vals = [ipd[k], mod[k]]
            ax.bar(["ipdSummary", "MODIFI"], vals, color=[OI["orange"], OI["blue"]], width=0.6, zorder=3)
            for i, v in enumerate(vals):
                ax.text(i, v, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
            ax.set_title(lab); ax.margins(y=0.2); ax.grid(axis="y", alpha=0.25)
            ax.spines[["top", "right"]].set_visible(False)
            ax.set_xlabel(f"{vals[0]/vals[1]:.1f}x lower", fontsize=9, color="#555")
        fig.suptitle("Subread efficiency on test_100 (327 contigs): MODIFI vs ipdSummary", fontsize=12)
        fig.tight_layout()
        fig.savefig(FIG / "fig_runtime_efficiency.png", bbox_inches="tight")
        fig.savefig(FIG / "fig_runtime_efficiency.pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {FIG}/fig_runtime_efficiency.png/.pdf")

    ft_t = sum_times([OUT / "fibertools/selected10.ft_predict.time", OUT / "fibertools/selected10.ft_pileup.time"])
    jas_t = sum_times([OUT / "jasmine/selected10.jasmine.time", OUT / "jasmine/selected10.cpgtools.time"])
    mh_t = sum_times([OUT / "modifi_hifi/selected10.modifi_hifi.time"])

    M = sets["MODIFI_hifi"]
    def line(b):
        B = sets[b]; sh = len(M & B); tot = len(M | B)
        return (f"- **MODIFI-HiFi vs {b}**: MODIFI={len(M)}, {b}={len(B)}, shared={sh}, "
                f"Jaccard={sh/tot:.2f}" + (f"  (shared: {sorted(M & B)})" if sh else ""))
    md = ["# MODIFI vs SOTA — motif-level comparison (10 circular soil genomes, test_100)\n",
          "Motifs per tool = union over contigs, filtered by MODIFI defaults (fraction≥0.4, nDetected≥30), reverse-complement pairs collapsed to one.\n",
          "Per-genome motif discovery; union of motifs across the 10 genomes. MODIFI uses its "
          "native motifs; ipdSummary/fibertools/jasmine use motifMaker on their per-contig calls.\n",
          "## Overlap with MODIFI-HiFi\n", line("ipdSummary"), line("fibertools"), line("jasmine"), ""]
    md.append("## Runtime / memory\n")
    md.append("| tool | read type | scope | wall (h) | CPU (h) | peak RSS (GB) |")
    md.append("|---|---|---|---|---|---|")
    if ipd: md.append(f"| ipdSummary | subreads | test_100 (327 ctg) | {ipd['wall_hr']:.2f} | {ipd['cpu_hr']:.1f} | {ipd['rss_gb']:.1f} |")
    if mod: md.append(f"| MODIFI | subreads | test_100 (327 ctg) | {mod['wall_hr']:.2f} | {mod['cpu_hr']:.1f} | {mod['rss_gb']:.1f} |")
    if mh_t["n"]: md.append(f"| MODIFI | HiFi | merged 10 genomes | {mh_t['wall_hr']:.2f} | {mh_t['cpu_hr']:.2f} | {mh_t['rss_gb']:.1f} |")
    if ft_t["n"]: md.append(f"| fibertools | HiFi | merged 10 genomes | {ft_t['wall_hr']:.2f} | {ft_t['cpu_hr']:.2f} | {ft_t['rss_gb']:.1f} |")
    if jas_t["n"]: md.append(f"| jasmine+pb-CpG | HiFi | merged 10 genomes | {jas_t['wall_hr']:.2f} | {jas_t['cpu_hr']:.2f} | {jas_t['rss_gb']:.1f} |")
    if ipd and mod:
        md.append(f"\n**Subread efficiency:** MODIFI is {ipd['wall_hr']/mod['wall_hr']:.1f}x faster wall-clock, "
                  f"{ipd['cpu_hr']/mod['cpu_hr']:.1f}x less CPU, {ipd['rss_gb']/mod['rss_gb']:.1f}x less peak memory than ipdSummary.\n")
    (OUT / "summary.md").write_text("\n".join(md))
    (FIG / "summary.md").write_text("\n".join(md))   # next to the figure, for reproduction
    print(f"  wrote {OUT}/summary.md and {OUT}/motif_presence.csv")
    print("\n".join(md))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "analyze"
    if mode == "gffs":
        build_gffs()
    elif mode == "analyze":
        analyze()
    else:
        sys.exit("usage: motif_comparison.py [gffs|analyze]")
