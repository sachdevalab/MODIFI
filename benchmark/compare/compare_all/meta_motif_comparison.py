#!/usr/bin/env python3
"""
Whole-metagenome (v2) motif-level comparison over the selected contigs.
Per-genome motif discovery already done; here we take the UNION per tool across
the selected contigs and compare, centred on MODIFI. Reuses draw_upset / parse_time
from motif_comparison.py.

MODIFI  = native de-novo motifs (modifi_hifi_full/motifs/<c>.motifs.csv)
ipdSummary / fibertools / jasmine = per-contig motifMaker outputs (percontig/<tool>/<c>.motifs.csv)
Tools with no outputs yet (e.g. ipdSummary before the subread run finishes) are skipped.

Outputs to compare_all_meta/ + tmp/rev_figs/compare_all_meta/.
"""
from __future__ import annotations
import csv, glob, os, sys, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
from motif_comparison import draw_upset, parse_time, OI  # reuse
from Bio.Seq import Seq  # same RC dependency the repo's derep_motifs.py uses

# MODIFI default motif-retention criteria (main.py defaults)
MIN_FRAC = 0.4
MIN_SITES = 30


def canon(m, cpos):
    """Canonical (motifString, centerPos) identity. A motif and its reverse
    complement are the same modification and are collapsed to one entry; when
    flipping to the RC the methylated position is flipped too
    (centerPos -> len-1-centerPos), so the same motif keeps the same centerPos.
    Two motifs are 'the same' only if BOTH string and centerPos match."""
    m = m.strip().upper()
    try:
        rc = str(Seq(m).reverse_complement())
    except Exception:
        rc = m
    rc_cpos = len(m) - 1 - cpos
    return min((m, cpos), (rc, rc_cpos))

OUT = "/home/shuaiw/borg/paper/ipdsummary/compare_all_meta"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/compare_all_meta"
MODIFI_HIFI = f"{OUT}/modifi_hifi_full"
MODIFI_SUB = f"{OUT}/modifi_subread_full"
PC = f"{OUT}/motifs_compare/percontig"
TIMES = f"{OUT}/times"
SEL = f"{OUT}/selected_contigs.tsv"


def selected_contigs():
    with open(SEL) as f:
        return [r["contig"] for r in csv.DictReader(f, delimiter="\t")]


def union_motifs(files):
    """Union of motifs across contigs, filtered by MODIFI defaults
    (fraction >= 0.4 AND nDetected >= 30) and RC-collapsed to canonical form."""
    s = set()
    for p in files:
        if not os.path.isfile(p):
            continue
        for row in csv.DictReader(open(p)):
            m = (row.get("motifString") or "").strip().upper()
            if not m:
                continue
            try:
                frac = float(row.get("fraction", 0) or 0)
                nd = int(float(row.get("nDetected", 0) or 0))
                cpos = int(float(row.get("centerPos", 0) or 0))
            except ValueError:
                continue
            if frac < MIN_FRAC or nd < MIN_SITES:
                continue
            s.add(canon(m, cpos))          # (motifString, centerPos) identity
    return s


def sum_tool_times(contigs, suffixes):
    cpu = wall = rss = 0.0
    n = 0
    for c in contigs:
        c_cpu = c_wall = c_rss = 0.0
        got = False
        for suf in suffixes:
            r = parse_time(f"{TIMES}/{c}.{suf}.time")
            if r:
                c_cpu += r["cpu_hr"]; c_wall += r["wall_hr"]; c_rss = max(c_rss, r["rss_gb"]); got = True
        if got:
            cpu += c_cpu; wall += c_wall; rss = max(rss, c_rss); n += 1
    return dict(cpu_hr=cpu, wall_hr=wall, rss_gb=rss, n=n)


def main():
    os.makedirs(FIG, exist_ok=True)
    contigs = selected_contigs()
    print(f"selected contigs: {len(contigs)}")

    MODIFI_SUB216 = f"{OUT}/modifi_subread_216"
    sets = {}
    sets["MODIFI-HiFi"] = union_motifs([f"{MODIFI_HIFI}/motifs/{c}.motifs.csv" for c in contigs])
    sets["MODIFI-sub"] = union_motifs([f"{MODIFI_SUB216}/motifs/{c}.motifs.csv" for c in contigs])
    sets["ipdSummary"] = union_motifs([f"{PC}/ipdSummary/{c}.motifs.csv" for c in contigs])
    sets["fibertools"] = union_motifs([f"{PC}/fibertools/{c}.motifs.csv" for c in contigs])
    sets["SMAC"] = union_motifs([f"{PC}/SMAC/{c}.motifs.csv" for c in contigs])   # single-molecule 6mA

    # include a tool if it was RUN (has per-contig motif files), even if 0 motifs survive.
    def ran(t):
        d = {"MODIFI-HiFi": f"{MODIFI_HIFI}/motifs", "MODIFI-sub": f"{MODIFI_SUB216}/motifs"}.get(t, f"{PC}/{t}")
        return len(glob.glob(f"{d}/*.motifs.csv")) > 0
    # jasmine excluded: it is a mammalian 5mC-CpG caller, not designed for microbes.
    # SMAC appears once its per-contig motifs exist.
    order = [t for t in ["MODIFI-HiFi", "MODIFI-sub", "ipdSummary", "SMAC", "fibertools"]
             if sets[t] or ran(t)]
    cols = {"MODIFI-HiFi": OI["blue"], "MODIFI-sub": OI["vermillion"], "ipdSummary": OI["orange"],
            "SMAC": OI["purple"], "fibertools": OI["green"]}
    for t in order:
        print(f"  {t:12s}: {len(sets[t])} motifs")

    # motif x tool presence table + the exact per-intersection data behind the UpSet,
    # written NEXT TO THE FIGURE for reproducibility (and to OUT).
    fmt = lambda mp: f"{mp[0]}@{mp[1]}"   # (motifString, centerPos) -> "MOTIF@pos"
    allm = sorted(set().union(*[sets[t] for t in order]))
    for d in (OUT, FIG):
        with open(f"{d}/fig_motif_upset.presence.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["motif@centerPos"] + order)
            for m in allm:
                w.writerow([fmt(m)] + [("1" if m in sets[t] else "0") for t in order])
    combo = {}
    for m in allm:
        combo.setdefault(tuple(t for t in order if m in sets[t]), []).append(m)
    for d in (OUT, FIG):
        with open(f"{d}/fig_motif_upset.intersections.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["tools_in_intersection", "n_motifs", "motifs"])
            for key, ms in sorted(combo.items(), key=lambda kv: -len(kv[1])):
                w.writerow(["+".join(key), len(ms), ";".join(fmt(x) for x in sorted(ms))])

    # UpSet
    fig = plt.figure(figsize=(12, 6))
    draw_upset(fig, [sets[t] for t in order], order, [cols[t] for t in order])
    tag = "" if "ipdSummary" in order else " (HiFi side; MODIFI-sub + ipdSummary pending)"
    fig.suptitle(f"Methylation-motif concordance across tools — {len(contigs)} contigs "
                 f"(≥100kb, HiFi≥10x); motifs filtered frac≥0.4 & ≥30 sites, RC-collapsed{tag}", fontsize=11)
    fig.savefig(f"{FIG}/fig_motif_upset.png", bbox_inches="tight")
    fig.savefig(f"{FIG}/fig_motif_upset.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG}/fig_motif_upset.png + source data (presence + intersections)")

    # resource table
    ft = sum_tool_times(contigs, ["ft_predict", "ft_pileup"])
    jas = sum_tool_times(contigs, ["jasmine", "cpgtools"])
    ipd = sum_tool_times(contigs, ["ipdSummary"])
    mh = parse_time(f"{MODIFI_HIFI}.time")
    ms = parse_time(f"{MODIFI_SUB216}.time")   # MODIFI subreads on the 216-contig sub-ref

    md = [f"# Whole-metagenome motif comparison — {len(contigs)} contigs (≥100kb, HiFi≥10x)\n",
          "Motifs per tool = union over contigs, filtered by MODIFI defaults (fraction≥0.4, nDetected≥30); a motif is identified by (motifString, centerPos), and reverse-complement pairs (with matched centerPos) collapsed to one.\n",
          "## Motif counts and overlap\n"]
    if sets.get("MODIFI-HiFi"):
        M = sets["MODIFI-HiFi"]
        for t in order:
            if t.startswith("MODIFI"):
                continue
            sh = len(M & sets[t]); tot = len(M | sets[t])
            md.append(f"- **MODIFI-HiFi vs {t}**: MODIFI-HiFi={len(M)}, {t}={len(sets[t])}, shared={sh}, Jaccard={sh/tot:.2f}")
    if sets.get("MODIFI-sub") and sets.get("ipdSummary"):   # key same-read-type comparison
        a, b = sets["MODIFI-sub"], sets["ipdSummary"]; sh = len(a & b)
        md.append(f"- **MODIFI-sub vs ipdSummary (both subreads)**: MODIFI-sub={len(a)}, ipdSummary={len(b)}, "
                  f"shared={sh}, Jaccard={sh/len(a | b):.2f}")
    md += ["", "## Running time / memory\n",
           "| tool | scope | wall | CPU-h | peak RSS (GB) |", "|---|---|---|---|---|"]
    if mh:
        md.append(f"| MODIFI (HiFi) | whole metagenome, 1 run | {mh['wall_hr']:.2f} h | {mh['cpu_hr']:.0f} | {mh['rss_gb']:.0f} |")
    if ms:
        md.append(f"| MODIFI (subread) | 216 contigs, 1 run | {ms['wall_hr']:.2f} h | {ms['cpu_hr']:.1f} | {ms['rss_gb']:.1f} |")
    if ipd["n"]:
        md.append(f"| ipdSummary | {ipd['n']} contigs (per-contig) | {ipd['wall_hr']:.2f} h | {ipd['cpu_hr']:.1f} | {ipd['rss_gb']:.1f} |")
    md.append(f"| fibertools | {ft['n']} contigs (per-contig) | {ft['wall_hr']:.2f} h | {ft['cpu_hr']:.1f} | {ft['rss_gb']:.1f} |")
    md.append("\n_Per-contig wall/contig_: "
              f"fibertools {ft['wall_hr']*3600/max(ft['n'],1):.0f}s"
              + (f", ipdSummary {ipd['wall_hr']*3600/max(ipd['n'],1):.0f}s" if ipd['n'] else "")
              + ". (fibertools CPU-h inflated by node oversubscription; wall is reliable.)\n")
    for d in (OUT, FIG):
        open(f"{d}/summary.md", "w").write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
