"""
Flanking-composition test within fixed 10-mers (reviewer mechanism).

The k-mer R^2 plateau implies that composition beyond the 10-mer is "absorbed". We test
this directly: for each well-sampled 10-mer in the modification-free WGA data, bin its
genomic instances by the GC of the surrounding +-25 bp (and +-100 bp), take the mean
(Schadt-normalized) log-IPD per bin, and record the WITHIN-10-mer spread (max-min of the
bin means) -- i.e. how much a *fixed* 10-mer's IPD moves with flanking GC. A GC-label
permutation gives the pure-sampling-noise null; observed - null is the real flanking
effect. All in log-IPD units, directly comparable to the methylation signal
(log of the IPD ratio; a modification is ~ln(2..6) = 0.7..1.8).

Run:  conda run -n modifi python benchmark/check_context/gc_bin_analysis.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import one_way_anova as owa

U, D = 7, 2
MININST = 10           # min genomic instances of a 10-mer (so its own mean is stable)
WINDOWS = [25, 50, 100, 250, 500, 1000, 2000]   # flanking half-widths (bp)
GC_EDGES = np.arange(0.0, 1.0001, 0.02)          # flanking-GC bins for variance decomp
COLORS = {"E. coli C227": "#3b7dd8", "C. israelensis DSM 3043": "#2ca25f",
          "H. pylori J99": "#d1495b"}

REVDIR = "/home/shuaiw/borg/revision/context"
REFDIR = "/home/shuaiw/methylation/data/published_data/fanggang/ref"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/gc_bin"
STRAINS = [
    ("E. coli C227", f"{REVDIR}/c224/c227_ipd_schadt.csv", f"{REFDIR}/C227.fa"),
    ("C. israelensis DSM 3043", f"{REVDIR}/dsm/dsm_ipd_schadt.csv", f"{REFDIR}/DSM-3043.fa"),
    ("H. pylori J99", f"{REVDIR}/j99/j99_ipd_schadt.csv", f"{REFDIR}/J99.fa"),
]


def flanking_gc(ref, pos, half):
    """GC fraction over [p-half, p+half] for each position p (vectorized)."""
    isgc = ((ref == 2) | (ref == 3)).astype(np.int64)      # C=2, G=3
    cg = np.concatenate([[0], np.cumsum(isgc)])
    L = ref.shape[0]
    lo = np.clip(pos - half, 0, L)
    hi = np.clip(pos + half + 1, 0, L)
    return (cg[hi] - cg[lo]) / (hi - lo)


def analyse(name, cache, ref_path):
    _, ref, comp = owa.load_reference(ref_path)
    df = pd.read_csv(cache)
    fwd = df["strand"].values == "+"
    pf = df["tpl"].values[fwd].astype(np.int64)
    yf = df["y"].values[fwd]
    pr = df["tpl"].values[~fwd].astype(np.int64)
    yr = df["y"].values[~fwd]
    pf, yf, pr, yr = owa.valid_sites(ref, comp, pf, yf, pr, yr, U, D)
    codes = owa.codes_at(ref, comp, pf, pr, U, D)
    y = np.concatenate([yf, yr])
    pos = np.concatenate([pf, pr])

    # residual = position IPD minus its 10-mer's own mean (only 10-mers with >=MININST
    # instances, so the mean is stable); isolates WITHIN-10-mer variation.
    uniq, inv, cnt = np.unique(codes, return_inverse=True, return_counts=True)
    mean_code = np.bincount(inv, weights=y) / cnt
    resid = (y - mean_code[inv])[cnt[inv] >= MININST]
    keep = cnt[inv] >= MININST
    total = resid.var()
    grand = resid.mean()

    eta = {}
    for w in WINDOWS:
        gc = flanking_gc(ref, pos, w)[keep]
        b = np.clip(np.digitize(gc, GC_EDGES) - 1, 0, GC_EDGES.size - 2)
        nb = np.bincount(b, minlength=GC_EDGES.size - 1)
        sb = np.bincount(b, weights=resid, minlength=GC_EDGES.size - 1)
        ok = nb > 0
        mb = np.zeros_like(sb); mb[ok] = sb[ok] / nb[ok]
        between = np.sum(nb * (mb - grand) ** 2) / resid.size   # cross-GC variance
        eta[w] = between / total                                # fraction explained by GC
        print(f"  {name} ±{w}bp: eta2_GC={eta[w]*100:.2f}%  "
              f"(cross-GC SD={np.sqrt(between):.3f}, within SD={np.sqrt(total-between):.3f})",
              flush=True)
    return eta


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    etas, summ = {}, []
    for name, cache, ref in STRAINS:
        print(name, flush=True)
        etas[name] = analyse(name, cache, ref)
        for w in WINDOWS:
            summ.append(dict(strain=name, flank_bp=w, eta2_gc_pct=round(etas[name][w] * 100, 3)))
    pd.DataFrame(summ).to_csv(os.path.join(OUTDIR, "gc_bin_summary.csv"), index=False)
    print("\n", pd.DataFrame(summ).to_string(index=False), flush=True)
    plot_eta_vs_window(etas)
    print("\nwrote", OUTDIR, flush=True)


def plot_eta_vs_window(etas):
    """Fraction of within-10-mer IPD variation explained by flanking GC vs flanking
    window size. Peaks ~3.7% at +-100 bp then declines -> GC is negligible at every
    scale, so the k-mer baseline is not confounded by flanking composition."""
    fig, ax = plt.subplots(figsize=(9, 6))
    for name, _, _ in STRAINS:
        v = [etas[name][w] * 100 for w in WINDOWS]
        ax.plot(WINDOWS, v, "-o", color=COLORS[name], lw=2, ms=6, label=name)
    ax.set_xscale("log")
    ax.set_xticks(WINDOWS)
    ax.set_xticklabels([f"±{w}" for w in WINDOWS])
    ax.set_ylim(0, 100)
    ax.axhline(100, color="#999", lw=1, ls="--")
    ax.text(WINDOWS[-1], 100, "100% = all within-10-mer IPD variation  ",
            va="bottom", ha="right", fontsize=10, color="#666")
    ax.set_xlabel("flanking window half-width (bp)", fontsize=13)
    ax.set_ylabel("within-10-mer IPD variation explained by flanking GC  (η², %)",
                  fontsize=12)
    ax.tick_params(labelsize=12)
    ax.minorticks_off()
    ax.legend(fontsize=11, frameon=False, loc="upper left")
    ax.grid(ls=":", alpha=0.5)
    # zoom inset so the (tiny) shape is legible while the main axis shows the 0-100 scale
    axin = ax.inset_axes([0.5, 0.42, 0.46, 0.45])
    for name, _, _ in STRAINS:
        v = [etas[name][w] * 100 for w in WINDOWS]
        axin.plot(WINDOWS, v, "-o", color=COLORS[name], lw=2, ms=5)
    axin.set_xscale("log"); axin.set_xticks(WINDOWS)
    axin.set_xticklabels([str(w) for w in WINDOWS], fontsize=8)
    axin.tick_params(labelsize=8); axin.minorticks_off()
    axin.set_title("zoom", fontsize=9)
    axin.grid(ls=":", alpha=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "gc_baseline_robustness.pdf"))
    plt.close(fig)


def plot_variance(decomps):
    """Bar comparison: within-10-mer variation vs GC-driven variation (SD, log-IPD)."""
    names = [s[0] for s in STRAINS]
    fig, axes = plt.subplots(1, len(WINDOWS), figsize=(12, 5.2), sharey=True)
    x = np.arange(len(names))
    for ax, w in zip(axes, WINDOWS):
        within = [decomps[n][w]["within_bin_sd"] for n in names]
        cross = [decomps[n][w]["cross_gc_sd"] for n in names]
        ax.bar(x - 0.2, within, 0.38, color="#7b9acc",
               label="within-GC-bin (10-mer internal)")
        ax.bar(x + 0.2, cross, 0.38, color="#e08214",
               label="cross-GC (GC-driven)")
        for xi, n in zip(x, names):
            e = decomps[n][w]["eta2_gc"]
            ax.text(xi + 0.2, decomps[n][w]["cross_gc_sd"] + 0.005,
                    f"η²={e*100:.1f}%", ha="center", va="bottom", fontsize=9, color="#8a4b0a")
        ax.set_xticks(x)
        ax.set_xticklabels([n.replace(" ", "\n", 1) for n in names], fontsize=11)
        ax.set_title(f"flanking window ±{w} bp", fontsize=13)
        ax.set_ylabel("within-10-mer IPD variation (SD, log-IPD)", fontsize=12)
        ax.tick_params(labelsize=11)
        ax.grid(axis="y", ls=":", alpha=0.5)
        ax.legend(fontsize=10, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "gc_variance_decomposition.pdf"))
    plt.close(fig)


def plot(results):
    names = [s[0] for s in STRAINS]
    fig, axes = plt.subplots(1, len(WINDOWS), figsize=(13, 5.2), sharey=True)
    for ax, w in zip(axes, WINDOWS):
        ax.axhline(0, color="#888", lw=1, ls="--", zorder=1)
        for name in names:
            c, m, se, n = results[name][w]
            if c.size == 0:
                continue
            col = COLORS[name]
            ax.fill_between(c, m - 1.96 * se, m + 1.96 * se, color=col, alpha=0.2, lw=0)
            ax.plot(c, m, "-o", color=col, lw=2, ms=4, label=name)
        ax.set_title(f"flanking window ±{w} bp", fontsize=13)
        ax.set_xlabel("flanking GC content", fontsize=13)
        ax.set_ylabel("within-10-mer IPD deviation (log-IPD)\n"
                      "position IPD − its 10-mer mean", fontsize=12)
        ax.tick_params(labelsize=11)
        ax.grid(ls=":", alpha=0.5)
        ax.legend(fontsize=10, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "gc_bin_spread.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
