"""
Cross-species reproducibility of per-10-mer control IPD.

Tests whether the same 10-mer context has the same (Schadt-normalized) IPD in different
genomes. Using the three WGA (modification-free) strains, we compute the mean normalized
IPD per MODIFI 10-mer (u=7 upstream, d=2 downstream; k=10) in each strain and, for each
species pair, plot standardized per-10-mer IPD of one species vs the other (Schadt-style)
and report the correlation.

Reuses the k-mer machinery from one_way_anova.py so the context/strand convention is
identical to the variance-decomposition analysis.

Run:  conda run -n modifi python benchmark/check_context/cross_genome_ipd.py
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, linregress
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import one_way_anova as owa

U, D = 7, 2                 # MODIFI 10-mer: 7 upstream, 2 downstream, k=10
MIN_POS = 3                 # min positions per 10-mer required in BOTH strains of a pair
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/cross_species_ipd"

REVDIR = "/home/shuaiw/borg/revision/context"
REFDIR = "/home/shuaiw/methylation/data/published_data/fanggang/ref"
STRAINS = [
    dict(key="C227", label="E. coli C227",
         cache=f"{REVDIR}/c224/c227_ipd_schadt.csv", ref=f"{REFDIR}/C227.fa"),
    dict(key="DSM", label="C. israelensis DSM 3043",
         cache=f"{REVDIR}/dsm/dsm_ipd_schadt.csv", ref=f"{REFDIR}/DSM-3043.fa"),
    dict(key="J99", label="H. pylori J99",
         cache=f"{REVDIR}/j99/j99_ipd_schadt.csv", ref=f"{REFDIR}/J99.fa"),
]
PAIRS = [("C227", "DSM"), ("C227", "J99"), ("DSM", "J99")]


def per_kmer_table(strain):
    """Return DataFrame: 10-mer code -> (mean IPD, n positions) for one strain."""
    _, ref, comp = owa.load_reference(strain["ref"])
    df = pd.read_csv(strain["cache"])
    fwd = df["strand"].values == "+"
    pos_f = df.loc[fwd, "tpl"].values.astype(np.int64)
    y_f = df.loc[fwd, "y"].values.astype(np.float64)
    pos_r = df.loc[~fwd, "tpl"].values.astype(np.int64)
    y_r = df.loc[~fwd, "y"].values.astype(np.float64)

    pf, yf, pr, yr = owa.valid_sites(ref, comp, pos_f, y_f, pos_r, y_r, U, D)
    codes = owa.codes_at(ref, comp, pf, pr, U, D)          # pooled incorporating-strand
    y = np.concatenate([yf, yr])

    uniq, inv = np.unique(codes, return_inverse=True)
    n = np.bincount(inv, minlength=uniq.size)
    s = np.bincount(inv, weights=y, minlength=uniq.size)
    mean = s / n
    out = pd.DataFrame({"code": uniq, f"mean_{strain['key']}": mean,
                        f"n_{strain['key']}": n})
    print(f"  {strain['key']}: {uniq.size} distinct 10-mers "
          f"({(n >= MIN_POS).sum()} with >= {MIN_POS} positions)", flush=True)
    return out


def zscore(x):
    return (x - x.mean()) / x.std(ddof=0)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("aggregating per-10-mer IPD per strain ...", flush=True)
    tabs = {s["key"]: per_kmer_table(s) for s in STRAINS}
    labels = {s["key"]: s["label"] for s in STRAINS}

    summary = []
    merged = {}
    for a, b in PAIRS:
        m = tabs[a].merge(tabs[b], on="code", how="inner")
        m = m[(m[f"n_{a}"] >= MIN_POS) & (m[f"n_{b}"] >= MIN_POS)].copy()
        m[f"z_{a}"] = zscore(m[f"mean_{a}"].values)
        m[f"z_{b}"] = zscore(m[f"mean_{b}"].values)
        r, _ = pearsonr(m[f"z_{a}"], m[f"z_{b}"])
        rho, _ = spearmanr(m[f"z_{a}"], m[f"z_{b}"])
        sl = linregress(m[f"z_{a}"], m[f"z_{b}"]).slope
        summary.append(dict(pair=f"{a}-{b}", n=len(m), pearson_r=r,
                            spearman_rho=rho, slope=sl))
        merged[(a, b)] = m
        m.to_csv(os.path.join(OUTDIR, f"cross_species_{a}_{b}.csv"), index=False)
        print(f"{a}-{b}: n={len(m):>6d}  Pearson r={r:.3f}  Spearman rho={rho:.3f}  "
              f"slope={sl:.3f}", flush=True)

    pd.DataFrame(summary).to_csv(
        os.path.join(OUTDIR, "cross_species_correlation_summary.csv"), index=False)
    plot(merged, labels, summary)
    print("wrote figure + tables to", OUTDIR, flush=True)


def _italic(name):
    # mathtext italic species label, e.g. "$\it{E.\ coli\ C227}$"
    esc = name.replace(" ", r"\ ")
    return "Standardized IPD in $\\it{" + esc + "}$"


def plot(merged, labels, summary):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))
    for ax, (a, b), srow in zip(axes, PAIRS, summary):
        m = merged[(a, b)]
        x, y = m[f"z_{a}"].values, m[f"z_{b}"].values
        ax.scatter(x, y, s=4, c="black", alpha=0.35, edgecolors="none")
        lim = [min(x.min(), y.min()) - 0.2, max(x.max(), y.max()) + 0.2]
        ax.plot(lim, lim, color="#333333", lw=1)              # identity line y = x
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(_italic(labels[a]), fontsize=9)
        ax.set_ylabel(_italic(labels[b]), fontsize=9)
        ax.text(0.05, 0.93, f"r = {srow['pearson_r']:.2f}\nn = {srow['n']:,}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9)
        ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "cross_species_ipd.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
