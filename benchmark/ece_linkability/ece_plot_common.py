#!/usr/bin/env python3
"""Shared helpers for the ECE linkage-ability figures (isolate + simulated)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

MASTER = "/home/shuaiw/borg/revision/ece_linkability/ece_linkability_master.csv"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_linkability"

plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3,
                     "axes.axisbelow": True, "figure.dpi": 110})
# Wong / Okabe-Ito colorblind-safe palette
BLUE, ORANGE, GREEN, GREY, PURPLE = "#0072B2", "#D55E00", "#009E73", "#888888", "#8256A5"
C_LINK, C_UNLINK, C_HOST = GREEN, ORANGE, BLUE


def load_master():
    return pd.read_csv(MASTER, low_memory=False)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, 0.0, 0.0
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return p, max(0.0, c - half), min(1.0, c + half)


def mci(v):
    v = np.array([x for x in v if x == x], float)
    if not len(v):
        return np.nan, 0.0
    return v.mean(), (1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def ecdf(x):
    x = np.sort(np.asarray([v for v in x if v == v], float))
    if not len(x):
        return np.array([]), np.array([])
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def mwu_p(a, b):
    a = [v for v in a if v == v]
    b = [v for v in b if v == v]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    try:
        return stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
    except ValueError:
        return np.nan


def p_star(p):
    if p != p:
        return "n/a"
    if p < 1e-3:
        return f"p={p:.1e}"
    return f"p={p:.3f}"


def quantile_bins(values, nbins=5):
    """Return bin edges by quantile over finite values (unique)."""
    v = np.asarray([x for x in values if x == x], float)
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbins + 1)))
    return edges


def binned_rate(df, xcol, ycol, edges):
    """Linked fraction + Wilson CI + bin center per quantile bin."""
    rows = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i == len(edges) - 2:
            m = (df[xcol] >= lo) & (df[xcol] <= hi)
        else:
            m = (df[xcol] >= lo) & (df[xcol] < hi)
        sub = df[m]
        n = len(sub)
        k = int(sub[ycol].sum())
        p, clo, chi = wilson(k, n)
        center = float(np.median(sub[xcol])) if n else (lo + hi) / 2
        rows.append(dict(bin_lo=lo, bin_hi=hi, center=center, n=n, k=k,
                         rate=p, ci_lo=clo, ci_hi=chi))
    return pd.DataFrame(rows)


def save(fig, stem):
    os.makedirs(OUT, exist_ok=True)
    pdf = os.path.join(OUT, stem + ".pdf")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {pdf} (+ .png)")
