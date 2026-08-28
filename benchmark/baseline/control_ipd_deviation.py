"""
Note S4 | Analysis of the control IPD inference.

A reviewer challenged MODIFI's assumption that most instances of a k-mer are
unmodified, which may not hold for low-diversity microbiomes, dominant organisms
with heavily modified motifs, or communities enriched for specific RM systems.

The control IPD deviation is
    sigma = a / u = 1 + phi * (r - 1),
where phi = Nmod / Ntotal is the fraction of all instances of a k-mer that are
modified in the assembly, and r = h / u is the modification effect ratio (the
factor by which the modified IPD deviates from the unmodified control). For a
single dominant genome of length t in a metagenome of total length X the fraction
reduces to phi = t / X, giving sigma = 1 + (t / X) * (r - 1).

An inflated baseline (sigma > 1) propagates to the estimated ipdRatio of the
modified motif: MODIFI's estimated control is a = u * sigma, so
    estimated ipdRatio = h / a = (h / u) / sigma = r / sigma,
which is deflated from r toward the detection floor of 1. This sensitivity loss
is confined to the shared-motif k-mers themselves; all other k-mers keep sigma = 1
and recover the true ipdRatio r.

This script produces four panels:
  A  sigma vs metagenome length X   (recap of Supplementary Figure S15)
  B  sigma vs dominant genome length t
  C  sigma vs modified k-mer fraction phi
  D  estimated ipdRatio (r / sigma) vs sigma

Pure numpy/matplotlib math; safe to run on the login node.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/baseline"

# Okabe-Ito colorblind-safe palette, one color per modification effect ratio.
R_VALUES = [1.5, 2, 2.5, 3]
COLORS = {1.5: "#0072B2", 2: "#D55E00", 2.5: "#009E73", 3: "#CC79A7"}

# Fixed points used across panels so the panels agree at the shared 0.05 fraction.
T_FIXED = 5.0     # representative dominant-genome length (Mbp), panel A
X_FIXED = 100.0   # representative metagenome length (Mbp), panel B

# Panel B: dominance experiment. A fixed community of N_GENOMES equal-length
# genomes is sequenced with a fixed total yield; relative abundances follow a
# log-normal distribution whose spread LN_SIGMA sets the dominance (small = even,
# large = one/few genomes dominate). Each genome's depth = abundance * TOTAL_SEQ /
# GENOME_LEN; a genome is covered when depth >= DEPTH_FLOOR. For a k-mer modified
# in a single genome among the n covered (equal-length) genomes, the pooled
# modified fraction is phi = 1/n, so the dominant motif's estimated IPD ratio is
# r / (1 + phi(r-1)) = r / (1 + (r-1)/n).
N_GENOMES = 100              # community size, equal-length genomes
GENOME_LEN = 5.0             # Mbp, each genome
TOTAL_SEQ = 10_000.0         # Mbp of sequencing (10 Gbp), fixed budget
MEAN_DEPTH = TOTAL_SEQ / (N_GENOMES * GENOME_LEN)   # 20x
DEPTH_FLOOR = 5.0            # MODIFI calling floor (min_ctg_cov, main.py:126)
LN_SIGMA_GRID = np.linspace(0.3, 6.0, 200)   # dominance axis (log-normal sigma)
N_REP = 2000                 # replicate communities per sigma
RNG_SEED = 42


def sigma_single(t, X, r):
    """Control IPD deviation for a single dominant genome: 1 + (t/X)(r-1)."""
    return 1 + (t / X) * (r - 1)


def sigma_phi(phi, r):
    """Control IPD deviation from the pooled modified fraction: 1 + phi(r-1)."""
    return 1 + phi * (r - 1)


def _legend_label(r):
    return f"modification effect ratio r = {r}"


def _save(fig, basename):
    png = os.path.join(OUTDIR, basename + ".png")
    pdf = os.path.join(OUTDIR, basename + ".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"  wrote {png}")
    print(f"  wrote {pdf}")


def _save_sourcedata(rows, cols, basename):
    df = pd.DataFrame(rows, columns=cols)
    path = os.path.join(OUTDIR, basename + "_sourcedata.csv")
    df.to_csv(path, index=False)
    print(f"  wrote {path}")
    return df


# --------------------------------------------------------------------------- #
# Panel A: sigma vs number of genomes in the metagenome. With one modified genome
# among n genomes of equal length, phi = t/X = 1/n, so sigma = 1 + (r-1)/n.
# --------------------------------------------------------------------------- #
def plot_panel_A(ax):
    n = np.arange(1, 601, dtype=float)   # number of genomes in the metagenome
    rows = []
    for r in R_VALUES:
        y = 1 + (r - 1) / n              # sigma = 1 + (t/X)(r-1) with t/X = 1/n
        ax.plot(n, y, linewidth=2, color=COLORS[r], label=_legend_label(r))
        rows.extend([(ni, r, yi) for ni, yi in zip(n, y)])
    ax.axhline(1, ls="--", color="#888888", lw=1)
    ax.set_xlabel("Number of genomes in metagenome")
    ax.set_ylabel("Control IPD deviation")
    ax.set_ylim(bottom=1)
    ax.set_title("a", loc="left", fontweight="bold")
    return _save_sourcedata(rows, ["n_genomes", "r", "sigma"],
                            "sigma_vs_metagenome_length")


# --------------------------------------------------------------------------- #
# Panel B: dominance experiment (estimated IPD ratio vs log-normal sigma)
# --------------------------------------------------------------------------- #
# N_GENOMES equal-length genomes are sequenced with a fixed total yield TOTAL_SEQ.
# Relative abundances are log-normal with spread LN_SIGMA (the dominance axis).
# Each genome's depth = abundance * TOTAL_SEQ / GENOME_LEN; n = # genomes with
# depth >= DEPTH_FLOOR are covered. For a k-mer modified in a single genome among
# the n covered (equal-length) genomes, phi = 1/n, so the dominant motif's
# estimated IPD ratio is r / (1 + phi(r-1)) = r / (1 + (r-1)/n). Larger sigma means
# stronger dominance -> smaller n -> larger phi -> lower estimated IPD ratio.

def compute_dominance_experiment(sigma_grid=LN_SIGMA_GRID):
    """Over N_REP log-normal communities per sigma, return the covered-genome count
    n (mean, 5th, 95th pct) and the mean estimated IPD ratio per r (phi = 1/n)."""
    rng = np.random.default_rng(RNG_SEED)
    S = sigma_grid.size
    n_mean = np.empty(S)
    n_p5 = np.empty(S)
    n_p95 = np.empty(S)
    ipd_mean = {r: np.empty(S) for r in R_VALUES}
    for k, sig in enumerate(sigma_grid):
        x = rng.lognormal(mean=0.0, sigma=sig, size=(N_REP, N_GENOMES))
        a = x / x.sum(axis=1, keepdims=True)
        depth = a * TOTAL_SEQ / GENOME_LEN               # (N_REP, N_GENOMES)
        n = (depth >= DEPTH_FLOOR).sum(axis=1)           # covered genomes per rep
        n_mean[k] = n.mean()
        n_p5[k], n_p95[k] = np.percentile(n, [5, 95])
        phi = 1.0 / n                                    # one modified genome of n
        for r in R_VALUES:
            ipd_mean[r][k] = (r / (1 + phi * (r - 1))).mean()
    return n_mean, n_p5, n_p95, ipd_mean


def plot_panel_B(ax):
    sig = LN_SIGMA_GRID
    n_mean, n_p5, n_p95, ipd_mean = compute_dominance_experiment(sig)

    # Left axis: dominant-motif estimated IPD ratio (phi = 1/n covered genomes).
    rows = []
    for r in R_VALUES:
        y = ipd_mean[r]
        ax.plot(sig, y, color=COLORS[r], lw=2, label=_legend_label(r))
        rows.extend([(si, nm, p5, p95, r, yi)
                     for si, nm, p5, p95, yi
                     in zip(sig, n_mean, n_p5, n_p95, y)])
    ax.axhline(1, ls="--", color="#888888", lw=1)
    ax.set_xlabel("Community unevenness (log-normal $\\sigma$)")
    ax.set_ylabel("Estimated IPD ratio")
    ax.set_ylim(1, max(R_VALUES) + 0.1)
    ax.set_xlim(sig[0], sig[-1])
    ax.set_title("c", loc="left", fontweight="bold")

    # Right axis: number of covered genomes (depth >= 5x).
    NEUTRAL = "#333333"
    ax2 = ax.twinx()
    ax2.fill_between(sig, n_p5, n_p95, color=NEUTRAL, alpha=0.15, linewidth=0)
    ax2.plot(sig, n_mean, color=NEUTRAL, lw=2.5,
             label="genomes with depth $\\geq$ 5x")
    ax2.set_ylabel(r"Genomes with depth $\geq$ 5x", color=NEUTRAL)
    ax2.tick_params(axis="y", labelcolor=NEUTRAL)
    ax2.set_ylim(0, N_GENOMES + 3)
    ax2.grid(False)
    ax2.legend(loc="upper right", fontsize=9, frameon=True)

    return _save_sourcedata(
        rows,
        ["ln_sigma", "n_covered_mean", "n_covered_p5", "n_covered_p95",
         "r", "dominant_est_ipdratio"],
        "dominance_experiment")


# --------------------------------------------------------------------------- #
# Panel C: estimated ipdRatio (r / sigma) vs modified k-mer fraction phi
# --------------------------------------------------------------------------- #
def plot_panel_C(ax):
    phi = np.linspace(0.0, 1.0, 201)  # full range: at phi=1 all instances modified
    rows = []
    for r in R_VALUES:
        s = sigma_phi(phi, r)
        y = r / s  # estimated ipdRatio = h/a = r/sigma
        ax.plot(phi, y, linewidth=2, color=COLORS[r], label=_legend_label(r))
        rows.extend([(pi, r, si, yi) for pi, si, yi in zip(phi, s, y)])
    ax.axhline(1, ls="--", color="#888888", lw=1)  # detection floor
    ax.set_xlabel(r"Community modification fraction $\phi$")
    ax.set_ylabel("Estimated IPD ratio")
    ax.set_title("b", loc="left", fontweight="bold")
    return _save_sourcedata(rows, ["phi", "r", "sigma", "estimated_ipdRatio"],
                            "est_ipdratio_vs_phi")


# --------------------------------------------------------------------------- #
# Standalone figure helper
# --------------------------------------------------------------------------- #
def standalone(plot_fn, basename):
    fig, ax = plt.subplots(figsize=(8, 5))
    plot_fn(ax)
    ax.legend(fontsize=9)
    fig.tight_layout()
    _save(fig, basename)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Console summary (copy-paste ready)
# --------------------------------------------------------------------------- #
def print_summary():
    line = "=" * 70
    print("\n" + line)
    print("Control IPD deviation sigma = 1 + phi*(r-1);  single genome phi = t/X")
    print(line)

    print(f"\nPanel a  sigma vs metagenome length X  (t = {T_FIXED:.0f} Mbp):")
    print(f"  at X = {X_FIXED:.0f} Mbp (phi = {T_FIXED / X_FIXED:.3f}):")
    for r in R_VALUES:
        print(f"    r = {r}:  sigma = {sigma_single(T_FIXED, X_FIXED, r):.4f}")

    print(f"\nPanel b  dominance experiment  (N = {N_GENOMES} genomes x "
          f"{GENOME_LEN:.0f} Mbp, seq = {TOTAL_SEQ / 1000:.0f} Gbp, "
          f"mean depth = {MEAN_DEPTH:.0f}x, floor = {DEPTH_FLOOR:.0f}x):")
    sig_pts = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
    n_mean, _, _, ipd_mean = compute_dominance_experiment(sig_pts)
    hdr = "  sigma  genomes>=5x  phi=1/n " + "".join(
        f"  ipd(r={r})" for r in R_VALUES)
    print(hdr)
    for i, sg in enumerate(sig_pts):
        vals = "".join(f"  {ipd_mean[r][i]:<8.3f}" for r in R_VALUES)
        print(f"  {sg:<6.1f} {n_mean[i]:<11.1f}  {1 / n_mean[i]:<7.4f}{vals}")

    print("\nPanel c  estimated IPD ratio = r/sigma vs community fraction phi:")
    header = "  phi   " + "".join(f"  r={r:<5}" for r in R_VALUES)
    print(header)
    for phi in [0.0, 0.05, 0.1, 0.25, 0.5, 1.0]:
        vals = "".join(f"  {r / sigma_phi(phi, r):<6.4f}" for r in R_VALUES)
        print(f"  {phi:<6.2f}{vals}")
    print(line + "\n")


def check_consistency():
    """Sanity checks on the analytic panels and the experiment."""
    # Paper anchor: X=100 gives sigma in [1.025, 1.10] for r in [1.5, 3].
    assert abs(sigma_single(T_FIXED, X_FIXED, 1.5) - 1.025) < 1e-9
    assert abs(sigma_single(T_FIXED, X_FIXED, 3) - 1.10) < 1e-9
    # Panel c at phi=0 recovers the true ratio r; at phi=1 collapses to 1.
    for r in R_VALUES:
        assert abs(r / sigma_phi(0.0, r) - r) < 1e-12
        assert abs(r / sigma_phi(1.0, r) - 1.0) < 1e-12
    # Dominance experiment: covered-genome count falls as dominance (sigma) rises,
    # and the estimated IPD ratio declines from near r toward 1.
    n_mean, _, _, ipd_mean = compute_dominance_experiment(LN_SIGMA_GRID)
    assert np.all(np.diff(n_mean) <= 1.0), "n(sigma) must trend down (MC jitter ok)"
    assert n_mean[0] - n_mean[-1] > 10, "covered genomes must drop as sigma grows"
    for r in R_VALUES:
        assert ipd_mean[r][0] > ipd_mean[r][-1], "IPD ratio must fall with sigma"
        assert 1.0 <= ipd_mean[r][-1] <= r + 1e-9


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    check_consistency()

    print("Writing standalone figures + source data:")
    standalone(plot_panel_A, "sigma_vs_metagenome_length")
    standalone(plot_panel_C, "est_ipdratio_vs_phi")

    # Panel B standalone (twin-axis; do not add the generic single legend).
    fig, ax = plt.subplots(figsize=(8, 5))
    plot_panel_B(ax)
    fig.tight_layout()
    _save(fig, "dominance_experiment")
    plt.close(fig)

    # Combined 3x1 panel with a single shared legend.
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    plot_panel_A(axes[0])
    plot_panel_C(axes[1])   # community modification fraction phi -> panel b
    plot_panel_B(axes[2])   # dominance experiment -> panel c
    axes[0].legend(loc="upper right", fontsize=10, frameon=True)
    fig.tight_layout()
    print("Writing combined 3-panel figure:")
    _save(fig, "control_ipd_deviation_3panel")
    plt.close(fig)

    print_summary()


if __name__ == "__main__":
    main()
