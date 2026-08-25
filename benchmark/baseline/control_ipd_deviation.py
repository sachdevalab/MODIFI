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
# Panel A: sigma vs metagenome length X (recap of Supplementary Figure S15)
# --------------------------------------------------------------------------- #
def plot_panel_A(ax):
    X = np.arange(5, 3005, 5, dtype=float)  # Mbp
    rows = []
    for r in R_VALUES:
        y = sigma_single(T_FIXED, X, r)
        ax.plot(X, y, linewidth=2, color=COLORS[r], label=_legend_label(r))
        rows.extend([(xi, r, yi) for xi, yi in zip(X, y)])
    ax.axhline(1, ls="--", color="#888888", lw=1)
    ax.set_xlabel("Metagenome length X (Mbp)")
    ax.set_ylabel(r"Control IPD deviation $\sigma$")
    ax.set_ylim(bottom=1)
    ax.set_title("a", loc="left", fontweight="bold")
    return _save_sourcedata(rows, ["metagenome_length_Mbp", "r", "sigma"],
                            "sigma_vs_metagenome_length")


# --------------------------------------------------------------------------- #
# Panel B: estimated ipdRatio (r / sigma) vs dominant genome length t
# --------------------------------------------------------------------------- #
def plot_panel_B(ax):
    t = np.linspace(1, 15, 281)  # Mbp
    rows = []
    for r in R_VALUES:
        s = sigma_single(t, X_FIXED, r)
        y = r / s  # estimated ipdRatio = h/a = r/sigma
        ax.plot(t, y, linewidth=2, color=COLORS[r], label=_legend_label(r))
        rows.extend([(ti, r, si, yi) for ti, si, yi in zip(t, s, y)])
    ax.axhline(1, ls="--", color="#888888", lw=1)  # detection floor
    ax.set_xlabel("Dominant genome length t (Mbp)")
    ax.set_ylabel("Estimated IPD ratio")
    ax.set_title("b", loc="left", fontweight="bold")
    return _save_sourcedata(rows, ["genome_length_Mbp", "r", "sigma",
                                    "estimated_ipdRatio"],
                            "est_ipdratio_vs_genome_length")


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
    ax.set_title("c", loc="left", fontweight="bold")
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

    print(f"\nPanel A/B anchor  (t = {T_FIXED:.0f} Mbp, X = {X_FIXED:.0f} Mbp,"
          f" phi = {T_FIXED / X_FIXED:.3f}):")
    for r in R_VALUES:
        print(f"  r = {r}:  sigma = {sigma_single(T_FIXED, X_FIXED, r):.4f}")

    print("\nPanel B  estimated ipdRatio = r/sigma vs dominant genome length t"
          "  (X = 100 Mbp):")
    header = "  t(Mbp) " + "".join(f"  r={r:<5}" for r in R_VALUES)
    print(header)
    for t in [1, 5, 10, 15]:
        vals = "".join(f"  {r / sigma_single(t, X_FIXED, r):<6.4f}"
                       for r in R_VALUES)
        print(f"  {t:<6d}{vals}")

    print("\nPanel C  estimated ipdRatio = r/sigma vs modified fraction phi:")
    print(header.replace("t(Mbp)", "phi   "))
    for phi in [0.0, 0.05, 0.1, 0.25, 0.5, 1.0]:
        vals = "".join(f"  {r / sigma_phi(phi, r):<6.4f}" for r in R_VALUES)
        print(f"  {phi:<6.2f}{vals}")
    print(line + "\n")


def check_consistency():
    """Panels B and C must agree at the shared 0.05 fraction point."""
    b = sigma_single(T_FIXED, X_FIXED, 3)      # panel B point at t=5
    c = sigma_phi(T_FIXED / X_FIXED, 3)        # panel C point at phi=0.05
    assert abs(b - c) < 1e-12, (b, c)
    # Paper anchor: X=100 gives sigma in [1.025, 1.10] for r in [1.5, 3].
    assert abs(sigma_single(T_FIXED, X_FIXED, 1.5) - 1.025) < 1e-9
    assert abs(sigma_single(T_FIXED, X_FIXED, 3) - 1.10) < 1e-9
    # Panels B/C at phi=0 recover the true ratio r; at phi=1 collapse to 1.
    for r in R_VALUES:
        assert abs(r / sigma_phi(0.0, r) - r) < 1e-12
        assert abs(r / sigma_phi(1.0, r) - 1.0) < 1e-12


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    check_consistency()

    print("Writing standalone figures + source data:")
    standalone(plot_panel_A, "sigma_vs_metagenome_length")
    standalone(plot_panel_B, "est_ipdratio_vs_genome_length")
    standalone(plot_panel_C, "est_ipdratio_vs_phi")

    # Combined 3x1 panel with a single shared legend.
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    plot_panel_A(axes[0])
    plot_panel_B(axes[1])
    plot_panel_C(axes[2])
    axes[0].legend(loc="upper right", fontsize=10, frameon=True)
    fig.tight_layout()
    print("Writing combined 3-panel figure:")
    _save(fig, "control_ipd_deviation_3panel")
    plt.close(fig)

    print_summary()


if __name__ == "__main__":
    main()
