"""
Does downstream (3') context add predictive signal for IPD, and at what k-mer
occurrence cutoff does it become detectable?

Fix the upstream context at u=U_FIX (=7, the informative extent). For each downstream
extent d and each occurrence cutoff c, compare -- ON THE SAME SITES -- two predictors:
    * full : the (U_FIX, d) k-mer mean   (upstream + d downstream bases)
    * up   : the (U_FIX, 0) k-mer mean   (upstream only)
scored on the subset of test sites whose full (U_FIX,d) k-mer has >= c training
occurrences (so the downstream-extended k-mer is actually estimable). The matched
subset removes the population/selection confound: the only difference between the two
predictors is whether the d downstream bases are used.

    dR2(c,d) = R2_full - R2_up   (incremental value of the downstream bases)

If downstream carries signal that we can estimate, dR2 > 0. Too-high a cutoff -> the
extended k-mers are never observed (can't test); too-low -> their means are noisy and
dR2 can go negative. The sweep shows which cutoff, if any, reveals a downstream effect.

Run:  conda run -n modifi python benchmark/check_context/downstream_cutoff_test.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import one_way_anova as owa

U_FIX = 7
D_RANGE = range(0, 6)
CUTOFFS = [2, 3, 5, 10, 20, 50]
SEED = owa.BASE_SEED
FIGDIR = owa.FIGDIR
OUTDIR = owa.OUTDIR


def main():
    contig, ref, comp = owa.load_reference(owa.c224_contig)
    df = pd.read_csv(owa.CACHE)
    fwd = df["strand"].values == "+"
    pos_f = df.loc[fwd, "tpl"].values.astype(np.int64)
    y_f = df.loc[fwd, "y"].values.astype(np.float64)
    pos_r = df.loc[~fwd, "tpl"].values.astype(np.int64)
    y_r = df.loc[~fwd, "y"].values.astype(np.float64)
    print(f"contig={contig}  fixing upstream u={U_FIX}", flush=True)

    rng = np.random.default_rng(SEED)
    rows = []
    # per-d: codes + train means depend on d (not cutoff) -> compute once, sweep cutoffs
    for d in D_RANGE:
        pf, yf, pr, yr = owa.valid_sites(ref, comp, pos_f, y_f, pos_r, y_r, U_FIX, d)
        y = np.concatenate([yf, yr])
        n = y.shape[0]
        idx = rng.permutation(n)
        tr, te = idx[: n // 2], idx[n // 2:]
        gmean = y[tr].mean()

        codes_full = owa.codes_at(ref, comp, pf, pr, U_FIX, d)
        codes_up = owa.codes_at(ref, comp, pf, pr, U_FIX, 0)
        uf, invf = np.unique(codes_full, return_inverse=True)
        uu, invu = np.unique(codes_up, return_inverse=True)
        sf = np.bincount(invf[tr], weights=y[tr], minlength=uf.size)
        cf = np.bincount(invf[tr], minlength=uf.size)
        su = np.bincount(invu[tr], weights=y[tr], minlength=uu.size)
        cu = np.bincount(invu[tr], minlength=uu.size)
        with np.errstate(invalid="ignore", divide="ignore"):
            mf = sf / cf
            mu = su / cu
        cnt_te = cf[invf[te]]
        y_te = y[te]
        pred_full = mf[invf[te]]
        pred_up = mu[invu[te]]

        for c in CUTOFFS:
            keep = cnt_te >= c
            ret = float(keep.mean())
            if keep.sum() < 2:
                r2_full = r2_up = np.nan
            else:
                yk = y_te[keep]
                mk = yk.mean()
                sstot = np.sum((yk - mk) ** 2)
                r2_full = 1.0 - np.sum((yk - pred_full[keep]) ** 2) / sstot
                r2_up = 1.0 - np.sum((yk - pred_up[keep]) ** 2) / sstot
            rows.append(dict(d=d, cutoff=c, retained=ret,
                             r2_full=r2_full, r2_up=r2_up,
                             dR2=r2_full - r2_up))
            print(f"d={d} cutoff={c:>2} retained={ret:.3f} "
                  f"R2_full={r2_full:.4f} R2_up={r2_up:.4f} "
                  f"dR2={r2_full - r2_up:+.4f}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUTDIR, "downstream_cutoff_test.csv"), index=False)
    plot(res)
    print("wrote downstream_cutoff_test.{png,pdf} to", FIGDIR, flush=True)


def plot(res):
    cutoffs = sorted(res["cutoff"].unique())
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(cutoffs)))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

    for c, col in zip(cutoffs, colors):
        sub = res[res["cutoff"] == c].sort_values("d")
        ax1.plot(sub["d"], sub["dR2"], "-o", color=col, lw=2, ms=5,
                 label=f"cutoff {c}")
        ax2.plot(sub["d"], sub["retained"], "-o", color=col, lw=2, ms=5,
                 label=f"cutoff {c}")

    ax1.axhline(0, color="0.5", lw=1, ls="--")
    ax1.set_xlabel("downstream bases added  d  (on top of u=7)")
    ax1.set_ylabel("$\\Delta R^2$ = R²(u=7,d) − R²(u=7,0)\n(same sites)")
    ax1.set_title("Incremental value of downstream context")
    ax1.legend(fontsize=8, title="k-mer occurrence\ncutoff", title_fontsize=8)

    ax2.set_xlabel("downstream bases added  d")
    ax2.set_ylabel("fraction of sites retained")
    ax2.set_title("Sites where the (7,d) k-mer meets the cutoff")
    ax2.set_yscale("log")

    fig.suptitle("Downstream-context contribution vs. k-mer occurrence cutoff (C227, u=7)",
                 fontsize=11)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, "downstream_cutoff_test.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
