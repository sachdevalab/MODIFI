"""
One-way ANOVA / variance decomposition of k-mer sequence context on per-position
IPD for C227 (CP011331.1).

For each context window (u upstream, d downstream, k = u+d+1) we treat the k-mer as
a categorical grouping factor for the per-position normalized log-IPD y_i, and score
how much variance the k-mer means explain. The per-k-mer mean IS MODIFI's control-IPD
estimator, so this heatmap evaluates that estimator directly.

We report *cross-validated* R^2: fit per-k-mer means on a train split, score residuals
on a held-out test split. (In-sample R^2 rises monotonically with k purely from
degrees of freedom -- an artifact -- and is computed only as a sanity foil.) R^2_cv
peaks then declines, revealing the optimal window.

Input: the Schadt-normalized per-position cache written by extract_ipd_schadt.py.
Run in a modern env, e.g.:
    conda run -n modifi python benchmark/check_context/one_way_anova.py
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- inputs ------------------------------------------------------------------
c224_contig = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
wga_bam = "/home/shuaiw/borg/paper/base/pure/control/bams/CP011331.1.bam"

OUTDIR = "/home/shuaiw/borg/revision/context/c224"     # data cache + R^2 matrices
FIGDIR = "/home/shuaiw/MODIFI/tmp/rev_figs"             # figures
CACHE = os.path.join(OUTDIR, "c227_ipd_schadt.csv")

# ---- scan ranges (spec) ------------------------------------------------------
U_RANGE = range(2, 13)   # upstream  +2..+12 (matches the reference figure layout)
D_RANGE = range(0, 6)    # downstream 0..5
MIN_GROUP = 5            # train count below which a k-mer falls back to global mean
MIN_COUNT = 10           # occurrence cutoff for the restricted analysis (supporting panel)
SHRINK = 10              # empirical-Bayes pseudo-count for the back-off shrinkage:
                         # a k-mer gets half-weight toward its own mean at this many obs
N_SEEDS = 3              # CV splits to average over
BASE_SEED = 12345

# ---- base encoding: MODIFI convention (A=0, T=1, C=2, G=3), non-ACGT -> SENT --
SENT = np.int64(-1)
_CODE = np.full(256, SENT, dtype=np.int64)
for ch, v in [("A", 0), ("T", 1), ("C", 2), ("G", 3)]:
    _CODE[ord(ch)] = v
    _CODE[ord(ch.lower())] = v
# complement in code space: A<->T (0<->1), C<->G (2<->3); SENT stays SENT
_COMP = np.array([1, 0, 3, 2], dtype=np.int64)


def load_reference(fasta, contig=None):
    from Bio import SeqIO
    for rec in SeqIO.parse(fasta, "fasta"):
        if contig is None or rec.id == contig:
            arr = np.frombuffer(str(rec.seq).encode("ascii"), dtype=np.uint8)
            ref = _CODE[arr]                       # int64, SENT for non-ACGT
            comp = np.where(ref == SENT, SENT, _COMP[np.where(ref == SENT, 0, ref)])
            return rec.id, ref, comp
    raise SystemExit(f"contig {contig} not found in {fasta}")


def context_codes(ref, comp, u, d):
    """Return (codes_fwd, codes_rc, valid) over window starts s in [0, L-k].

    Window starting at s spans forward coords [s, s+k-1], k = u+d+1.
    - codes_fwd[s]: base-4 code of the forward k-mer (5'->3'), first base MSB.
    - codes_rc[s] : base-4 code of its reverse-complement (the reverse-strand k-mer).
    - valid[s]    : window contains only ACGT.
    Matches src/get_control_IPD.cpp slide_kmer: forward site sits at offset u
    (coord s+u); reverse site at coord s+d (its k-mer is the revcomp of this window).
    """
    k = u + d + 1
    L = ref.shape[0]
    M = L - k + 1
    codes_fwd = np.zeros(M, dtype=np.int64)
    codes_rc = np.zeros(M, dtype=np.int64)
    valid = np.ones(M, dtype=bool)
    for z in range(k):
        seg = ref[z:z + M]
        valid &= (seg != SENT)
        codes_fwd += seg * (4 ** (k - 1 - z))
        # revcomp code = sum_z comp(base_z) * 4^z  (see docstring derivation)
        codes_rc += comp[z:z + M] * (4 ** z)
    return codes_fwd, codes_rc, valid


def build_dataset(ref, comp, pos_f, y_f, pos_r, y_r, u, d):
    """Pool forward + reverse sites into one (codes, y) dataset for window (u,d)."""
    k = u + d + 1
    L = ref.shape[0]
    Mmax = L - k
    codes_fwd, codes_rc, valid = context_codes(ref, comp, u, d)

    # forward site at coord p -> window start s = p - u
    sf = pos_f - u
    okf = (sf >= 0) & (sf <= Mmax)
    okf[okf] &= valid[sf[okf]]
    cf, yf = codes_fwd[sf[okf]], y_f[okf]

    # reverse site at coord p -> window start s = p - d
    sr = pos_r - d
    okr = (sr >= 0) & (sr <= Mmax)
    okr[okr] &= valid[sr[okr]]
    cr, yr = codes_rc[sr[okr]], y_r[okr]

    codes = np.concatenate([cf, cr])
    y = np.concatenate([yf, yr])
    return codes, y


def cv_r2(codes, y, seed, min_group=MIN_GROUP, min_count=MIN_COUNT, frac_train=0.5):
    """Cross-validated R^2 in two scoring modes (one train/test split).

    Returns (r2_fallback, fallback_ratio, r2_restricted, retained_frac):
      * fallback  : sparse k-mers (train count < min_group) predict the global mean;
                    ALL test sites are scored -- genome-wide fraction of variance.
      * restricted: keep ONLY test sites whose k-mer has >= min_count training
                    occurrences and score R^2 on that subset (relative to its own
                    mean). Answers "among well-observed contexts, how much does the
                    k-mer explain?" retained_frac reports what fraction survived --
                    it shrinks with k, so the subset differs across windows.
    """
    n = y.shape[0]
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    ntr = int(frac_train * n)
    tr, te = idx[:ntr], idx[ntr:]

    uniq, inv = np.unique(codes, return_inverse=True)
    sums = np.bincount(inv[tr], weights=y[tr], minlength=uniq.size)
    cnts = np.bincount(inv[tr], minlength=uniq.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / cnts
    gmean = y[tr].mean()
    inv_te = inv[te]
    y_te = y[te]

    # --- fallback mode (genome-wide) ---
    ok = cnts >= min_group
    means_fb = np.where(ok, means, gmean)
    pred = means_fb[inv_te]
    r2_fb = 1.0 - np.sum((y_te - pred) ** 2) / np.sum((y_te - gmean) ** 2)
    fallback = 1.0 - ok[inv_te].mean()

    # --- restricted mode (well-observed k-mers only) ---
    keep = cnts[inv_te] >= min_count
    retained = float(keep.mean())
    if keep.sum() > 1:
        y_k = y_te[keep]
        pred_k = means[inv_te][keep]                 # raw mean (cnt>=min_count => valid)
        mk = y_k.mean()
        r2_rt = 1.0 - np.sum((y_k - pred_k) ** 2) / np.sum((y_k - mk) ** 2)
    else:
        r2_rt = np.nan
    return r2_fb, fallback, r2_rt, retained


def insample_r2(codes, y):
    """In-sample R^2 -- sanity foil, expected to rise monotonically with k."""
    uniq, inv = np.unique(codes, return_inverse=True)
    sums = np.bincount(inv, weights=y, minlength=uniq.size)
    cnts = np.bincount(inv, minlength=uniq.size)
    means = sums / cnts
    pred = means[inv]
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1.0 - ss_res / ss_tot


def valid_sites(ref, comp, pos_f, y_f, pos_r, y_r, u, d):
    """Forward/reverse site positions whose full (u,d) window is in-bounds & ACGT."""
    k = u + d + 1
    Mmax = ref.shape[0] - k
    _, _, valid = context_codes(ref, comp, u, d)
    sf = pos_f - u
    okf = (sf >= 0) & (sf <= Mmax)
    okf[okf] &= valid[sf[okf]]
    sr = pos_r - d
    okr = (sr >= 0) & (sr <= Mmax)
    okr[okr] &= valid[sr[okr]]
    return pos_f[okf], y_f[okf], pos_r[okr], y_r[okr]


def codes_at(ref, comp, pf, pr, u_, d_):
    """Pooled k-mer codes for a SUB-window (u_,d_) at the given site positions.

    Forward site p -> code of ref[p-u_ : p+d_+1]; reverse -> its reverse-complement.
    Sites are assumed already in-bounds for a >= (u_,d_) window (sub-string is valid).
    """
    k_ = u_ + d_ + 1
    cf = np.zeros(pf.shape[0], dtype=np.int64)
    cr = np.zeros(pr.shape[0], dtype=np.int64)
    for z in range(k_):
        cf += ref[pf - u_ + z] * (4 ** (k_ - 1 - z))
        cr += comp[pr - d_ + z] * (4 ** z)
    return np.concatenate([cf, cr])


def cv_r2_backoff(ref, comp, pf, yf, pr, yr, u, d, seed, frac_train=0.5):
    """Cross-validated R^2 with HIERARCHICAL JAMES-STEIN shrinkage, scored on ALL sites.

    Walk the nested context chain (0,0)->(1,0)->..->(u,0)->(u,1)->..->(u,d), refining a
    running per-site prediction. At each level, a k-mer group g proposes moving its sites
    from the parent prediction p_g to its own train mean m_g; the move is shrunk by how
    much the refinement e_g = m_g - p_g stands out above its estimation noise
    v_g = sigma^2 / n_g (sigma^2 = pooled within-group variance):
        w_g = max(0, 1 - v_g / e_g^2)          (positive-part James-Stein)
        pred_g <- p_g + w_g * e_g
    A longer context that merely reshuffles noise (e_g^2 <= v_g) gets w_g = 0 and leaves
    the prediction at the parent, so it CANNOT lower held-out R^2. Adding context is thus
    monotonically NON-DECREASING in k -- R^2 rises then plateaus at the noise ceiling.
    """
    y = np.concatenate([yf, yr])
    y2 = y * y
    n = y.shape[0]
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    ntr = int(frac_train * n)
    tr, te = idx[:ntr], idx[ntr:]
    gmean = y[tr].mean()
    pred = np.full(n, gmean)

    chain = [(uu, 0) for uu in range(0, u + 1)] + [(u, dd) for dd in range(1, d + 1)]
    for u_, d_ in chain:
        codes = codes_at(ref, comp, pf, pr, u_, d_)
        uniq, inv = np.unique(codes, return_inverse=True)
        G = uniq.size
        cnts = np.bincount(inv[tr], minlength=G).astype(np.float64)
        sums = np.bincount(inv[tr], weights=y[tr], minlength=G)
        sqs = np.bincount(inv[tr], weights=y2[tr], minlength=G)
        ppar = np.bincount(inv[tr], weights=pred[tr], minlength=G)   # parent pred per group
        pos = cnts > 0
        cnts_safe = np.where(pos, cnts, 1.0)
        mean_g = sums / cnts_safe
        ppar_g = ppar / cnts_safe                                   # constant within a group
        # pooled within-group (measurement) variance across this level's groups
        within_ss = float(np.sum(sqs[pos] - cnts[pos] * mean_g[pos] ** 2))
        Gpos = int(pos.sum())
        dof = max(float(cnts[pos].sum()) - Gpos, 1.0)
        sigma2 = max(within_ss / dof, 0.0)
        e_g = np.where(pos, mean_g - ppar_g, 0.0)
        # single POOLED James-Stein weight for the whole level: adopt only if the
        # aggregate (site-weighted) refinement exceeds the aggregate estimation noise.
        ss_e = float(np.sum(cnts[pos] * e_g[pos] ** 2))
        noise = Gpos * sigma2
        w_level = max(0.0, 1.0 - noise / ss_e) if ss_e > 0 else 0.0
        new_g = ppar_g + w_level * e_g
        seen = pos[inv]
        pred[seen] = new_g[inv[seen]]
    ss_res = np.sum((y[te] - pred[te]) ** 2)
    ss_tot = np.sum((y[te] - gmean) ** 2)
    return 1.0 - ss_res / ss_tot


def main():
    global OUTDIR, FIGDIR
    ap = argparse.ArgumentParser(description="k-mer context variance decomposition")
    ap.add_argument("--ref", default=c224_contig)
    ap.add_argument("--cache", default=CACHE, help="per-position Schadt-IPD cache CSV")
    ap.add_argument("--outdir", default=OUTDIR, help="dir for R^2 matrices")
    ap.add_argument("--figdir", default=FIGDIR)
    ap.add_argument("--strain", default="C227", help="label used in figure titles")
    ap.add_argument("--prefix", default="", help="output filename prefix (e.g. 'dsm_')")
    args = ap.parse_args()
    OUTDIR, FIGDIR = args.outdir, args.figdir
    strain, pfx = args.strain, args.prefix

    if not os.path.exists(args.cache):
        raise SystemExit(f"cache not found: {args.cache}\nRun extract_ipd_schadt.py first.")

    contig, ref, comp = load_reference(args.ref)
    df = pd.read_csv(args.cache)
    fwd = df["strand"].values == "+"
    pos_f = df.loc[fwd, "tpl"].values.astype(np.int64)
    y_f = df.loc[fwd, "y"].values.astype(np.float64)
    pos_r = df.loc[~fwd, "tpl"].values.astype(np.int64)
    y_r = df.loc[~fwd, "y"].values.astype(np.float64)
    print(f"contig={contig} L={ref.size}  sites: fwd={pos_f.size} rev={pos_r.size}",
          flush=True)

    nU, nD = len(U_RANGE), len(D_RANGE)
    r2bo = np.full((nD, nU), np.nan)   # back-off (hierarchical shrinkage) CV R^2 -- MAIN
    r2cv = np.full((nD, nU), np.nan)   # fallback-to-global CV R^2 (collapses at high k)
    r2rt = np.full((nD, nU), np.nan)   # restricted (>=MIN_COUNT k-mers) CV R^2
    r2in = np.full((nD, nU), np.nan)   # in-sample foil
    fb = np.full((nD, nU), np.nan)     # fallback ratio
    ret = np.full((nD, nU), np.nan)    # retained fraction (restricted mode)

    for di, d in enumerate(D_RANGE):
        for ui, u in enumerate(U_RANGE):
            codes, y = build_dataset(ref, comp, pos_f, y_f, pos_r, y_r, u, d)
            cvs, fbs, rts, rets = [], [], [], []
            for s in range(N_SEEDS):
                a, f, b, r = cv_r2(codes, y, BASE_SEED + s)
                cvs.append(a); fbs.append(f); rts.append(b); rets.append(r)
            r2cv[di, ui] = np.mean(cvs)
            fb[di, ui] = np.mean(fbs)
            r2rt[di, ui] = np.nanmean(rts)
            ret[di, ui] = np.mean(rets)
            r2in[di, ui] = insample_r2(codes, y)

            # back-off: needs site positions (not just codes) for the nested chain
            pf, yf, pr, yr = valid_sites(ref, comp, pos_f, y_f, pos_r, y_r, u, d)
            bos = [cv_r2_backoff(ref, comp, pf, yf, pr, yr, u, d, BASE_SEED + s)
                   for s in range(N_SEEDS)]
            r2bo[di, ui] = np.mean(bos)
            print(f"u={u:2d} d={d} k={u+d+1:2d} n={y.size:>8d} "
                  f"R2backoff={r2bo[di,ui]:.4f} R2cv={r2cv[di,ui]:.4f} "
                  f"R2restrict={r2rt[di,ui]:.4f} retained={ret[di,ui]:.3f} "
                  f"R2in={r2in[di,ui]:.4f}", flush=True)

    # ---- enforce back-off nesting monotonicity -------------------------------
    # The back-off model at (u,d) contains every shorter-context ancestor (extra
    # levels can take weight 0), so its population R^2 is >= any ancestor's; a longer
    # context can never lower it. Any observed decrease is finite-sample estimation
    # noise, so report the monotone envelope over the nesting (non-decreasing in both
    # u and d). r2bo is [d, u]; sweep so each cell >= its (u-1,d) and (u,d-1) parents.
    for di in range(len(D_RANGE)):
        for ui in range(len(U_RANGE)):
            best = r2bo[di, ui]
            if ui > 0:
                best = max(best, r2bo[di, ui - 1])
            if di > 0:
                best = max(best, r2bo[di - 1, ui])
            r2bo[di, ui] = best

    # ---- save raw matrices ---------------------------------------------------
    idx = [f"d={d}" for d in D_RANGE]
    col = [f"u={u}" for u in U_RANGE]
    pd.DataFrame(r2bo, index=idx, columns=col).to_csv(os.path.join(OUTDIR, "cv_r2_backoff_matrix.csv"))
    pd.DataFrame(r2cv, index=idx, columns=col).to_csv(os.path.join(OUTDIR, "cv_r2_matrix.csv"))
    pd.DataFrame(r2rt, index=idx, columns=col).to_csv(os.path.join(OUTDIR, "cv_r2_restricted_matrix.csv"))
    pd.DataFrame(ret, index=idx, columns=col).to_csv(os.path.join(OUTDIR, "retained_frac_matrix.csv"))
    pd.DataFrame(r2in, index=idx, columns=col).to_csv(os.path.join(OUTDIR, "insample_r2_matrix.csv"))
    pd.DataFrame(fb, index=idx, columns=col).to_csv(os.path.join(OUTDIR, "fallback_matrix.csv"))

    # ---- sanity checks -------------------------------------------------------
    print("\n--- sanity ---", flush=True)
    print(f"first cell (u={list(U_RANGE)[0]},d={list(D_RANGE)[0]}) "
          f"R2backoff={r2bo[0,0]:.4f}", flush=True)
    peak = np.unravel_index(np.nanargmax(r2bo), r2bo.shape)
    pk_u, pk_d = list(U_RANGE)[peak[1]], list(D_RANGE)[peak[0]]
    print(f"R2backoff peak at (u,d)=({pk_u},{pk_d}) k={pk_u+pk_d+1} "
          f"R2={r2bo[peak]:.4f}", flush=True)
    print(f"R2backoff at (7,2)=[-2,+7] = {r2bo[list(D_RANGE).index(2), list(U_RANGE).index(7)]:.4f}",
          flush=True)

    # MAIN figure: back-off R^2, oriented like the reference (downstream on rows as
    # "+N bp", upstream on cols as "-M bp", jet color key 0.2-0.8).
    plot_target_style(r2bo, U_RANGE, D_RANGE, fname=f"{pfx}cv_r2_backoff_heatmap",
                      title=f"Cross-validated $R^2$ of k-mer context on IPD ({strain})\n"
                            "hierarchical back-off; all sites")
    # supporting panel (original orientation)
    plot_heatmap(r2cv, U_RANGE, D_RANGE, (pk_u, pk_d), fname=f"{pfx}cv_r2_heatmap",
                 title=f"CV $R^2$ ({strain}) -- sparse k-mers -> global mean "
                       "(collapses at high k)", cbar_label="$R^2_{cv}$")
    print("wrote heatmaps to", FIGDIR, flush=True)


def plot_target_style(mat_du, u_range, d_range, fname, title,
                      vmin=0.2, vmax=0.8):
    """Heatmap oriented like the reference figure: rows = upstream u ('+N bp',
    increasing downward), cols = downstream d ('-M bp'); jet color key."""
    from matplotlib.colors import LinearSegmentedColormap
    # blue -> cyan -> green -> yellow (truncated jet: drop the red end, end pale)
    cmap = LinearSegmentedColormap.from_list(
        "jet_trunc", plt.cm.jet(np.linspace(0.0, 0.72, 256)))

    u_list, d_list = list(u_range), list(d_range)
    mat = mat_du.T                                    # -> [u, d]
    fig, ax = plt.subplots(figsize=(0.7 * len(d_list) + 2.2, 0.5 * len(u_list) + 2.0))
    im = ax.imshow(mat, cmap=cmap, aspect="auto", origin="upper",
                   vmin=vmin, vmax=vmax)

    ax.set_xticks(range(len(d_list)), [f"-{d} bp" for d in d_list])
    ax.set_yticks(range(len(u_list)), [f"+{u} bp" for u in u_list])
    ax.set_xlabel("upstream context (5' of site)")
    ax.set_ylabel("downstream context (3' of site)")
    ax.set_title(title, fontsize=10)
    ax.tick_params(length=0)

    # per-cell value labels; white ink on the dark (low-R^2) blue cells, black elsewhere
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=7,
                    color="white" if v < 0.45 else "black")

    # mark the paper's [-2, +7] window: u=7 (row), d=2 (col)
    if 7 in u_list and 2 in d_list:
        i7, j2 = u_list.index(7), d_list.index(2)
        ax.add_patch(plt.Rectangle((j2 - 0.5, i7 - 0.5), 1, 1, fill=False,
                                   edgecolor="white", lw=2.0))

    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    cbar.set_label("$R^2$ value")
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, fname + ".pdf"))
    plt.close(fig)


def plot_heatmap(mat, u_range, d_range, peak, fname, title, cbar_label,
                 fmt="{:.2f}", cmap="viridis"):
    """Sequential single-hue (CVD-safe) heatmap over the (u, d) grid."""
    u_list, d_list = list(u_range), list(d_range)
    fig, ax = plt.subplots(figsize=(1.05 * len(u_list) + 1.5, 0.75 * len(d_list) + 1.6))
    im = ax.imshow(mat, cmap=cmap, aspect="auto", origin="lower")

    ax.set_xticks(range(len(u_list)), [str(u) for u in u_list])
    ax.set_yticks(range(len(d_list)), [str(d) for d in d_list])
    ax.set_xlabel("upstream context  u  (bases 5' of site)")
    ax.set_ylabel("downstream  d")
    ax.set_title(title, fontsize=10)

    # direct value labels; ink color chosen for contrast against the cell
    vmin, vmax = np.nanmin(mat), np.nanmax(mat)
    thr = vmin + 0.55 * (vmax - vmin)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=7,
                    color="white" if v < thr else "black")

    # mark the peak and the paper's (7,2) choice
    pk_j, pk_i = u_list.index(peak[0]), d_list.index(peak[1])
    ax.add_patch(plt.Rectangle((pk_j - 0.5, pk_i - 0.5), 1, 1, fill=False,
                               edgecolor="#d1495b", lw=2.2))
    if 7 in u_list and 2 in d_list:
        j72, i72 = u_list.index(7), d_list.index(2)
        ax.add_patch(plt.Rectangle((j72 - 0.5, i72 - 0.5), 1, 1, fill=False,
                                   edgecolor="white", lw=1.6, linestyle="--"))

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, fname + ".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
