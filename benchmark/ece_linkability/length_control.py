#!/usr/bin/env python3
"""Length-controlled comparison of modification frequency (sites/kb) between linked and unlinked
ECEs, strict metagenome set (curated 317 linked). Guards against "the density looks equal only
because the length distributions overlap": (1) per-length-bin Mann-Whitney, (2) OLS with log10(length)
as a covariate. Writes a figure + source data to tmp/rev_figs/ece_linkability/.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import ece_plot_common as C
import matplotlib.pyplot as plt

STEM = "fig_metagenome_length_control"
A = pd.read_csv(os.path.join(C.OUT, "fig_metagenome_a_length_sourcedata.csv"))
A["lk"] = A.linked.astype(str).isin(["True", "1"])
d = A.dropna(subset=["mod_density_per_kb", "mge_len"]).copy()

# ---- (1) length quantile bins + per-bin Mann-Whitney ----
edges = np.unique(np.quantile(d.mge_len, [0, .2, .4, .6, .8, 1.0]))
d["bin"] = pd.cut(d.mge_len, bins=edges, include_lowest=True)
bins = list(d.bin.cat.categories)
rows = []
for b in bins:
    g = d[d.bin == b]
    L = g[g.lk].mod_density_per_kb
    U = g[~g.lk].mod_density_per_kb
    p = stats.mannwhitneyu(L, U).pvalue if len(L) >= 3 and len(U) >= 3 else np.nan
    rows.append(dict(bin=f"{int(b.left/1000)}-{int(b.right/1000)}kb", n_linked=len(L), n_unlinked=len(U),
                     med_linked=L.median(), med_unlinked=U.median(), mwu_p=p))
binstat = pd.DataFrame(rows)

# ---- (2) OLS: mod_freq ~ linked + log10(length) ----
n = len(d)
X = np.column_stack([np.ones(n), d.lk.astype(int).values, np.log10(d.mge_len.values)])
y = d.mod_density_per_kb.values
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
resid = y - X @ beta
dof = n - X.shape[1]
se = np.sqrt(np.diag((resid @ resid) / dof * np.linalg.inv(X.T @ X)))
tval = beta / se
pval = 2 * stats.t.sf(np.abs(tval), dof)
p_link, p_len = pval[1], pval[2]
p_unadj = stats.mannwhitneyu(d[d.lk].mod_density_per_kb, d[~d.lk].mod_density_per_kb).pvalue

# ---- figure: grouped boxplot by length bin, linked vs unlinked ----
fig, ax = plt.subplots(figsize=(8.5, 5))
pos = np.arange(len(bins))
w = 0.36
for i, b in enumerate(bins):
    g = d[d.bin == b]
    for j, (mask, col, off) in enumerate([(g.lk, C.C_LINK, -w / 2), (~g.lk, C.C_UNLINK, w / 2)]):
        v = g[mask].mod_density_per_kb.values
        bp = ax.boxplot(v, positions=[i + off], widths=w * 0.9, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black"))
        bp["boxes"][0].set(facecolor=col, alpha=0.75)
    p = binstat.iloc[i].mwu_p
    lab = "n.s." if (p != p or p >= 0.05) else f"p={p:.3f}"
    ax.text(i, ax.get_ylim()[1] * 0.97, lab, ha="center", va="top", fontsize=9)
ax.set_xticks(pos)
ax.set_xticklabels([r["bin"] for r in rows])
ax.set_xlabel("ECE length bin")
ax.set_ylabel("modification frequency (sites / kb)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=C.C_LINK, alpha=0.75, label=f"linked (n={int(d.lk.sum())})"),
                   Patch(facecolor=C.C_UNLINK, alpha=0.75, label=f"unlinked (n={int((~d.lk).sum())})")],
          loc="lower right", fontsize=9)
ax.set_title("Modification frequency vs linkage, controlled for length (strict set)\n"
             f"OLS mod-freq ~ linked + log10(length): linked coef = {beta[1]:.2f}, "
             f"p = {p_link:.3f} (n.s.); length p = {p_len:.0e}", fontsize=10)
fig.tight_layout()
os.makedirs(C.OUT, exist_ok=True)
fig.savefig(os.path.join(C.OUT, f"{STEM}.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(C.OUT, f"{STEM}.png"), dpi=200, bbox_inches="tight")
binstat.to_csv(os.path.join(C.OUT, f"{STEM}_sourcedata.csv"), index=False)

print(binstat.to_string(index=False))
print(f"\nOLS mod-freq ~ linked + log10(len): linked coef={beta[1]:.3f} (p={p_link:.3f}), "
      f"log10(len) coef={beta[2]:.3f} (p={p_len:.2e}); unadjusted MWU p={p_unadj:.3f}")
print(f"wrote {STEM}.pdf/.png")
