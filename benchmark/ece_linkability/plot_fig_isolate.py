#!/usr/bin/env python3
"""ECE linkage-ability -- ISOLATE dataset: data + stats producer (plasmids and viruses POOLED,
since viruses are too few to analyse independently). Writes the Source Data CSVs read by
plot_fig_isolate.R; renders no figure.

Evaluable ECE universe = full gate (depth >= 10x AND length >= 5 kb) with a usable motif profile
(1,120 ECEs: 993 plasmids + 127 viruses). linked = confident (final_score>0.5 & specificity<0.01)
AND correct host. Motif density = host chromosome recognition-site occurrences/kb.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import ece_plot_common as C

DATASET = "isolate"
STEM = "fig_isolate_linkability"


def stars(p):
    if p != p:
        return "ns"
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"


def main():
    df = C.load_master()
    d = df[df.dataset == DATASET].copy()
    ev = d[d.type.isin(["plasmid", "virus"]) & d.passes_full_gate
           & d.motif_density_per_kb.notna()].copy()
    print(f"evaluable ECEs: {len(ev)}  " + str(ev.type.value_counts().to_dict()))

    lk = ev[ev.linked]
    ul = ev[~ev.linked]
    host = ev.drop_duplicates("sample")["host_own_density_per_kb"].dropna()  # one per isolate

    # ---- Source Data (pooled plasmid + virus; type kept for reference) ----
    length_df = pd.concat([
        pd.DataFrame({"group": "linked", "type": lk.type.values, "ece_len_bp": lk.ece_len.values}),
        pd.DataFrame({"group": "unlinked", "type": ul.type.values, "ece_len_bp": ul.ece_len.values}),
    ], ignore_index=True)
    density_df = pd.concat([
        pd.DataFrame({"group": "host_chromosome", "density_per_kb": host.values}),
        pd.DataFrame({"group": "linked_ECE", "density_per_kb": lk.motif_density_per_kb.values}),
        pd.DataFrame({"group": "unlinked_ECE", "density_per_kb": ul.motif_density_per_kb.values}),
    ], ignore_index=True)
    # modification (GFF modified_base, score cutoff 30) grouped by linkage score, NOT by the
    # confident linked/unlinked call: score>0 (any host-motif modification signal) vs score=0 (none)
    pos = ev[ev.final_score > 0]
    zero = ev[~(ev.final_score > 0)]   # score == 0 or missing = no linkage signal
    mod_df = pd.concat([
        pd.DataFrame({"group": "score>0", "n_mod_sites": pos.n_mod_sites.values,
                      "mod_density_per_kb": pos.mod_density_per_kb.values}),
        pd.DataFrame({"group": "score=0", "n_mod_sites": zero.n_mod_sites.values,
                      "mod_density_per_kb": zero.mod_density_per_kb.values}),
    ], ignore_index=True)
    os.makedirs(C.OUT, exist_ok=True)
    length_df.to_csv(os.path.join(C.OUT, f"{STEM}_a_length_sourcedata.csv"), index=False)
    density_df.to_csv(os.path.join(C.OUT, f"{STEM}_b_density_sourcedata.csv"), index=False)
    mod_df.to_csv(os.path.join(C.OUT, f"{STEM}_cd_modification_sourcedata.csv"), index=False)

    # ---- stats ----
    p_len = stats.ttest_ind(np.log10(lk.ece_len), np.log10(ul.ece_len), equal_var=False).pvalue
    p_hl = stats.ttest_ind(host, lk.motif_density_per_kb, equal_var=False).pvalue
    p_lu = stats.ttest_ind(lk.motif_density_per_kb, ul.motif_density_per_kb, equal_var=False).pvalue
    p_hu = stats.ttest_ind(host, ul.motif_density_per_kb, equal_var=False).pvalue
    # modification grouped by linkage score: log1p for count (heavy-tailed), raw for density
    p_modn = stats.ttest_ind(np.log1p(pos.n_mod_sites), np.log1p(zero.n_mod_sites), equal_var=False).pvalue
    p_modd = stats.ttest_ind(pos.mod_density_per_kb, zero.mod_density_per_kb, equal_var=False).pvalue
    print(f"pooled: linked={len(lk)} unlinked={len(ul)} isolates(host)={len(host)}")
    print(f" length median linked {int(lk.ece_len.median()):,} vs unlinked {int(ul.ece_len.median()):,}  "
          f"{stars(p_len)} p={p_len:.2e}")
    print(f" recog density host {host.median():.2f} / linked {lk.motif_density_per_kb.median():.2f} / "
          f"unlinked {ul.motif_density_per_kb.median():.2f}  h-l {stars(p_hl)} l-u {stars(p_lu)} h-u {stars(p_hu)}")
    print(f" [by score] score>0 (n={len(pos)}) vs score=0 (n={len(zero)}):")
    print(f"   mod SITES  {pos.n_mod_sites.median():.0f} vs {zero.n_mod_sites.median():.0f}  "
          f"{stars(p_modn)} p={p_modn:.2e}")
    print(f"   mod DENSITY {pos.mod_density_per_kb.median():.2f} vs {zero.mod_density_per_kb.median():.2f}/kb  "
          f"{stars(p_modd)} p={p_modd:.2e}")
    print(f"wrote Source Data to {C.OUT}")


if __name__ == "__main__":
    main()
