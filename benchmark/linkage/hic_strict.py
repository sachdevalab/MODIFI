#!/usr/bin/env python3
"""Hi-C validation of MODIFI ECE-host linkages, restricted to the STRICT high-confidence ECE dataset.
Reuses compare_bin3c functions. Runs on the Hi-C samples (cow_bioreactor_1/4/5), keeps only linkages whose
(sample, ECE) is in the final strict linkage table, and compares their bin3c contact values to random pairs.
Output -> tmp/rev_figs/ece_anno/hic/ ; linkage category = "MODIFI linkages supported by Hi-C", hue = ECE type."""
import sys, os, csv
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/linkage")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from scipy.stats import ttest_ind
import compare_bin3c as cb
from compare_bin3c import assess_linage, read_contact_values, compare_contact, random_test
from sample_object import My_sample

ALL_DIR = "/home/shuaiw/borg/paper/run2/"
LK = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/hic"
os.makedirs(OUT, exist_ok=True)
LABEL = "MODIFI linkages supported by Hi-C"
HIC_SAMPLES = ["cow_bioreactor_1", "cow_bioreactor_4", "cow_bioreactor_5"]

strict = {(r["sample"], r["MGE"]): r["type"]
          for r in csv.DictReader(open(LK))}   # (sample, ECE) -> plasmid/virus

rows = []
for prefix in HIC_SAMPLES:
    so = My_sample(prefix, ALL_DIR, methy_v=4)
    our = assess_linage(so.read_bin3c(), so)
    cv = read_contact_values(so.contact_value_file)
    n = 0
    for mge in our:
        if (prefix, mge) not in strict:
            continue                              # strict-set filter
        bin_name = our[mge][0]
        _, _, _, _, skip = compare_contact(mge, bin_name, cv)
        if skip:
            continue
        rows.append([prefix, LABEL, cv[mge][bin_name], mge, strict[(prefix, mge)]]); n += 1
    cb.contact_values = cv                        # random_test/get_info_ctg read this module global
    for v in random_test(cv, iterations=1000):
        rows.append([prefix, "Random Pairs", v, "", ""])
    print(f"{prefix}: strict linkages tested = {n}")

df = pd.DataFrame(rows, columns=["Sample", "ECE type", "Contact Value", "ECE", "MGE_type"])
df.to_csv(f"{OUT}/contact_values.csv", index=False)

ol = df.loc[df["ECE type"] == LABEL, "Contact Value"].values
rp = df.loc[df["ECE type"] == "Random Pairs", "Contact Value"].values
t, p = ttest_ind(ol, rp, equal_var=False)
print(f"\nSTRICT Hi-C: linkages n={len(ol)}  mean={np.mean(ol):.1f} median={np.median(ol):.1f} | "
      f"random mean={np.mean(rp):.2f}  Welch t={t:.3f} p={p:.2e}")
print(f"strict linkages with positive Hi-C contact: {(ol > 0).sum()}/{len(ol)} "
      f"({100*(ol > 0).mean():.0f}%)")

plt.figure(figsize=(5, 5))
sns.boxplot(x="Sample", y="Contact Value", hue="ECE type", data=df)
plt.yscale("log"); plt.ylabel("Contact Value (log scale)"); plt.xlabel("")
plt.xticks(rotation=30, ha="right"); plt.legend(loc="best", fontsize=8, title="ECE type")
plt.tight_layout()
plt.savefig(f"{OUT}/bin3c_boxplot.pdf"); plt.savefig(f"{OUT}/bin3c_boxplot.png", dpi=200)
print("wrote", OUT, "/bin3c_boxplot.{pdf,png} + contact_values.csv")
