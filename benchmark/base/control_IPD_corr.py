import pandas as pd
from scipy.stats import spearmanr

native_csv = "/home/shuaiw/borg/paper/base/pure/native/ipd_ratio/CP011331.1.ipd3.csv"
control_csv = "/home/shuaiw/borg/paper/base/pure/control/ipd_ratio/CP011331.1.ipd3.csv"

meta_native_csv = "/home/shuaiw/borg/paper/base/meta/native/ipd_ratio/CP011331.1.ipd3.csv"
meta_control_csv = "/home/shuaiw/borg/paper/base/meta/control/ipd_ratio/CP011331.1.ipd3.csv"


DTYPES = {"tpl": "int32", "strand": "int8", "tMean": "float32", "coverage": "int32", "control": "float32"}
DEPTH_MIN = 20


def load_ipd(path, extra_cols=None):
    cols = ["tpl", "strand", "tMean", "coverage"] + (extra_cols or [])
    dtypes = {k: DTYPES[k] for k in cols if k in DTYPES}
    df = pd.read_csv(path, usecols=cols, dtype=dtypes, engine="c")
    return df[df["coverage"] > DEPTH_MIN].drop(columns="coverage")


def compute_corr(s_a, s_b, label_a, label_b):
    merged = s_a.merge(s_b, on=["tpl", "strand"], suffixes=("_a", "_b"))
    r, p = spearmanr(merged["tMean_a"], merged["tMean_b"])
    print(f"{label_a} vs {label_b}: n={len(merged)}, r={r:.4f}, p={p:.2e}")
    return r


# Correlation: native vs control (pure reads, ground truth)
control = load_ipd(control_csv)
native = load_ipd(native_csv)
compute_corr(native, control, "native (pure)", "control (ground truth)")

# Correlation: MODIFI estimated control (meta_native["control"]) vs ground truth control tMean
meta_native = load_ipd(meta_native_csv, extra_cols=["control"])
meta_est_ctrl = meta_native[["tpl", "strand", "control"]].rename(columns={"control": "tMean"})
compute_corr(meta_est_ctrl, control, "meta estimated control (MODIFI)", "control (ground truth)")
