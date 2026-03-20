#!/usr/bin/env python3
"""
Per sample: proportion of non-novel MGEs (plasmid + virus) that have at least one
host linkage in mge_host_gc_cov.csv. Novel MGEs are excluded from numerator and
denominator. Uses the same depth >= 5 filter and sample list as count_mge_per_sample.py.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "isolation"))
from sample_object import My_sample

DEFAULT_ALL_DIR = "/home/shuaiw/borg/paper/run2/"
DEFAULT_META_FILE = "/home/shuaiw/MODIFI/assembly_pipe/prefix_table.tab"
DEFAULT_LINKAGE_CSV = (
    "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
)


def read_metadata(meta_file):
    sample_env_dict = {}
    with open(meta_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            sample = parts[1]
            env = parts[3]
            sample_env_dict[sample] = env
    return sample_env_dict


def linked_mge_by_sample(linkage_csv):
    """sample -> set of MGE contig IDs with a host link; novel rows dropped."""
    import pandas as pd

    df = pd.read_csv(linkage_csv)
    if "MGE_type" not in df.columns or "sample" not in df.columns or "MGE" not in df.columns:
        raise ValueError(f"Unexpected columns in {linkage_csv}: {list(df.columns)}")
    df = df[df["MGE_type"].astype(str).str.lower() != "novel"]
    df = df[df["MGE_type"].isin(["plasmid", "virus"])]
    out = {}
    for sample, g in df.groupby("sample", sort=False):
        out[str(sample)] = set(g["MGE"].astype(str))
    return out


def main(
    all_dir=None,
    meta_file=None,
    linkage_csv=None,
    out_csv=None,
):
    if all_dir is None:
        all_dir = DEFAULT_ALL_DIR
    if meta_file is None:
        meta_file = DEFAULT_META_FILE
    if linkage_csv is None:
        linkage_csv = DEFAULT_LINKAGE_CSV
    if out_csv is None:
        out_csv = os.path.join(
            os.path.dirname(linkage_csv),
            "mge_linked_fraction_per_sample.csv",
        )

    sample_env_dict = read_metadata(meta_file)
    linked_by_sample = linked_mge_by_sample(linkage_csv)

    rows = []
    for my_dir in sorted(os.listdir(all_dir)):
        prefix = my_dir
        if prefix in ("ERR5621427_sludge", "ERR5621429_sludge", "ERR5621430_sludge"):
            continue
        if prefix not in sample_env_dict:
            continue
        sample_obj = My_sample(prefix, all_dir)
        mge_dict = sample_obj.read_MGE()
        if mge_dict is None:
            continue
        sample_obj.read_depth()
        depth_dict = getattr(sample_obj, "depth_dict", None) or {}
        mge_dict = {k: v for k, v in mge_dict.items() if k in depth_dict and depth_dict[k] >= 5}
        non_novel = {k for k, t in mge_dict.items() if t in ("plasmid", "virus")}
        n_non_novel = len(non_novel)
        linked_set = linked_by_sample.get(prefix, set())
        n_linked = len(non_novel & linked_set)
        prop = (n_linked / n_non_novel) if n_non_novel > 0 else float("nan")
        rows.append({
            "sample": prefix,
            "environment": sample_env_dict[prefix],
            "n_plasmid_virus": n_non_novel,
            "n_linked": n_linked,
            "prop_linked": prop,
        })

    import pandas as pd

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} rows -> {out_csv}")
    return out_csv


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(
        description="Fraction of non-novel MGEs with a host link per sample"
    )
    p.add_argument("--all-dir", default=DEFAULT_ALL_DIR, help="Base dir with sample folders")
    p.add_argument("--meta", default=DEFAULT_META_FILE, help="prefix_table.tab path")
    p.add_argument(
        "--linkage",
        default=DEFAULT_LINKAGE_CSV,
        help="mge_host_gc_cov.csv from network linkage",
    )
    p.add_argument("-o", "--out", default=None, help="Output CSV path")
    args = p.parse_args()
    main(
        all_dir=args.all_dir,
        meta_file=args.meta,
        linkage_csv=args.linkage,
        out_csv=args.out,
    )
