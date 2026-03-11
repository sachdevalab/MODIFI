#!/usr/bin/env python3
"""
Count MGEs (plasmid, virus, novel) per sample using the same logic as
benchmark/specificity/profile_good_ctgs.py: My_sample(prefix, all_dir).read_MGE(),
with sample -> environment from prefix_table.tab.
Writes CSV: sample, environment, n_plasmid, n_virus, n_novel for use by count_linkages.R.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "isolation"))
from sample_object import My_sample

# Same as profile_good_ctgs.py
DEFAULT_ALL_DIR = "/home/shuaiw/borg/paper/run2/"
DEFAULT_META_FILE = "/home/shuaiw/MODIFI/assembly_pipe/prefix_table.tab"


def read_metadata(meta_file):
    """sample -> environment (same as profile_good_ctgs.read_metadata)."""
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


def main(all_dir=None, meta_file=None, out_csv=None):
    if all_dir is None:
        all_dir = DEFAULT_ALL_DIR
    if meta_file is None:
        meta_file = DEFAULT_META_FILE
    if out_csv is None:
        out_csv = os.path.join(
            os.path.dirname(__file__),
            "../../tmp/figures/multi_env_linkage/network_99/mge_counts_per_sample.csv",
        )

    sample_env_dict = read_metadata(meta_file)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)

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
        # Discard MGEs not in depth_dict; keep only those in depth_dict with depth >= 5
        mge_dict = {k: v for k, v in mge_dict.items() if k in depth_dict and depth_dict[k] >= 5}
        n_plasmid = sum(1 for t in mge_dict.values() if t == "plasmid")
        n_virus = sum(1 for t in mge_dict.values() if t == "virus")
        n_novel = sum(1 for t in mge_dict.values() if t == "novel")
        rows.append({
            "sample": prefix,
            "environment": sample_env_dict[prefix],
            "n_plasmid": n_plasmid,
            "n_virus": n_virus,
            "n_novel": n_novel,
        })

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} rows -> {out_csv}")
    return out_csv


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Count MGEs per sample for linkage normalization")
    p.add_argument("--all-dir", default=DEFAULT_ALL_DIR, help="Base dir with sample folders")
    p.add_argument("--meta", default=DEFAULT_META_FILE, help="prefix_table.tab path")
    p.add_argument("-o", "--out", default=None, help="Output CSV path")
    args = p.parse_args()
    main(all_dir=args.all_dir, meta_file=args.meta, out_csv=args.out)
