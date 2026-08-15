#!/usr/bin/env python3
"""build_c3_cov.py — build one C3 coverage-titration community (fixed composition, varying
donor depth). Reproduces the bg_80 composition (40 donor species x1 strain + 40 background,
seed 42) at a given DONOR depth so recall/precision can be mapped vs coverage. Background
depth stays at the usual ~10x. bg_80 (30x) is the high-depth anchor already built.

Usage: build_c3_cov.py --depth 5   ->  community label cov_d5 under simu_meta_dir/C1/
"""
import argparse
import build_community as bc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depth", type=int, required=True, help="donor target depth (5/10/20/40)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=32)
    a = ap.parse_args()
    bc.DEPTH_DONOR = a.depth                          # monkeypatch the donor target depth
    print(f"[cov_d{a.depth}] donor depth={a.depth}x (bg={bc.DEPTH_BG}x), bg_80 composition seed {a.seed}")
    bc.build_community(n_species=40, strains_per_species=1, n_background=40,
                       label=f"cov_d{a.depth}", seed=a.seed, threads=a.threads,
                       keep_prepped=False)


if __name__ == "__main__":
    main()
