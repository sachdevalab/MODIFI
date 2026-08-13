#!/usr/bin/env python3
"""
build_c2_strain.py — build one C2 strain-depth community (nested on bg_300).

C2 varies max strains/species K in {1,2,3,4,all} on the fixed bg_300 structure
(58 donor species + 246 background). NESTED by construction: the base is bg_300's
*actual* donor strains (so K=1 == bg_300 exactly), and higher K ADDS the first (K-1)
additional con-specific strains per multi-strain donor species, chosen deterministically
(sorted). Background (the same 246 isolates as bg_300) and the donor-species set are held
fixed. Thus bg_300 = K=1 ⊂ K=2 ⊂ K=3 ⊂ K=4 ⊂ K=all, and strain depth is the only variable.

Each added strain carries its own ECE (all donors have curated ECEs).

Usage:  python build_c2_strain.py --K 2 [--threads 32]
        python build_c2_strain.py --K all
Output community label: bg300_k<K>  (under simu_meta_dir/C1/)
"""
import argparse
import pandas as pd
import build_community as bc

BG300_MANIFEST = f"{bc.OUT_ROOT}/bg_300/bg_300.manifest.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--K", required=True, help="2, 3, 4, or all")
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--keep-prepped", action="store_true")
    a = ap.parse_args()

    donors, background, mge = bc.load_pool()
    man = pd.read_csv(BG300_MANIFEST)
    base_don = set(man[man["role"] == "donor"]["sample"])
    base_bg = set(man[man["role"] == "background"]["sample"])

    base_donor_rows = donors[donors["Sample"].isin(base_don)].drop_duplicates("Sample")
    bg_rows = background[background["Sample"].isin(base_bg)].drop_duplicates("Sample")
    assert len(base_donor_rows) == len(base_don), "some bg_300 donors not found in pool"
    assert len(bg_rows) == len(base_bg), "some bg_300 background not found in pool"

    # additional strains per donor species (deterministic, nested by sorted Strain)
    added = []
    for sp, grp in donors.groupby("Species"):
        if not isinstance(sp, str) or sp.strip() == "":
            base_sp = base_donor_rows[base_donor_rows["Species"] == sp]
        else:
            base_sp = base_donor_rows[base_donor_rows["Species"] == sp]
        if len(base_sp) == 0:
            continue
        base_strain = base_sp["Strain"].iloc[0]
        reps = grp.sort_values("Sample").groupby("Strain", as_index=False).first()
        others = reps[reps["Strain"] != base_strain].sort_values("Strain")
        k_add = len(others) if a.K == "all" else min(int(a.K) - 1, len(others))
        if k_add > 0:
            added.append(others.head(k_add))

    donor_rows = (pd.concat([base_donor_rows] + added, ignore_index=True)
                  if added else base_donor_rows)
    label = f"bg300_k{a.K}"
    print(f"[{label}] donors={len(donor_rows)} ({donor_rows['Species'].nunique()} species) "
          f"+ background={len(bg_rows)}  => {len(donor_rows)+len(bg_rows)} genomes")
    bc.build_one(label, donor_rows, bg_rows, mge, seed=42,
                 threads=a.threads, keep_prepped=a.keep_prepped)


if __name__ == "__main__":
    main()
