#!/usr/bin/env python3
"""build_strain_mix.py — clean C2 rebuild with a DISJOINT donor/background species split.

Unlike the old C2 (bg300_k*, whose background leaked donor-species con-specific strains),
strain_mix's background is drawn ONLY from non-donor species (bc.select_background with
exclude_donor=True), so K = max con-specific strains per DONOR species is the sole
con-specific variable and K=1 is truly 1 strain/species -> strain accuracy should be ~100%.

Self-contained & seed-deterministic: all 58 donor species are used; per species the strain
order is shuffled once by seed and the first-K strains are taken, so levels are nested
(k1 subset k2 ... kall) at a given seed. Background is fixed across K at a given seed.

Usage: build_strain_mix.py --K 2 [--tag rep2 --seed 43 --n-background 242]
Output label: strain_mix_k<K>[_<tag>]  under simu_meta_dir/C1/
"""
import argparse
import pandas as pd
import build_community as bc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--K", required=True, help="1, 2, 3, 4, or all")
    ap.add_argument("--tag", default="", help="replicate suffix, e.g. rep2")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-background", type=int, default=242)
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--keep-prepped", action="store_true")
    a = ap.parse_args()

    donors, background, mge = bc.load_pool()
    # Keep all donor species, INCLUDING the blank-GTDB-name group (6 isolates without a
    # species-level assignment) so K=1 = 58 donor species, matching old C2. groupby fills
    # NaN->"" so those 6 form one "" species group (they resolve fine, being distinct orgs).
    donors = donors.copy()
    donors["Species"] = donors["Species"].fillna("")

    # all donor species; per species: 1 isolate/strain, seed-shuffled, take first-K (nested)
    picked = []
    for sp, grp in donors.groupby("Species"):
        reps = (grp.sort_values("Sample").groupby("Strain", as_index=False).first()
                .sample(frac=1, random_state=a.seed))
        k = len(reps) if a.K == "all" else min(int(a.K), len(reps))
        picked.append(reps.head(k))
    donor_rows = pd.concat(picked, ignore_index=True)
    donor_species = set(donor_rows["Species"].dropna())

    # clean background: non-donor species only
    bg_rows = bc.select_background(background, a.n_background, donor_species, a.seed,
                                   exclude_donor=True)

    label = f"strain_mix_k{a.K}" + (f"_{a.tag}" if a.tag else "")
    print(f"[{label}] donor strains={len(donor_rows)} ({len(donor_species)} species) + "
          f"background={len(bg_rows)} ({bg_rows['Species'].nunique() if len(bg_rows) else 0} "
          f"non-donor species) => {len(donor_rows)+len(bg_rows)} genomes")
    # sanity: background must share NO species with donors
    overlap = set(bg_rows["Species"]) & donor_species if len(bg_rows) else set()
    assert not overlap, f"background leaked donor species: {overlap}"
    bc.build_one(label, donor_rows, bg_rows, mge, seed=a.seed,
                 threads=a.threads, keep_prepped=a.keep_prepped)


if __name__ == "__main__":
    main()
