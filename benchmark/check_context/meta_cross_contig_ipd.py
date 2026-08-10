"""
Metagenome cross-contig reproducibility of per-10-mer control IPD, vs taxonomic distance
and environment.

Restricted to the high-quality (near-complete) contigs listed in
tmp/figures/multi_env_linkage/motif_num_all_samples.csv. For each metagenome (sample),
we compute the per-MODIFI-10-mer mean IPD of each HQ contig (reusing MODIFI's precomputed
ipd/<contig>.ipd1.csv), then correlate 10-mer IPD between contig pairs WITHIN THE SAME
sample. Each pair is tagged by taxonomic distance (lowest common GTDB rank) and the
sample's environment.

Two stages (control with --stage):
  cache : build & cache a per-contig 10-mer profile (code, mean log-IPD, n) -> meta_profiles/
  pairs : within-sample pairwise Pearson r from caches -> pairs.csv + two figures
  all   : both (default)

Run:  conda run -n modifi python benchmark/check_context/meta_cross_contig_ipd.py
"""

import os
import argparse
import numpy as np
import pandas as pd

import one_way_anova as owa

U, D = 7, 2                # MODIFI 10-mer
COVMIN = 10                # min per-position read coverage (HiFi metagenome ~16x median)
MIN_POS = 3                # min positions per 10-mer within a contig
MIN_SHARED = 30            # min shared 10-mers to report a pair correlation
MAXCTG = 10**9             # no contig cap: use ALL HQ contigs per sample
MAXROWS = 3_000_000        # I/O cap: read at most this many ipd1 rows per contig
                           # (~1.5 Mb of both-strand positions -> ~2e5 well-sampled 10-mers;
                           #  bounds reads on giant contigs without dropping any contig)
KPARTNER = 40              # for pairing, each contig is matched to up to this many random
                           # partners within its sample (bounds pairs to ~linear, not n^2)

RUN2 = "/home/shuaiw/borg/paper/run2"
HQ_CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/motif_num_all_samples.csv"
PROFDIR = "/home/shuaiw/borg/revision/context/meta_profiles"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/meta_cross_contig"

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def paths(sample, contig):
    base = f"{RUN2}/{sample}/{sample}_methylation4"
    return (f"{base}/ipd/{contig}.ipd1.csv", f"{base}/contigs/{contig}.fa")


def esti_relation(lin1, lin2):
    """Lowest common GTDB rank between two lineages (mirrors analyze_motif_share)."""
    if "Unclassified" in lin1 or "Unclassified" in lin2:
        return "different_lineage"
    t1, t2 = lin1.split(";"), lin2.split(";")
    for i in range(len(RANKS) - 1, -1, -1):
        if i < len(t1) and i < len(t2) and t1[i] == t2[i] and len(t1[i]) > 3:
            return f"same_{RANKS[i]}"
    return "different_lineage"


# ---------------- stage 1: per-contig 10-mer profile cache --------------------
def build_profile(task):
    sample, contig = task
    ipd, fa = paths(sample, contig)
    cache = os.path.join(PROFDIR, f"{contig}.npz")
    if os.path.exists(cache):
        return (contig, -1, -1)
    try:
        df = pd.read_csv(ipd, usecols=["strand", "tpl", "coverage", "tMean"],
                         dtype={"strand": np.int8, "tpl": np.int64,
                                "coverage": np.int32, "tMean": np.float64},
                         nrows=MAXROWS)   # I/O cap: representative prefix (both strands)
        df = df[(df["coverage"] >= COVMIN) & (df["tMean"] > 0)]
        _, ref, comp = owa.load_reference(fa)
        y_all = np.log(df["tMean"].values)
        fwd = df["strand"].values == 1
        pos_f = df["tpl"].values[fwd].astype(np.int64)
        y_f = y_all[fwd]
        pos_r = df["tpl"].values[~fwd].astype(np.int64)
        y_r = y_all[~fwd]
        pf, yf, pr, yr = owa.valid_sites(ref, comp, pos_f, y_f, pos_r, y_r, U, D)
        codes = owa.codes_at(ref, comp, pf, pr, U, D)
        y = np.concatenate([yf, yr])
        uniq, inv = np.unique(codes, return_inverse=True)
        n = np.bincount(inv, minlength=uniq.size)
        s = np.bincount(inv, weights=y, minlength=uniq.size)
        mean = s / n
        np.savez(cache, code=uniq, mean=mean.astype(np.float32), n=n.astype(np.int32))
        return (contig, int(df.shape[0]), int(uniq.size))
    except Exception as e:  # noqa
        return (contig, -2, str(e)[:80])


def stage_cache(hq, threads):
    os.makedirs(PROFDIR, exist_ok=True)
    tasks = []
    for sample, g in hq.groupby("sample"):
        if len(g) < 2:
            continue
        top = g.sort_values("ctg_len", ascending=False).head(MAXCTG)
        tasks += [(sample, c) for c in top["contig"]]
    todo = [t for t in tasks if not os.path.exists(os.path.join(PROFDIR, f"{t[1]}.npz"))]
    print(f"stage cache: {len(tasks)} HQ contigs selected, {len(todo)} to build "
          f"({len(tasks) - len(todo)} cached)", flush=True)
    from concurrent.futures import ProcessPoolExecutor
    done = 0
    with ProcessPoolExecutor(max_workers=threads) as ex:
        for contig, npos, ncodes in ex.map(build_profile, todo):
            done += 1
            if npos == -2:
                print(f"  ERROR {contig}: {ncodes}", flush=True)
            elif done % 25 == 0 or done == len(todo):
                print(f"  {done}/{len(todo)} built (last {contig}: "
                      f"{npos} pos, {ncodes} 10-mers)", flush=True)
    return tasks


# ---------------- stage 2: within-sample pairwise correlation -----------------
def load_arr(contig):
    """Return (code, mean) numpy arrays for 10-mers with >= MIN_POS positions.
    code is sorted-unique (np.unique output), enabling fast set intersection."""
    cache = os.path.join(PROFDIR, f"{contig}.npz")
    if not os.path.exists(cache):
        return None
    d = np.load(cache)
    keep = d["n"] >= MIN_POS
    return d["code"][keep], d["mean"][keep].astype(np.float64)


def zscore(x):
    return (x - x.mean()) / x.std(ddof=0)


def stage_pairs(hq):
    os.makedirs(OUTDIR, exist_ok=True)
    lineage = dict(zip(hq["contig"], hq["lineage"]))
    env = dict(zip(hq["contig"], hq["environment"]))
    rng = np.random.default_rng(12345)
    rows = []
    for sample, g in hq.groupby("sample"):
        if len(g) < 2:
            continue
        contigs = list(g.sort_values("ctg_len", ascending=False).head(MAXCTG)["contig"])
        arrs = {}
        for c in contigs:
            a = load_arr(c)
            if a is not None and a[0].size > 0:
                arrs[c] = a
        cs = list(arrs)
        if len(cs) < 2:
            continue
        # pairs: each contig matched to up to KPARTNER random partners in the sample
        pairset = set()
        for i in range(len(cs)):
            others = [j for j in range(len(cs)) if j != i]
            k = min(KPARTNER, len(others))
            for j in rng.choice(others, size=k, replace=False):
                pairset.add((min(i, int(j)), max(i, int(j))))
        for i, j in pairset:
            a, b = cs[i], cs[j]
            ca, ma = arrs[a]
            cb, mb = arrs[b]
            common, ia, ib = np.intersect1d(ca, cb, assume_unique=True,
                                            return_indices=True)
            if common.size < MIN_SHARED:
                continue
            x, y = ma[ia], mb[ib]
            if x.std() == 0 or y.std() == 0:
                continue
            r = float(np.corrcoef(x, y)[0, 1])   # Pearson (scale/shift invariant)
            rows.append(dict(sample=sample, environment=env[a], contigA=a, contigB=b,
                             tax_rank=esti_relation(lineage[a], lineage[b]),
                             pearson_r=r, n_shared=int(common.size)))
        print(f"  {sample}: {len(cs)} contigs, {len(pairset)} pairs", flush=True)

    pairs = pd.DataFrame(rows)
    pairs.to_csv(os.path.join(OUTDIR, "pairs.csv"), index=False)
    print(f"\ntotal pairs: {len(pairs)}  median r={pairs['pearson_r'].median():.3f}  "
          f"median n_shared={pairs['n_shared'].median():.0f}", flush=True)
    summarize(pairs)
    plot(pairs)
    print("wrote outputs to", OUTDIR, flush=True)


def summarize(pairs):
    rank_order = ["same_species", "same_genus", "same_family", "same_order",
                  "same_class", "same_phylum", "same_domain"]
    by_rank = (pairs.groupby("tax_rank")["pearson_r"]
               .agg(["count", "median",
                     lambda s: s.quantile(.25), lambda s: s.quantile(.75)]))
    by_rank.columns = ["n_pairs", "median_r", "q25", "q75"]
    by_rank = by_rank.reindex([r for r in rank_order if r in by_rank.index])
    by_rank.to_csv(os.path.join(OUTDIR, "summary_by_rank.csv"))
    by_env = (pairs.groupby("environment")["pearson_r"]
              .agg(["count", "median",
                    lambda s: s.quantile(.25), lambda s: s.quantile(.75)]))
    by_env.columns = ["n_pairs", "median_r", "q25", "q75"]
    by_env.to_csv(os.path.join(OUTDIR, "summary_by_environment.csv"))
    print("\n--- by taxonomic rank ---\n", by_rank, flush=True)
    print("\n--- by environment ---\n", by_env, flush=True)


def plot(pairs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def box(ax, groups, labels, title, xlabel, ticklabels=None):
        data = [groups[g] for g in labels]
        bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6)
        for p in bp["boxes"]:
            p.set(facecolor="#9ecae1", edgecolor="#333")
        for med in bp["medians"]:
            med.set(color="#d1495b", lw=2)
        ax.set_xticks(range(1, len(labels) + 1))
        ax.set_xticklabels(ticklabels or labels, rotation=30, ha="right", fontsize=13)
        for k, g in enumerate(labels, 1):
            ax.text(k, 0.02, f"n={len(groups[g])}", ha="center", va="bottom",
                    fontsize=11, color="#555", transform=ax.get_xaxis_transform())
        ax.set_ylim(0, 1)
        ax.tick_params(axis="y", labelsize=13)
        ax.set_ylabel("Pearson r (per-10-mer IPD, contig pair)", fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.grid(axis="y", ls=":", alpha=0.5)

    fig, (ax2, ax1) = plt.subplots(1, 2, figsize=(18, 5.5))

    # Panel a: r vs taxonomic distance
    rank_order = ["same_species", "same_genus", "same_family", "same_order",
                  "same_class", "same_phylum", "same_domain"]
    # lowest common rank X means the pair differs at the rank just below X
    differs = {"same_genus": "different_species", "same_family": "different_genus",
               "same_order": "different_family", "same_class": "different_order",
               "same_phylum": "different_class", "same_domain": "different_phylum"}
    grp = {r: pairs.loc[pairs["tax_rank"] == r, "pearson_r"].values for r in rank_order}
    labs = [r for r in rank_order if len(grp[r]) > 0]
    ticklabs = [f"{r}\n({differs[r]})" if r in differs else r for r in labs]
    box(ax1, grp, labs, "", "lowest common taxonomic rank", ticklabels=ticklabs)

    # Panel b: r across environments (ordered by median)
    envs = (pairs.groupby("environment")["pearson_r"].median()
            .sort_values(ascending=False).index.tolist())
    genv = {e: pairs.loc[pairs["environment"] == e, "pearson_r"].values for e in envs}
    box(ax2, genv, envs, "", "habitat")

    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "cross_contig_ipd.pdf"))
    plt.close(fig)


def main():
    global MAXCTG
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["cache", "pairs", "plot", "all"], default="all")
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--maxctg", type=int, default=MAXCTG)
    ap.add_argument("--env", default=None,
                    help="comma-separated environment(s) to restrict to, e.g. 'soil'")
    args = ap.parse_args()
    MAXCTG = args.maxctg

    hq = pd.read_csv(HQ_CSV)
    hq = hq[hq["lineage"].str.startswith("d__", na=False)].copy()  # classified only
    if args.env:
        envs = set(args.env.split(","))
        hq = hq[hq["environment"].isin(envs)].copy()
    print(f"HQ contigs (classified): {len(hq)} across {hq['sample'].nunique()} samples, "
          f"{hq['environment'].nunique()} environments"
          + (f"  [env filter: {args.env}]" if args.env else ""), flush=True)

    if args.stage == "plot":
        # quick reproduce: rebuild figures/summaries from the saved intermediate table
        pairs = pd.read_csv(os.path.join(OUTDIR, "pairs.csv"))
        print(f"loaded {len(pairs)} pairs from pairs.csv", flush=True)
        summarize(pairs)
        plot(pairs)
        print("wrote outputs to", OUTDIR, flush=True)
        return
    if args.stage in ("cache", "all"):
        stage_cache(hq, args.threads)
    if args.stage in ("pairs", "all"):
        stage_pairs(hq)


if __name__ == "__main__":
    main()
