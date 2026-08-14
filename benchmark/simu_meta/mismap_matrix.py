#!/usr/bin/env python3
"""mismap_matrix.py — cross-strain read mis-mapping matrix (Part E, concern #4).

For a community that preserves each isolate's original @RG (build_community does), align
the merged reads to the concatenated reference and, for every aligned read, compare its
TRUE origin isolate (from @RG) to the isolate whose contig it mapped to (SRA prefix). A
read that maps to a different isolate is mis-mapped; among con-specific strains this mixes
methylation signal. Reported:
  - overall mis-mapping rate, split by MAPQ (<20 dropped multi-mappers = coverage loss;
    >=20 wrong-strain near-ties = IPD contamination)
  - strain x strain mis-mapping matrix among the focal (donor) strains
  - mis-mapping rate vs pairwise ANI (skani) between the focal strains

Usage: mismap_matrix.py <community_dir>   (e.g. .../C1/ecoli_panel)
"""
import os, sys, subprocess, collections
import numpy as np
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/simu_meta")
import build_community as bc

PBMM2 = os.path.join(bc.MODIFI_ENV_BIN, "pbmm2")
SKANI = "/shared/software/bin/skani"
OUTFIG = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/strain_het"
T = "32"


def rg_to_sample(man):
    """Map each isolate's original @RG ID -> sample by reading its CCS BAM header."""
    rg2s = {}
    for s in man["sample"]:
        try:
            b = bc.ccs_bam_path(s)
        except FileNotFoundError:
            continue
        hdr = subprocess.run(["samtools", "view", "-H", b], capture_output=True, text=True).stdout
        for line in hdr.splitlines():
            if line.startswith("@RG"):
                for f in line.split("\t"):
                    if f.startswith("ID:"):
                        rg2s[f[3:]] = s
    return rg2s


def align(D, label):
    aln = f"{D}/{label}.mismap.bam"
    if not os.path.exists(aln):
        subprocess.run([PBMM2, "align", "--preset", "CCS", "-j", T,
                        f"{D}/{label}.ref.fa", f"{D}/{label}.bam", aln, "--sort"], check=True)
        subprocess.run(["samtools", "index", "-@", T, aln], check=True)
    return aln


def main():
    D = sys.argv[1].rstrip("/")
    label = os.path.basename(D)
    man = pd.read_csv(f"{D}/{label}.manifest.csv")
    donors = man[man["role"] == "donor"]["sample"].tolist()
    strain_of = man.set_index("sample")["strain"].astype(str).to_dict()
    rg2s = rg_to_sample(man)
    print(f"[mismap] {len(rg2s)} RG->sample; {len(donors)} donor strains")
    aln = align(D, label)

    # stream aligned reads: origin (RG) vs mapped (contig SRA prefix), by MAPQ bin
    tot = collections.Counter()                       # (mapq>=20) -> total reads
    mis = collections.Counter()                       # mis-mapped reads
    pair = collections.Counter()                      # (origin, mapped) among donors, mapq>=20
    p = subprocess.Popen(["samtools", "view", aln], stdout=subprocess.PIPE, text=True, bufsize=1)
    donor_set = set(donors)
    no_rg = 0
    for line in p.stdout:
        f = line.split("\t", 5)                        # QNAME FLAG RNAME POS MAPQ | rest(with tags)
        rname, mapq = f[2], int(f[4])
        if rname == "*":
            continue
        i = f[5].find("RG:Z:")                          # RG may be any optional field
        rg = f[5][i + 5:].split("\t", 1)[0].rstrip() if i >= 0 else None
        origin = rg2s.get(rg)
        mapped = rname.split("_")[0]
        if origin is None:
            no_rg += 1
            continue
        hi = mapq >= 20
        tot[hi] += 1
        if origin != mapped:
            mis[hi] += 1
        if hi and origin in donor_set and mapped in donor_set:
            pair[(origin, mapped)] += 1
    p.wait()
    if no_rg:
        print(f"  [warn] {no_rg} reads had an @RG not matched to a sample (RG collision on merge?)")

    for hi in (True, False):
        t, m = tot[hi], mis[hi]
        tag = "MAPQ>=20 (contamination)" if hi else "MAPQ<20 (dropped multimappers)"
        print(f"  {tag}: {m}/{t} reads mis-mapped = {100*m/t:.3f}%" if t else f"  {tag}: no reads")

    # strain x strain matrix among donors (fraction of a strain's reads landing on each strain)
    ds = sorted(donor_set, key=lambda s: strain_of.get(s, s))
    idx = {s: i for i, s in enumerate(ds)}
    Mx = np.zeros((len(ds), len(ds)))
    for (o, mp), c in pair.items():
        Mx[idx[o], idx[mp]] += c
    row = Mx.sum(1, keepdims=True); frac = np.divide(Mx, row, where=row > 0)
    labels = [strain_of.get(s, s) for s in ds]
    pd.DataFrame(Mx, index=labels, columns=labels).to_csv(f"{D}/mismap_counts.csv")
    pd.DataFrame(frac, index=labels, columns=labels).round(4).to_csv(f"{D}/mismap_frac.csv")

    # pairwise ANI (skani) between donor genomes
    gmap = man.set_index("sample")["genome"].to_dict()
    glist = f"{D}/donor_genomes.txt"
    with open(glist, "w") as o:
        for s in ds:
            if isinstance(gmap.get(s), str) and os.path.exists(gmap[s]):
                o.write(gmap[s] + "\n")
    ani_out = f"{D}/donor_ani.tsv"
    subprocess.run([SKANI, "triangle", "-l", glist, "-o", ani_out, "--full-matrix", "-t", T],
                   check=True, capture_output=True)

    _plot(D, label, labels, frac, ani_out, ds, gmap)
    print(f"[mismap] wrote {D}/mismap_frac.csv and figure")


def _plot(D, label, labels, frac, ani_out, ds, gmap):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(OUTFIG, exist_ok=True)
    # ANI matrix aligned to ds order
    ani = pd.read_csv(ani_out, sep="\t", index_col=0)
    base = {os.path.basename(gmap[s]): s for s in ds if isinstance(gmap.get(s), str)}
    # map skani row/col files -> strain label
    def lab(fn):
        s = base.get(os.path.basename(fn)); return None if s is None else labels[ds.index(s)]
    ani.index = [lab(i) for i in ani.index]; ani.columns = [lab(c) for c in ani.columns]

    fig, ax = plt.subplots(1, 2, figsize=(13, 5.4))
    im0 = ax[0].imshow(frac, cmap="magma", vmin=0, vmax=max(0.05, np.nanmax(frac - np.eye(len(frac)))))
    ax[0].set(xticks=range(len(labels)), yticks=range(len(labels)),
              title="A. cross-strain read mis-mapping\n(fraction of a strain's reads -> each strain, MAPQ>=20)",
              xlabel="mapped-to strain", ylabel="true origin strain")
    ax[0].set_xticklabels(labels, rotation=90, fontsize=7); ax[0].set_yticklabels(labels, fontsize=7)
    fig.colorbar(im0, ax=ax[0], label="read fraction")
    # off-diagonal mis-map vs ANI scatter
    xs, ys = [], []
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            if i != j and a in ani.index and b in ani.columns:
                v = ani.loc[a, b]
                if pd.notna(v):
                    xs.append(float(v)); ys.append(frac[i, j] * 100)
    ax[1].scatter(xs, ys, s=30, color="#0072B2", edgecolor="k", lw=0.3)
    ax[1].set(xlabel="pairwise ANI (%)", ylabel="mis-mapping rate (% of origin reads)",
              title="B. cross-strain mis-mapping vs strain similarity")
    fig.suptitle(f"{label}: cross-strain read mis-mapping among con-specific strains", y=1.02, fontsize=12)
    fig.tight_layout()
    out = f"{OUTFIG}/{label}_mismapping.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
