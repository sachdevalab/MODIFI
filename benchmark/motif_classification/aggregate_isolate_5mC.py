#!/usr/bin/env python3
"""
Type-agnostic 5mC sweep over the isolate MTase-linker outputs.

For every isolate under isolate_mtl/<ACC>/output/mtase_assignment_table.tsv:
  - reference MTase predictions give (motif_pred, mod_type_pred): mod_type_pred 'm' = 5mC (C5),
    'ac' = 6mA/4mC (adenine / N4-cytosine, indistinguishable at the protein level).
  - match each DETECTED motif (all.motifs.csv, fraction>=0.5 & nDetected>=100) to those motif_pred
    sequences, IUPAC- and strand-aware. Type is read ONLY from the reference + the motif itself:
      * modified base A                       -> 6mA  (adenine cannot be 5mC/4mC)
      * modified base C, matched to 'm' MTase -> 5mC  (the target of this sweep)
      * modified base C, matched to 'ac'      -> 4mC
      * modified base C, no MTase match       -> unclassified
No MicrobeMod, no chemistry prior, no assumption about any motif's type.

Writes:
  isolate_motif_types.tsv   one row per detected motif (all isolates)
  isolate_5mC_hits.tsv      the subset typed 5mC (base C matched to an m5C MTase)
  isolate_5mC_summary.txt   per-type totals + isolate/species counts
"""
import os, re, sys, glob
import pandas as pd
from Bio.Seq import Seq

MTL_ROOT = "/home/shuaiw/borg/revision/motif_class/isolate_mtl"
BATCH = "/home/shuaiw/borg/paper/isolation/batch2_results"
OUT_DIR = "/home/shuaiw/borg/revision/motif_class"
FIG_SRC = "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"
MIN_FRAC, MIN_SITES = 0.5, 100

IUPAC = {'A':'A','C':'C','G':'G','T':'T','U':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]',
         'K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}


def rx(m):
    return re.compile('^' + ''.join(IUPAC.get(b.upper(), b) for b in m) + '$')


def rc(m):
    return str(Seq(m).reverse_complement())


def same_recognition(d, p):
    """True if detected motif d and predicted motif p denote the same site (either strand),
    allowing IUPAC-degeneracy differences (one pattern fullmatches the other)."""
    if len(d) != len(p):
        # allow strand/length-equal only; different lengths are not the same site here
        return False
    for a, b in ((d, p), (d, rc(p))):
        if rx(a).match(b) or rx(b).match(a):
            return True
    return False


def load_species():
    """ACC -> ScientificName from any sra_metadata_*.csv, if available."""
    sp = {}
    for f in glob.glob("/home/shuaiw/borg/paper/isolation/aws_methylation/sra_metadata_*.csv") + \
             glob.glob("/home/shuaiw/borg/paper/isolation/**/sra_metadata_*.csv", recursive=True):
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        acc_col = next((c for c in df.columns if c.lower() in ("run", "accession", "run_accession")), None)
        sp_col = next((c for c in df.columns if c.lower() in ("scientificname", "scientific_name", "organism")), None)
        if acc_col and sp_col:
            for a, s in zip(df[acc_col], df[sp_col]):
                sp.setdefault(str(a), str(s))
    return sp


def main():
    species = load_species()
    tables = sorted(glob.glob(os.path.join(MTL_ROOT, "*", "output", "mtase_assignment_table.tsv")))
    print(f"MTase tables found: {len(tables)}", file=sys.stderr)

    rows = []
    n_iso = 0
    for tpath in tables:
        acc = tpath.split("/")[-3]
        mt = pd.read_csv(tpath, sep="\t")
        # reference-derived MTase recognition motifs with a resolved sequence
        mt = mt[mt["motif_pred"].notna() & (mt["motif_pred"].astype(str).str.len() > 0)]
        m_motifs = mt[mt["mod_type_pred"] == "m"][["motif_pred", "REbase_ID"]].values.tolist()
        ac_motifs = mt[mt["mod_type_pred"] == "ac"][["motif_pred", "REbase_ID"]].values.tolist()

        mfile = os.path.join(BATCH, acc, f"{acc}_methylation4", "all.motifs.csv")
        if not os.path.exists(mfile):
            continue
        dd = pd.read_csv(mfile)
        dd = dd[(dd.fraction >= MIN_FRAC) & (dd.nDetected >= MIN_SITES)]
        if dd.empty:
            continue
        n_iso += 1
        # revcomp-dedup detected motifs
        kept, seen = [], set()
        for _, r in dd.iterrows():
            m = str(r.motifString)
            if m in seen or rc(m) in seen:
                continue
            seen.add(m)
            kept.append(r)

        for r in kept:
            m = str(r.motifString); cp = int(r.centerPos)
            base = m[cp - 1] if 1 <= cp <= len(m) else "?"
            mtype, matched, rebid = "unclassified", "", ""
            if base == "A":
                mtype = "6mA"
            elif base == "C":
                hit_m = next(((pm, rid) for pm, rid in m_motifs if same_recognition(m, str(pm))), None)
                hit_ac = next(((pm, rid) for pm, rid in ac_motifs if same_recognition(m, str(pm))), None)
                if hit_m:
                    mtype, matched, rebid = "5mC", hit_m[0], hit_m[1]
                elif hit_ac:
                    mtype, matched, rebid = "4mC", hit_ac[0], hit_ac[1]
                else:
                    mtype = "unclassified"
            else:
                mtype = "6mA" if base == "A" else "unclassified"
            rows.append(dict(accession=acc, species=species.get(acc, ""),
                             motif=m, centerPos=cp, modified_base=base,
                             fraction=r.fraction, nDetected=r.nDetected,
                             assigned_type=mtype, matched_MTase_motif=matched, REbase_ID=rebid))

    out = pd.DataFrame(rows)
    os.makedirs(FIG_SRC, exist_ok=True)
    all_path = os.path.join(OUT_DIR, "isolate_motif_types.tsv")
    out.to_csv(all_path, sep="\t", index=False)

    hits = out[out.assigned_type == "5mC"]
    hits_path = os.path.join(OUT_DIR, "isolate_5mC_hits.tsv")
    hits.to_csv(hits_path, sep="\t", index=False)

    summ = []
    summ.append(f"isolates with typed detected motifs: {n_iso}")
    summ.append(f"MTase tables processed: {len(tables)}")
    summ.append(f"total detected motifs typed: {len(out)}")
    vc = out.assigned_type.value_counts()
    for t in ["6mA", "4mC", "5mC", "unclassified"]:
        summ.append(f"  {t}: {int(vc.get(t, 0))}")
    summ.append(f"5mC detections: {len(hits)} in {hits.accession.nunique()} isolates")
    if len(hits):
        summ.append("5mC hit species: " + ", ".join(sorted(set(hits.species.dropna())) or ["(no metadata)"]))
    text = "\n".join(summ)
    with open(os.path.join(OUT_DIR, "isolate_5mC_summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)
    print(f"\nWrote:\n  {all_path}\n  {hits_path}")


if __name__ == "__main__":
    main()
