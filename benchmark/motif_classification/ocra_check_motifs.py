#!/usr/bin/env python3
"""
Check whether MODIFI recovered the known Oceanispirochaeta crateris K2 methylome
(Fomenkov et al. 2020) in each run (merged + per-cell):

  GCNGC  -> 5mC   (modified C at position 2)   <- the target
  GGATCC -> 4mC   (modified C at position 5)
  ATTAAT -> 6mA   (modified A at position 5)

For each run's all.motifs.csv, IUPAC- + revcomp-aware match each expected recognition sequence and
report: detected?, the detected motifString, its centerPos + modified base, fraction, nDetected,
nGenome, and the run's mean depth. Writes one row per (run x expected motif).
"""
import os, re, sys
import pandas as pd
from Bio.Seq import Seq

ROOT = "/home/shuaiw/borg/revision/ocra_5mC"
FIGDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"
os.makedirs(FIGDIR, exist_ok=True)

RUNS = [("merged", "modifi_merged"),
        ("SRR8580034", "modifi_SRR8580034"),
        ("SRR8580035", "modifi_SRR8580035")]

# expected recognition, known type, and the modified-base position (1-based) per the paper
EXPECT = [("GCNGC", "5mC", 2), ("GGATCC", "4mC", 5), ("ATTAAT", "6mA", 5)]

IU = {'A':'A','C':'C','G':'G','T':'T','U':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]',
      'K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
def rx(m): return re.compile('^' + ''.join(IU.get(b.upper(), b) for b in m) + '$')
def rc(m): return str(Seq(m).reverse_complement())

def same(det, exp):
    """det (detected) matches exp (expected) if equal length and IUPAC-compatible on either strand."""
    if len(det) != len(exp):
        return False
    for a, b in ((det, exp), (det, rc(exp))):
        if rx(a).match(b) or rx(b).match(a):
            return True
    return False

def mean_depth(rundir):
    f = os.path.join(ROOT, rundir, "mean_depth.csv")
    if not os.path.exists(f): return None
    d = pd.read_csv(f)
    # single-chromosome genome; report its depth
    return round(float(d["depth"].max()), 1) if "depth" in d.columns and len(d) else None

rows = []
for run_label, rundir in RUNS:
    mf = os.path.join(ROOT, rundir, "all.motifs.csv")
    dep = mean_depth(rundir)
    if not os.path.exists(mf):
        for exp, typ, pos in EXPECT:
            rows.append(dict(run=run_label, expected=exp, known_type=typ, mean_depth=dep,
                             detected="run-missing", detected_motif="", centerPos="",
                             modified_base="", fraction="", nDetected="", nGenome=""))
        continue
    df = pd.read_csv(mf)
    for exp, typ, pos in EXPECT:
        hit = None
        for _, r in df.iterrows():
            if same(str(r.motifString), exp):
                hit = r; break
        if hit is None:
            rows.append(dict(run=run_label, expected=exp, known_type=typ, mean_depth=dep,
                             detected="NO", detected_motif="", centerPos="",
                             modified_base="", fraction="", nDetected="", nGenome=""))
        else:
            m = str(hit.motifString); cp = int(hit.centerPos)
            base = m[cp-1] if 1 <= cp <= len(m) else "?"
            rows.append(dict(run=run_label, expected=exp, known_type=typ, mean_depth=dep,
                             detected="YES", detected_motif=m, centerPos=cp,
                             modified_base=base, fraction=round(float(hit.fraction), 3),
                             nDetected=int(hit.nDetected),
                             nGenome=int(hit.nGenome) if "nGenome" in hit else ""))

out = pd.DataFrame(rows, columns=["run","expected","known_type","mean_depth","detected",
                                  "detected_motif","centerPos","modified_base","fraction",
                                  "nDetected","nGenome"])
outpath = os.path.join(FIGDIR, "ocra_motif_recovery.csv")
out.to_csv(outpath, index=False)
print(out.to_string(index=False))
print(f"\nwrote {outpath}")
