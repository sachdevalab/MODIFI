#!/usr/bin/env python3
"""
JF8 WGA-control motif x contig matrix. Collapse IUPAC/fragmentation variants but KEEP
reverse-complement partners as SEPARATE rows (both strands, as in the original paper figure).
Each detected motif is assigned to one of the 16 ground-truth families AND a strand, by a
position-wise IUPAC-compatibility score against the family rep and its reverse complement.
Row value per contig = max profile over same-(family,strand) motifs.

Writes: jf8_wga_family_matrix.csv, jf8_wga_boldpos.csv (in borg/revision/motif_class)
"""
import pandas as pd
from Bio.Seq import Seq

D = "/home/shuaiw/borg/revision/motif_class"
J = f"{D}/jf8_wga"

# Strand-specific rows in display order: symmetric/palindromic families = 1 row; asymmetric
# families = 2 rows (both strands, forward rep then reverse-complement rep). (label, centerPos of 6mA)
ROWS = [
    ("GATC", 2), ("AATCC", 2), ("GGAGC", 3),
    ("CCATC", 3), ("GATGG", 2),
    ("CAGGAG", 5), ("CGMAGG", 4), ("GATGNAG", 6),
    ("CAGNNNNNGGA", 2), ("TCCNNNNNCTG", 2),
    ("CATNNNNNGTG", 2), ("CACNNNNNRTG", 2),
    ("CAGNNNNNCTG", 2),
    ("GACNNNNNRGAC", 2), ("GTCYNNNNNGTC", 3),
    ("CCANNNNNNCAT", 3), ("ATGNNNNNNTGG", 1),
    ("CAAYNNNNNNGTC", 3), ("GACNNNNNNRTTG", 2),
    ("GCACNNNNNNGTT", 3), ("AACNNNNNNGTGC", 2),
    ("TAANNNNNNCTTG", 3), ("CAAGNNNNNNTTA", 3),
    ("GAAGNNNNNNNTCC", 3), ("GGANNNNNNNCTTC", 3),
]
row_reps = [r for r, _ in ROWS]

SETS = {'A':set('A'),'C':set('C'),'G':set('G'),'T':set('T'),'R':set('AG'),'Y':set('CT'),
        'S':set('GC'),'W':set('AT'),'K':set('GT'),'M':set('AC'),'B':set('CGT'),'D':set('AGT'),
        'H':set('ACT'),'V':set('ACG'),'N':set('ACGT')}
def compat(x, y): return len(SETS[x] & SETS[y]) > 0
def rc(m): return str(Seq(m).reverse_complement())

def cover_score(m, ref):
    """Fraction of ref's INFORMATIVE positions that align (some shift) to a NON-N, IUPAC-compatible
    base of m. Requires all of ref's informative positions to fall within m. Aligning ref's core to
    m's N-region scores 0 (prevents a short rep hiding in a long motif's Ns); extra informative
    flanks on m (e.g. R/Y in RGATCY) are ignored, so length variants still match."""
    ref_info = [i for i in range(len(ref)) if ref[i] != 'N']
    if not ref_info: return 0.0
    best = 0.0
    for shift in range(-(len(m) - 1), len(ref)):          # ref[i] pairs with m[i-shift]
        if any(not (0 <= i - shift < len(m)) for i in ref_info): continue
        ok = sum(1 for i in ref_info
                 if m[i - shift] != 'N' and compat(ref[i], m[i - shift]))
        best = max(best, ok / len(ref_info))
    return best

def assign(motif):
    """Best-matching strand-specific row rep (direct, no revcomp -- each row is one strand)."""
    best = (None, 0.0, 0)
    for rep in row_reps:
        ninfo = sum(1 for b in rep if b != 'N')
        s = cover_score(motif, rep)
        if (s, ninfo) > (best[1], best[2]):
            best = (rep, s, ninfo)
    return best[0] if best[1] >= 0.9 else None

prof = pd.read_csv(f"{J}/motif_profile.csv", index_col=0)
prof.index = [i.rsplit("_", 1)[0] for i in prof.index]
prof = prof.groupby(level=0).max()

rows, unmapped = {}, []
for motif, row in prof.iterrows():
    rep = assign(motif)
    if rep is None: unmapped.append(motif); continue
    rows[rep] = row if rep not in rows else pd.concat([rows[rep], row], axis=1).max(axis=1)

order = [r for r in row_reps if r in rows]     # keep defined order, drop empty rows
mat = pd.DataFrame({r: rows[r] for r in order}).T
mat.index = order

ann = pd.read_csv(f"{D}/jf8_heatmap_colann.csv")
mat = mat[[c for c in ann["contig"] if c in mat.columns]]
mat.index.name = "motif"
mat.to_csv(f"{D}/jf8_wga_family_matrix.csv")

cpmap = dict(ROWS)
pd.DataFrame([(r, cpmap[r], 'A') for r in order], columns=["motif","centerPos","modified_base"]
             ).to_csv(f"{D}/jf8_wga_boldpos.csv", index=False)

print(f"rows (strand-split): {mat.shape[0]}  contigs: {mat.shape[1]}")
for r in order: print(f"  {r:<16} cp{cpmap[r]}  max={mat.loc[r].max():.3f}")
print(f"unmapped ({len(unmapped)}): {unmapped}")
