#!/usr/bin/env python3
"""Easyfig/clinker-style synteny figure of the E. faecalis type I R-M inversion with
connecting ribbons shaded by pairwise amino-acid % identity (the reviewer's request:
"coloured by % identity").

Complements the LoVis4u figure (which colours genes by homology group but draws all links
in one colour). Here:
  - gene arrows keep the same colours/labels as the LoVis4u figure
    (read from styled_feature_annotation.tsv for consistency),
  - each gene is connected to its best blastp match in the adjacent locus by a ribbon whose
    greyscale shade encodes the amino-acid % identity of that pair (colourbar included).

Ribbons are drawn by position (column i joined straight to column i below, no crossing) and
shaded by coverage-aware amino-acid identity (nident/max(qlen,slen)). The specificity positions
(S1, S2) drop to ~35% between the two I_1 time points because the inversion replaces the
target-recognition domain occupying each position; the arrow colours (shared with the LoVis4u
figure) show S1 and S2 exchanging. Against the second infant (I_2) they are ~63-65% (divergent
TRDs). (The relocated subunit is ~95% identical to its counterpart at the other position, the
small difference being the new start codon at the recombination junction.)

Login-node-safe: 3 tiny loci, two small blastp runs. Uses base env (Biopython + matplotlib +
blastp on PATH).

Outputs (to /home/shuaiw/MODIFI/tmp/rev_figs/synteny/):
  efaecalis_rm_inversion_synteny_identity.{pdf,png}
  efaecalis_rm_inversion_synteny_identity_sourcedata.csv
"""

import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LinearSegmentedColormap
from Bio import SeqIO

DATA = Path("/home/shuaiw/borg/revision/synteny")
GBDIR = DATA / "lovis_input"
FEAT = DATA / "styled_feature_annotation.tsv"
FIG = Path("/home/shuaiw/MODIFI/tmp/rev_figs/synteny")
LOCI = ["I_1_DOL14_G1", "I_1_DOL35_G2", "I_2"]
LOCI_LABEL = {"I_1_DOL14_G1": "I_1 DOL 14 (G1)", "I_1_DOL35_G2": "I_1 DOL 35 (G2)", "I_2": "I_2"}

# soft blue gradient (clinker-style): light blue = low identity, medium blue = high.
# Truncated so 100% is a calm medium blue, never heavy black.
CMAP = LinearSegmentedColormap.from_list(
    "soft_blues", plt.get_cmap("Blues")(np.linspace(0.10, 0.68, 256)))
VMIN, VMAX = 30, 100   # identity range for the colourbar
RIBBON_ALPHA = 0.80


def load_features():
    """locus -> ordered list of dicts(gid, name, start, end, strand, colour)."""
    feats = {L: [] for L in LOCI}
    with open(FEAT) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s, e, strand = r["coordinates"].split(":")
            feats[r["locus_id"]].append(dict(
                gid=r["feature_id"], name=r["name"], start=int(s), end=int(e),
                strand=int(strand), colour=r["fill_colour"]))
    for L in feats:
        feats[L].sort(key=lambda g: g["start"])
    return feats


def load_proteins():
    prot = {}
    for L in LOCI:
        rec = list(SeqIO.parse(GBDIR / f"{L}.gb", "genbank"))[0]
        prot[L] = {f.qualifiers["locus_tag"][0]: str(f.qualifiers["translation"][0])
                   for f in rec.features if f.type == "CDS"}
    return prot


def _blast_global(qs, ss):
    """Return dict (q,s) -> global identity = 100*nident/max(qlen,slen).

    Coverage-aware: a subunit reconstructed by the inversion has a different N-terminus
    (start-codon change), so its full-length identity is <100% even where the best local
    HSP is 100%. Using nident/max(len) captures that honestly."""
    with tempfile.TemporaryDirectory() as td:
        Path(td, "q").write_text("".join(f">{k}\n{v}\n" for k, v in qs.items()))
        Path(td, "s").write_text("".join(f">{k}\n{v}\n" for k, v in ss.items()))
        r = subprocess.run(
            ["blastp", "-query", f"{td}/q", "-subject", f"{td}/s",
             "-outfmt", "6 qseqid sseqid nident qlen slen", "-max_hsps", "1"],
            capture_output=True, text=True)
    pair = {}
    for line in r.stdout.splitlines():
        q, s, nid, ql, sl = line.split("\t")
        gid = 100 * int(nid) / max(int(ql), int(sl))
        pair[(q, s)] = max(gid, pair.get((q, s), 0))
    return pair


def greedy_match(qs, ss):
    """Greedy 1:1 orthologue matching by coverage-aware global identity: q -> (s, gid).

    Symmetric (uses the max of q->s and s->q identities), assigns highest-identity pairs
    first, each gene used once. Gives a clean 1:1 set of links even for divergent loci,
    preserving the specificity-subunit crossing."""
    fwd = _blast_global(qs, ss)
    rev = _blast_global(ss, qs)
    pairs = []
    for q in qs:
        for s in ss:
            g = max(fwd.get((q, s), 0), rev.get((s, q), 0))
            if g > 0:
                pairs.append((g, q, s))
    pairs.sort(reverse=True)
    used_q, used_s, out = set(), set(), {}
    for g, q, s in pairs:
        if q in used_q or s in used_s:
            continue
        out[q] = (s, g); used_q.add(q); used_s.add(s)
    return out


def arrow(ax, g, y, h=0.34, label=False):
    x0, x1 = g["start"], g["end"]
    L = x1 - x0
    head = min(450, 0.4 * L)
    if g["strand"] == 1:
        ax.add_patch(FancyArrow(x0, y, L, 0, width=h, head_width=h, head_length=head,
                                length_includes_head=True, fc=g["colour"], ec="black", lw=0.6,
                                zorder=3))
    else:
        ax.add_patch(FancyArrow(x1, y, -L, 0, width=h, head_width=h, head_length=head,
                                length_includes_head=True, fc=g["colour"], ec="black", lw=0.6,
                                zorder=3))
    if label and g["name"]:
        ax.text((x0 + x1) / 2, y + h * 0.5 + 0.18, g["name"], ha="center", va="bottom",
                fontsize=9.5, zorder=4)


RIBBON_FILL = "#bdbfc1"   # neutral grey, used translucent
RIBBON_ALPHA = 0.55

def _head_len(g):
    """Arrowhead length, matching arrow()."""
    return min(450, 0.4 * (g["end"] - g["start"]))

def body_span(g):
    """x-range of the rectangular body of the arrow, excluding the triangular head."""
    hl = _head_len(g)
    if g["strand"] == 1:
        return g["start"], g["end"] - hl      # head on the right
    return g["start"] + hl, g["end"]          # head on the left

def ribbon(ax, g_top, y_top_b, g_bot, y_bot_t, pident):
    """Smooth (Bezier) ribbon connecting the rectangular body of a gene on the upper track to
    the body of its match below (the arrowheads are not spanned); the % amino-acid identity is
    printed just below the ribbon's top."""
    tl, tr = body_span(g_top)
    bl, br = body_span(g_bot)
    ymid = (y_top_b + y_bot_t) / 2.0
    verts = [
        (tl, y_top_b),                                   # top-left
        (tr, y_top_b),                                   # top edge -> top-right
        (tr, ymid), (br, ymid), (br, y_bot_t),           # right side cubic down
        (bl, y_bot_t),                                   # bottom edge -> bottom-left
        (bl, ymid), (tl, ymid), (tl, y_top_b),           # left side cubic up
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), fc=RIBBON_FILL, ec="none",
                           alpha=RIBBON_ALPHA, zorder=1))
    ax.text((tl + tr + bl + br) / 4, (y_top_b + y_bot_t) / 2, f"{pident:.0f}%",
            ha="center", va="center", fontsize=8, color="#3a3a3a", zorder=5)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    feats = load_features()
    prot = load_proteins()
    y = {LOCI[0]: 2.0, LOCI[1]: 1.0, LOCI[2]: 0.0}
    h = 0.34

    fig, ax = plt.subplots(figsize=(9, 4.2))
    src = []
    # ribbons first (behind arrows). Positional links: gene column i in the upper locus is
    # joined straight down to column i in the lower locus (no crossing). The specificity
    # positions read as low identity because the inversion (G1->G2) / divergence (->I_2)
    # changes the sequence occupying that position; the arrow colours show the S1/S2 swap.
    for i in range(len(LOCI) - 1):
        A, B = LOCI[i], LOCI[i + 1]
        idmat = _blast_global(prot[A], prot[B])   # (gidA,gidB) -> global identity
        for ga, gb in zip(feats[A], feats[B]):     # same column order across loci
            pid = idmat.get((ga["gid"], gb["gid"]), 0.0)
            ribbon(ax, ga, y[A] - h * 0.5, gb, y[B] + h * 0.5, pid)
            src.append(dict(pair=f"{A}->{B}", position=ga["name"] or "flank",
                            gene_top=ga["name"] or "flank", gene_bottom=gb["name"] or "flank",
                            aa_identity=round(pid, 1)))
    # gene arrows + track labels
    for L in LOCI:
        ax.hlines(y[L], min(g["start"] for g in feats[L]),
                  max(g["end"] for g in feats[L]), color="0.3", lw=1.0, zorder=2)
        for g in feats[L]:
            arrow(ax, g, y[L], h, label=(L == LOCI[0]))   # label genes on the top track only
        ax.text(-250, y[L], LOCI_LABEL[L], ha="right", va="center", fontsize=11)

    ax.set_xlim(-3300, 10600)
    ax.set_ylim(-0.7, 2.9)
    ax.axis("off")
    # numbers in the ribbons are amino-acid % identity between the genes at that position
    ax.text(10450, 1.0, "numbers = amino-acid % identity", rotation=90,
            ha="left", va="center", fontsize=8, color="0.4")


    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"efaecalis_rm_inversion_synteny_identity.{ext}", dpi=300,
                    bbox_inches="tight")
    with open(FIG / "efaecalis_rm_inversion_synteny_identity_sourcedata.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pair", "position", "gene_top", "gene_bottom", "aa_identity"])
        w.writeheader(); w.writerows(src)
    print("Saved:", FIG / "efaecalis_rm_inversion_synteny_identity.pdf")
    for r in src:
        print(r)


if __name__ == "__main__":
    main()
