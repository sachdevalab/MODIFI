#!/usr/bin/env python3
"""Build LoVis4u GenBank input for the E. faecalis type I R-M inversion synteny plot.

Three loci (all subject INF1240040 / INF1240119, longitudinal E. faecalis):
  I_1 DOL 14 (G1) = infant_2_3_C   (week 2)
  I_1 DOL 35 (G2) = infant_14_31_C (week 5, same infant, later)
  I_2             = infant_26_3_C  (different infant)

Each is already extracted as an RM-locus fragment under
/home/shuaiw/borg/paper/E_faecalis/detail/<sample>_frag.{fa,faa} with the
MicrobeMod R-M call in <sample>_frag.rm.genes.tsv. Gene order in the fragment is
  _2 REase  _3 MTase  _4 S1(SP)  _5 XerC(recombinase)  _6 S2(SP)
XerC is the CDS between the two specificity (SP) subunits that MicrobeMod does not
call; we confirm it against a known XerC protein.

For each locus we slice the RM core + one flanking gene each side, build a Biopython
SeqRecord with labelled CDS features (REase/MTase/S1/XerC/S2/flank + translation), and
write one GenBank file per locus for `lovis4u -gb`. Also writes a source-data CSV.

Login-node-safe: reads 3 small fragments; one tiny blastp for XerC confirmation.
"""

import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

DETAIL = Path("/home/shuaiw/borg/paper/E_faecalis/detail")
XERC_REF = Path("/home/shuaiw/borg/paper/E_faecalis/prokka/CFMMJMDL_02264.faa")  # Tyrosine recombinase XerC
OUTDIR = Path("/home/shuaiw/borg/revision/synteny")
GBDIR = OUTDIR / "lovis_input"
SRC_CSV = OUTDIR / "efaecalis_rm_inversion_genes.csv"

# figure label -> fragment sample prefix, ordered top-to-bottom as in the original figure
LOCI = [
    ("I_1_DOL14_G1", "infant_2_frag"),
    ("I_1_DOL35_G2", "infant_14_frag"),
    ("I_2", "infant_26_frag"),
]

FLANK_GENES = 1  # number of flanking CDS to keep on each side of the RM core

FAA_HDR = re.compile(r"^>(?P<gid>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>-?1)")

PRODUCTS = {
    "REase": "Type I restriction enzyme R subunit (REase)",
    "MTase": "Type I restriction enzyme M subunit (MTase)",
    "S1": "Type I restriction enzyme S subunit (specificity, S1)",
    "S2": "Type I restriction enzyme S subunit (specificity, S2)",
    "XerC": "Tyrosine recombinase XerC",
}


def load_faa(sample):
    """gene_id -> (start, end, strand, protein_seq) from prodigal .faa headers/records."""
    coords, seqs = {}, {}
    faa = DETAIL / f"{sample}.faa"
    for rec in SeqIO.parse(faa, "fasta"):
        seqs[rec.id] = str(rec.seq).rstrip("*")
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                m = FAA_HDR.match(line)
                if m:
                    coords[m.group("gid")] = (
                        int(m.group("start")), int(m.group("end")), int(m.group("strand")))
    return {g: (coords[g][0], coords[g][1], coords[g][2], seqs[g]) for g in coords}


def load_rm(sample):
    """gene_id -> (gene_type RE/MT/SP, motif) from MicrobeMod rm.genes.tsv."""
    out = {}
    with open(DETAIL / f"{sample}.rm.genes.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r["Gene"]] = (r["Gene type"], r.get("Homolog motif", "").strip())
    return out


def load_xerc_ref():
    for rec in SeqIO.parse(XERC_REF, "fasta"):
        return str(rec.seq).rstrip("*")
    return None


def blastp_ident(query_seq, subject_seq):
    """Return best % identity of query vs subject via blastp; 0.0 if no hit / no blastp."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "q.faa").write_text(f">q\n{query_seq}\n")
        (td / "s.faa").write_text(f">s\n{subject_seq}\n")
        try:
            r = subprocess.run(
                ["blastp", "-query", str(td / "q.faa"), "-subject", str(td / "s.faa"),
                 "-outfmt", "6 pident", "-max_hsps", "1"],
                capture_output=True, text=True)
        except FileNotFoundError:
            return None  # blastp not available; caller falls back to positional logic
        vals = [float(x) for x in r.stdout.split()]
        return max(vals) if vals else 0.0


def classify_locus(sample, xerc_ref):
    """Return ordered list of gene dicts for the RM core + flanks of one fragment."""
    genes = load_faa(sample)         # gid -> (start, end, strand, prot)
    rm = load_rm(sample)             # gid -> (type, motif)
    # order genes by start position
    ordered = sorted(genes.items(), key=lambda kv: kv[1][0])

    # assign RM roles
    re_ids = [g for g, (t, _) in rm.items() if t == "RE"]
    mt_ids = [g for g, (t, _) in rm.items() if t == "MT"]
    sp_ids = sorted((g for g, (t, _) in rm.items() if t == "SP"),
                    key=lambda g: genes[g][0])  # by position -> S1 (first), S2 (second)
    role = {}
    for g in re_ids:
        role[g] = "REase"
    for g in mt_ids:
        role[g] = "MTase"
    if len(sp_ids) >= 1:
        role[sp_ids[0]] = "S1"
    if len(sp_ids) >= 2:
        role[sp_ids[-1]] = "S2"

    # XerC = CDS located between the two SP subunits, not an RM gene; confirm by identity
    xerc_id = None
    if len(sp_ids) >= 2:
        s1_start = genes[sp_ids[0]][0]
        s2_end = genes[sp_ids[-1]][1]
        lo, hi = min(s1_start, genes[sp_ids[-1]][0]), max(genes[sp_ids[0]][1], s2_end)
        between = [g for g, (s, e, st, _) in genes.items()
                   if g not in role and s >= lo - 5 and e <= hi + 5]
        # prefer the one best matching the XerC reference
        best, best_id = None, -1.0
        for g in between:
            pid = blastp_ident(genes[g][3], xerc_ref) if xerc_ref else None
            pid = pid if pid is not None else 0.0
            if pid > best_id:
                best, best_id = g, pid
        if between:
            xerc_id = best if best is not None else between[0]
            role[xerc_id] = "XerC"
            print(f"  {sample}: XerC = {xerc_id} (best blastp identity to ref = {best_id:.1f}%)")

    # core = RM genes + XerC; window = core +/- FLANK_GENES CDS
    core_ids = [g for g in role]
    core_positions = [i for i, (g, _) in enumerate(ordered) if g in core_ids]
    lo_i = max(0, min(core_positions) - FLANK_GENES)
    hi_i = min(len(ordered) - 1, max(core_positions) + FLANK_GENES)
    window = ordered[lo_i:hi_i + 1]

    out = []
    for idx, (g, (s, e, st, prot)) in enumerate(window, start=1):
        r = role.get(g, "flank")
        out.append(dict(gid=g, start=s, end=e, strand=st, prot=prot,
                        role=r, motif=rm.get(g, ("", ""))[1], idx=idx))
    return out


def build_record(label, genes):
    """Slice the fragment sequence to the gene window and build a labelled SeqRecord."""
    sample = dict(LOCI)[label]
    frag = next(SeqIO.parse(DETAIL / f"{sample}.fa", "fasta"))
    win_start = min(g["start"] for g in genes) - 100
    win_end = max(g["end"] for g in genes) + 100
    win_start = max(1, win_start)
    win_end = min(len(frag.seq), win_end)
    sub = frag.seq[win_start - 1:win_end]

    rec = SeqRecord(Seq(str(sub)), id=label, name=label,
                    description=f"E. faecalis type I R-M inversion locus ({label})",
                    annotations={"molecule_type": "DNA", "organism": "Enterococcus faecalis"})
    for g in genes:
        s = g["start"] - win_start        # 0-based within window
        e = g["end"] - win_start + 1
        strand = 1 if g["strand"] == 1 else -1
        name = g["role"]
        product = PRODUCTS.get(name, "hypothetical protein")
        if name == "flank":
            name = ""  # leave flanking genes unlabelled so the RM core stands out
        # clean, colon-free feature id (lovis4u -faf cannot parse ':' in ids)
        clean_id = f"{label}_g{g['idx']}"
        g["clean_id"] = clean_id
        feat = SeqFeature(
            FeatureLocation(max(0, s), min(len(sub), e), strand=strand),
            type="CDS",
            qualifiers={
                "gene": [name],
                "label": [name],
                "product": [product],
                "locus_tag": [clean_id],
                "protein_id": [clean_id],   # lovis4u -faf matches features by protein_id
                "translation": [g["prot"]],
            },
        )
        rec.features.append(feat)
    return rec, win_start, win_end


def main():
    GBDIR.mkdir(parents=True, exist_ok=True)
    xerc_ref = load_xerc_ref()
    if xerc_ref is None:
        print("WARNING: XerC reference protein not found; using positional logic only")

    src_rows = []
    for label, sample in LOCI:
        print(f"[{label}] {sample}")
        genes = classify_locus(sample, xerc_ref)
        rec, win_start, win_end = build_record(label, genes)
        gb = GBDIR / f"{label}.gb"
        SeqIO.write(rec, gb, "genbank")
        roles = ", ".join(f"{g['role']}({'+' if g['strand']==1 else '-'})" for g in genes)
        print(f"    window {win_start}-{win_end} ({win_end-win_start+1} bp); genes: {roles}")
        print(f"    wrote {gb}")
        for g in genes:
            src_rows.append(dict(
                locus=label, sample=sample, gene=g["clean_id"], orig_gene=g["gid"],
                role=g["role"],
                win_rel_start=g["start"] - win_start + 1,
                win_rel_end=g["end"] - win_start + 1,
                strand="+" if g["strand"] == 1 else "-",
                length_aa=len(g["prot"]), motif=g["motif"]))

    with open(SRC_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(src_rows[0].keys()))
        w.writeheader()
        w.writerows(src_rows)
    print(f"\nWrote source data: {SRC_CSV} ({len(src_rows)} genes)")
    print(f"GenBank input dir: {GBDIR}")


if __name__ == "__main__":
    sys.exit(main())
