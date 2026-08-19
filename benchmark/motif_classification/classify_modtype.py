#!/usr/bin/env python3
"""
Modification-type classification for MODIFI motifs (6mA / 4mC / 5mC).

The type of a detected motif is decided from two pieces of evidence:

  1. The modified base, taken from MODIFI's own detection: base = motifString[centerPos-1]
     (strand-aware). A modified A -> 6mA; a modified C -> 4mC or 5mC.
  2. Strain/species-level alignment against REBASE: the detected recognition sequence is matched
     (IUPAC-aware, both strands) ONLY against REBASE records whose organism matches the contig's
     strain/species, and the modification type is read from REBASE's <4> METHYLATION SITE field
     ((6)=6mA, (5)=5mC, (4)=4mC). This resolves 4mC vs 5mC for cytosine motifs and confirms 6mA.

Cross-organism REBASE matches are NOT trusted for cytosine motifs (the same recognition sequence
is used with different chemistries in different organisms); they are recorded as advisory only.
An optional strain methylome table (a SMRT *_motif_summary.csv / *.motifs.gff carrying a real
modType) can override, and is treated as authoritative when present.

Usage
-----
  # per-contig organism map (contig<TAB>organism), e.g. built for JF8 or PB24
  python classify_modtype.py \
      --motifs   <all.motifs.csv | all.motifs.merged.csv> \
      --org-map  <contig_species.tsv> \
      --rebase   /home/shuaiw/MODIFI/benchmark/rebase/withrefm.txt \
      --out      <out.tsv>

  # single organism for the whole file (e.g. an isolate genome)
  python classify_modtype.py --motifs all.motifs.csv --organism "Escherichia coli" \
      --rebase .../withrefm.txt --out out.tsv

The motif file must contain columns motifString and centerPos (MODIFI all.motifs.csv or the
merged file with a leading `contig` column). If a `contig` column is absent, --organism is
required.
"""

import argparse
import csv
import re
import sys
import os

# ---------------------------------------------------------------------------
# IUPAC handling (kept in-sync with scripts/motif_profile.py IUPAC_CODES)
# ---------------------------------------------------------------------------
IUPAC_SETS = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}
COMPLEMENT = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'S': 'S', 'W': 'W',
    'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N',
}
TYPE_CODE = {'6': '6mA', '5': '5mC', '4': '4mC'}


def reverse_complement(motif):
    try:
        return ''.join(COMPLEMENT[b] for b in reversed(motif))
    except KeyError:
        return None


def iupac_overlap(a, b):
    """True if IUPAC symbols a and b can represent the same nucleotide."""
    return bool(set(IUPAC_SETS.get(a, a)) & set(IUPAC_SETS.get(b, b)))


def recseq_equal(a, b):
    """IUPAC-aware equality of two equal-length recognition sequences."""
    return len(a) == len(b) and all(iupac_overlap(x, y) for x, y in zip(a, b))


# ---------------------------------------------------------------------------
# REBASE withrefm.txt parsing
# ---------------------------------------------------------------------------
def parse_rebase(withrefm_path):
    """Return a list of records: dict(enzyme, recseq, org(lower), types(set), meth_positions).

    Only records with a clean recognition sequence (no cleavage coordinates) and a parseable
    <4> METHYLATION SITE type are kept.
    """
    records = []
    cur = {}
    with open(withrefm_path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("<1>"):
                cur = {"enzyme": line[3:]}
            elif line.startswith("<3>"):
                cur["recseq"] = line[3:].replace("^", "")
            elif line.startswith("<4>"):
                cur["meth"] = line[3:]
            elif line.startswith("<5>"):
                cur["org"] = line[3:]
            elif line.startswith("<8>"):
                rec = _finalize_rebase(cur)
                if rec:
                    records.append(rec)
                cur = {}
    return records


def _finalize_rebase(cur):
    recseq = cur.get("recseq", "")
    meth = cur.get("meth", "")
    if not recseq or "(" in recseq or reverse_complement(recseq) is None:
        return None
    types = set()
    positions = []  # (position_1based_signed, type)
    for m in re.finditer(r"(-?\d+)\((\d)\)", meth):
        t = TYPE_CODE.get(m.group(2))
        if not t:
            continue
        types.add(t)
        positions.append((int(m.group(1)), t))
    if not types:
        return None
    return {
        "enzyme": cur.get("enzyme", ""),
        "recseq": recseq,
        "org": cur.get("org", "").lower(),
        "types": types,
        "positions": positions,
    }


def genus_species(name):
    """Lowercased 'Genus species' key from a scientific name / strain string."""
    parts = name.replace("[", "").replace("]", "").split()
    return " ".join(parts[:2]).lower() if len(parts) >= 2 else name.lower()


# extra organism-name aliases so strain-scoped matching survives taxonomic renames
ALIASES = {
    "bacteroides vulgatus": ["phocaeicola vulgatus"],
    "phocaeicola vulgatus": ["bacteroides vulgatus"],
    "clostridium bolteae": ["enterocloster bolteae"],
    "enterocloster bolteae": ["clostridium bolteae"],
    "ruminococcus gnavus": ["mediterraneibacter gnavus"],
    "mediterraneibacter gnavus": ["ruminococcus gnavus"],
}


def org_keys(organism):
    """Return the set of organism substrings to scope REBASE matching to this strain/species."""
    gs = genus_species(organism)
    keys = {gs}
    keys.update(ALIASES.get(gs, []))
    return keys


# ---------------------------------------------------------------------------
# Classification core
# ---------------------------------------------------------------------------
def modified_base(motif, center_pos):
    """Return (base, strand) of the modified position.

    center_pos is 1-based on the given motif strand. If it points to A/C -> forward strand.
    If it points to T/G (a real modification is always on A or C), the modified base is on the
    complementary strand; report the complementary base and strand '-'.
    """
    if not (1 <= center_pos <= len(motif)):
        return "?", "?"
    b = motif[center_pos - 1]
    if b in ("A", "C"):
        return b, "+"
    if b in ("T", "G"):
        return COMPLEMENT[b], "-"
    # ambiguous IUPAC at the modified position -> take its allowed bases
    allowed = IUPAC_SETS.get(b, "")
    if "A" in allowed or "C" in allowed:
        return b, "+"
    return b, "?"


def rebase_types_for(motif, records, org_key_set):
    """Return (types_set, evidence) matching `motif` against REBASE records scoped to org_key_set.

    Matches the motif and its reverse complement against each record's recognition sequence,
    IUPAC-aware. If org_key_set is None, matches all organisms (advisory).
    """
    types = set()
    evidence = []
    rc = reverse_complement(motif)
    for rec in records:
        if org_key_set is not None and not any(k in rec["org"] for k in org_key_set):
            continue
        for cand in (motif, rc):
            if cand and recseq_equal(cand, rec["recseq"]):
                types |= rec["types"]
                evidence.append((rec["enzyme"], rec["recseq"], ",".join(sorted(rec["types"])), rec["org"]))
                break
    return types, evidence


def classify(motif, center_pos, organism, records, methylome=None):
    """Classify one motif. Returns a dict with the call, confidence, and evidence.

    Priority of evidence:
      1. Strain methylome modType (authoritative) if provided and it has this motif.
      2. Strain/species-scoped REBASE type combined with the modified base.
      3. Modified base alone (A -> 6mA) when REBASE has no strain-scoped record.
    Cross-organism REBASE is advisory only and never sets a cytosine type by itself.
    """
    base, strand = modified_base(motif, center_pos)
    org_key_set = org_keys(organism)

    # 1. strain methylome override
    if methylome is not None:
        mt = methylome.get((organism_key(organism), motif))
        if mt:
            return _result(motif, center_pos, base, strand, mt, "high", "strain_methylome",
                           evidence=f"{organism} methylome modType={mt}")

    # 2. strain-scoped REBASE
    scoped_types, scoped_ev = rebase_types_for(motif, records, org_key_set)
    any_types, any_ev = rebase_types_for(motif, records, None)

    ev_str = "; ".join(f"{e[0]}:{e[1]}({e[2]})" for e in (scoped_ev or any_ev)[:3])

    if base == "A":
        # adenine modification is unambiguous: 6mA. REBASE (scoped) should confirm.
        conf = "high" if "6mA" in scoped_types else "high"
        src = "base=A+rebase" if "6mA" in scoped_types else "base=A"
        return _result(motif, center_pos, base, strand, "6mA", conf, src, ev_str)

    if base == "C":
        c_scoped = {t for t in scoped_types if t in ("4mC", "5mC")}
        if len(c_scoped) == 1:
            return _result(motif, center_pos, base, strand, next(iter(c_scoped)),
                           "high", "strain_rebase", ev_str)
        if len(c_scoped) == 2:
            return _result(motif, center_pos, base, strand, "4mC/5mC",
                           "ambiguous", "strain_rebase_conflict", ev_str)
        # no strain-scoped cytosine type: report advisory + flag for methylome resolution
        c_any = {t for t in any_types if t in ("4mC", "5mC")}
        advisory = "/".join(sorted(c_any)) if c_any else ""
        return _result(motif, center_pos, base, strand, "mod-C(unresolved)",
                       "needs_strain_methylome", "base=C",
                       ev_str, advisory=advisory)

    # ambiguous / non-AC modified base
    return _result(motif, center_pos, base, strand,
                   "/".join(sorted(scoped_types)) if scoped_types else "unassigned",
                   "low", "rebase_only" if scoped_types else "none", ev_str)


def organism_key(organism):
    return genus_species(organism)


def _result(motif, center_pos, base, strand, mod_type, confidence, source, evidence, advisory=""):
    return {
        "motifString": motif,
        "centerPos": center_pos,
        "modified_base": base,
        "strand": strand,
        "mod_type": mod_type,
        "confidence": confidence,
        "source": source,
        "advisory_type": advisory,
        "rebase_evidence": evidence,
    }


# ---------------------------------------------------------------------------
# Optional strain methylome loader (SMRT motif_summary.csv with modificationType)
# ---------------------------------------------------------------------------
def load_methylomes(methylome_dir):
    """Load any *_motif_summary.csv in methylome_dir that has a modificationType column.

    Returns {(genus_species_key, motifString): modType}. File names are expected to embed the
    organism; a companion manifest.tsv (file<TAB>organism) is used if present.
    """
    if not methylome_dir or not os.path.isdir(methylome_dir):
        return None
    manifest = {}
    man_path = os.path.join(methylome_dir, "manifest.tsv")
    if os.path.exists(man_path):
        for row in csv.reader(open(man_path), delimiter="\t"):
            if len(row) >= 2:
                manifest[os.path.basename(row[0])] = row[1]
    out = {}
    for fn in os.listdir(methylome_dir):
        if not fn.endswith(".csv"):
            continue
        path = os.path.join(methylome_dir, fn)
        org = manifest.get(fn) or _org_from_filename(fn)
        if not org:
            continue
        key = genus_species(org)
        try:
            reader = csv.DictReader(open(path))
        except Exception:
            continue
        if not reader.fieldnames or "modificationType" not in reader.fieldnames:
            continue
        for r in reader:
            mo = r.get("motifString")
            mt = _norm_type(r.get("modificationType", ""))
            if mo and mt:
                out[(key, mo)] = mt
    return out or None


def _org_from_filename(fn):
    # e.g. GSM1711676_Methanocorpusculum_labreanum_Z_native.motif_summary.csv
    m = re.search(r"GSM\d+_([A-Za-z]+_[a-z]+)", fn)
    return m.group(1).replace("_", " ") if m else None


def _norm_type(v):
    v = v.strip().lower()
    return {"m6a": "6mA", "6ma": "6mA", "m4c": "4mC", "4mc": "4mC",
            "m5c": "5mC", "5mc": "5mC"}.get(v, "")


# ---------------------------------------------------------------------------
# Organism map + motif loading
# ---------------------------------------------------------------------------
def load_org_map(path):
    m = {}
    for row in csv.reader(open(path), delimiter="\t"):
        if len(row) >= 2:
            m[row[0]] = row[1]
    return m


def iter_motifs(motif_file):
    """Yield (contig_or_None, motifString, centerPos) from a MODIFI motif CSV."""
    reader = csv.DictReader(open(motif_file))
    fields = reader.fieldnames or []
    has_contig = "contig" in fields
    for r in reader:
        mo = r.get("motifString")
        cp = r.get("centerPos")
        if not mo or cp in (None, ""):
            continue
        try:
            cp = int(cp)
        except ValueError:
            continue
        contig = r.get("contig") if has_contig else None
        yield contig, mo, cp


def main():
    ap = argparse.ArgumentParser(description="Classify MODIFI motifs into 6mA/4mC/5mC.")
    ap.add_argument("--motifs", required=True, help="MODIFI all.motifs.csv or merged csv")
    ap.add_argument("--rebase", required=True, help="REBASE withrefm.txt")
    ap.add_argument("--org-map", help="TSV contig<TAB>organism (for merged/per-contig files)")
    ap.add_argument("--organism", help="single organism name for the whole file")
    ap.add_argument("--methylome-dir", help="dir of strain SMRT motif_summary.csv (+manifest.tsv)")
    ap.add_argument("--out", required=True, help="output TSV")
    args = ap.parse_args()

    if not args.org_map and not args.organism:
        ap.error("provide either --org-map or --organism")

    records = parse_rebase(args.rebase)
    sys.stderr.write(f"[classify] loaded {len(records)} REBASE records with a methylation type\n")
    methylomes = load_methylomes(args.methylome_dir)
    if methylomes:
        sys.stderr.write(f"[classify] loaded {len(methylomes)} strain-methylome motif types\n")
    org_map = load_org_map(args.org_map) if args.org_map else None

    seen = set()
    rows = []
    for contig, motif, cp in iter_motifs(args.motifs):
        organism = args.organism or (org_map.get(contig) if org_map else None)
        if not organism:
            continue  # contig not assigned to a strain/species -> skip
        key = (contig, motif, cp, organism)
        if key in seen:
            continue
        seen.add(key)
        res = classify(motif, cp, organism, records, methylomes)
        res["contig"] = contig or ""
        res["organism"] = organism
        rows.append(res)

    cols = ["organism", "contig", "motifString", "centerPos", "modified_base", "strand",
            "mod_type", "confidence", "source", "advisory_type", "rebase_evidence"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # brief summary to stderr
    from collections import Counter
    dist = Counter(r["mod_type"] for r in rows)
    sys.stderr.write(f"[classify] wrote {len(rows)} rows -> {args.out}\n")
    sys.stderr.write(f"[classify] type distribution: {dict(dist)}\n")


if __name__ == "__main__":
    main()
