#!/usr/bin/env python3
"""Assembly-only SV scan in R-M loci across time points (no methylation used).

For each longitudinal strain, take every R-M operon at its earliest time point
and blastn the operon region (+/- flank) against each later time point's
assembled contig(s). Classify the structural relationship:
  colinear   - full-length hit, same strand         (no SV)
  inversion  - full-length hit, opposite strand      (TRD/segment inversion)
  rearranged - hit split across strands / distant blocks
  indel      - partial hit (0.3-0.8 of region aligns) (insertion/deletion)
  absent     - little/no hit (<0.3)                   (gene loss / big change)

Quick: blastn -subject (no db build), small operon queries.

Output: /home/shuaiw/borg/revision/phase_variation/rm_sv_assembly.tsv
"""

import csv
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

RUN2 = Path("/groups/banfield/projects/multienv/methylation/data/borg/paper/run2")
BASE = Path("/home/shuaiw/borg/revision/phase_variation")
RM_LOCI = BASE / "rm_loci.tsv"
FLANK = 500


def meth_dir(sample):
    for m in ("methylation4", "methylation3"):
        d = RUN2 / sample / f"{sample}_{m}"
        if d.exists():
            return d
    return None


def contig_fa(sample, contig):
    d = meth_dir(sample)
    return d / "contigs" / f"{contig}.fa" if d else None


def sample_of(contig):
    m = re.match(r"^(.+)_\d+_[CL]$", contig)
    return m.group(1) if m else None


def subseq(fa, contig, start, end, out):
    with open(out, "w") as o:
        subprocess.run(["seqkit", "subseq", "-r", f"{start}:{end}", str(fa)],
                       check=True, stdout=o, stderr=subprocess.DEVNULL)


def blastn(query_fa, subject_fa):
    cmd = ["blastn", "-query", str(query_fa), "-subject", str(subject_fa),
           "-outfmt", "6 qseqid qlen qstart qend sstart send length pident sstrand",
           "-max_hsps", "50", "-evalue", "1e-10"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    hsps = []
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 9:
            continue
        qs, qe = int(f[2]), int(f[3])
        hsps.append(dict(qlen=int(f[1]), qstart=min(qs, qe), qend=max(qs, qe),
                         length=int(f[6]), pident=float(f[7]), sstrand=f[8]))
    return hsps


def _merge(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def _len(ivs):
    return sum(e - s + 1 for s, e in ivs)


def _subtract(a, b):
    """interval list a minus interval list b."""
    res = []
    for s, e in a:
        cur = [(s, e)]
        for bs, be in b:
            nxt = []
            for cs, ce in cur:
                if be < cs or bs > ce:
                    nxt.append((cs, ce))
                else:
                    if cs < bs:
                        nxt.append((cs, bs - 1))
                    if ce > be:
                        nxt.append((be + 1, ce))
            cur = nxt
        res.extend(cur)
    return _merge([list(x) for x in res])


def classify(hsps):
    """Query-coordinate aware classification.

    Single-strand full-length hit => colinear (contig global orientation is
    arbitrary between assemblies, so minus-only is NOT a biological inversion).
    A genuine localized inversion = an internal query segment covered by the
    OPPOSITE strand while the flanks are covered by the other strand.
    plus and minus covering the SAME query region => inverted-repeat/palindrome
    (not a rearrangement between the two assemblies).
    """
    if not hsps:
        return "absent", 0.0, ""
    qlen = hsps[0]["qlen"]
    plus = _merge([[h["qstart"], h["qend"]] for h in hsps if h["sstrand"] == "plus"])
    minus = _merge([[h["qstart"], h["qend"]] for h in hsps if h["sstrand"] == "minus"])
    cov = min(1.0, _len(_merge(plus + minus)) / qlen)
    note = f"cov={cov:.2f};plus_bp={_len(plus)};minus_bp={_len(minus)}"
    if cov < 0.3:
        return "absent", cov, note
    # strand that is dominant becomes the "reference" orientation; the other
    # strand, where it covers query NOT covered by the dominant strand, is an
    # inverted internal segment.
    if _len(plus) >= _len(minus):
        maj, mino = plus, minus
    else:
        maj, mino = minus, plus
    inv_only = _subtract(mino, maj)          # minority strand, distinct region
    inv_only = [iv for iv in inv_only if iv[1] - iv[0] + 1 >= 400]
    if inv_only:
        # is it internal (major strand present on both sides in query)?
        lo = min(s for s, e in inv_only)
        hi = max(e for s, e in inv_only)
        left = any(e < lo for s, e in maj)
        right = any(s > hi for s, e in maj)
        note += f";inv_seg={_len(inv_only)}bp;internal={left and right}"
        return ("inversion" if (left and right) else "rearranged"), cov, note
    # both strands over the same region => palindrome/IR, not a real change
    if _len(mino) >= 0.3 * qlen and _len(_subtract(mino, maj)) < 400:
        return ("palindrome" if cov >= 0.8 else "indel"), cov, note
    if cov < 0.8:
        return "indel", cov, note
    return "colinear", cov, note


def main():
    rows = list(csv.DictReader(open(RM_LOCI), delimiter="\t"))
    by_strain = defaultdict(list)
    for r in rows:
        by_strain[(r["cluster"], r["subject"])].append(r)

    out = BASE / "rm_sv_assembly.tsv"
    n_sv = 0
    with open(out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["cluster", "subject", "species", "operon", "system_types",
                    "n_SP", "ref_tp", "ref_contig", "ref_start", "ref_end",
                    "query_tp", "query_contig", "sv_class", "coverage", "note"])
        for (cluster, subject), rrows in by_strain.items():
            species = rrows[0]["species"]
            tps = sorted({int(r["timepoint"]) for r in rrows})
            ref_tp = tps[0]
            ref_rows = [r for r in rrows if int(r["timepoint"]) == ref_tp
                        and r["operon_start"] and r["operon_end"]]
            # query contigs per later timepoint
            later = {}
            for tp in tps[1:]:
                contigs = sorted({r["contig"] for r in rrows if int(r["timepoint"]) == tp})
                later[tp] = contigs
            if not later:
                continue
            with tempfile.TemporaryDirectory() as td:
                td = Path(td)
                for rr in ref_rows:
                    fa = contig_fa(sample_of(rr["contig"]), rr["contig"])
                    if not fa or not fa.exists():
                        continue
                    s = max(1, int(rr["operon_start"]) - FLANK)
                    e = int(rr["operon_end"]) + FLANK
                    q = td / "op.fa"
                    subseq(fa, rr["contig"], s, e, q)
                    for tp, contigs in later.items():
                        subj = td / f"subj_{tp}.fa"
                        with open(subj, "w") as so:
                            for c in contigs:
                                cf = contig_fa(sample_of(c), c)
                                if cf and cf.exists():
                                    so.write(cf.read_text())
                        hsps = blastn(q, subj)
                        sv, cov, note = classify(hsps)
                        if sv not in ("colinear",):
                            n_sv += 1
                        w.writerow([cluster, subject, species, rr["operon"],
                                    rr["system_types"], rr["n_SP"], ref_tp,
                                    rr["contig"], s, e, tp,
                                    ";".join(contigs), sv, f"{cov:.2f}", note])
    print(f"Wrote {out}")
    print(f"  non-colinear R-M locus comparisons: {n_sv}")


if __name__ == "__main__":
    main()
