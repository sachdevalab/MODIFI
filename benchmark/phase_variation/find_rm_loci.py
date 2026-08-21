#!/usr/bin/env python3
"""Stage 1 of the phase-variation screen: locate every R-M operon in each
longitudinal strain at each time point.

For each strain in longitudinal_strains.tsv, for each time point's member
contig, read that sample's MicrobeMod R-M gene table and the prodigal
coordinates from the companion .faa headers, group genes into operons, and
emit one row per (strain, timepoint, contig, operon).

Login-node-safe: reads one R-M table + one .faa per longitudinal sample
(~14 samples), all small.

Output: /home/shuaiw/borg/revision/phase_variation/rm_loci.tsv
"""

import csv
import re
from pathlib import Path
from collections import defaultdict

RUN2 = Path("/groups/banfield/projects/multienv/methylation/data/borg/paper/run2")
OUTDIR = Path("/home/shuaiw/borg/revision/phase_variation")
STRAINS = OUTDIR / "longitudinal_strains.tsv"

FAA_HDR = re.compile(r"^>(?P<gid>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>-?1)")


def rm_dir(sample):
    for m in ("methylation4", "methylation3"):
        d = RUN2 / sample / f"{sample}_{m}" / "RM_systems"
        if (d / "all_ctgs_RM.rm.genes.tsv").exists():
            return d
    return None


def load_gene_coords(faa_path):
    """gene_id -> (start, end, strand) from prodigal .faa headers."""
    coords = {}
    with open(faa_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            m = FAA_HDR.match(line)
            if m:
                coords[m.group("gid")] = (
                    int(m.group("start")), int(m.group("end")), int(m.group("strand")))
    return coords


def load_rm_table(tsv_path):
    """list of gene dicts."""
    rows = []
    with open(tsv_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(r)
    return rows


# cache per sample so repeated timepoints/contigs don't re-read
_cache = {}


def sample_data(sample):
    if sample in _cache:
        return _cache[sample]
    d = rm_dir(sample)
    if d is None:
        _cache[sample] = (None, None)
        return None, None
    genes = load_rm_table(d / "all_ctgs_RM.rm.genes.tsv")
    coords = load_gene_coords(d / "all_ctgs_RM.rm.genes.faa")
    _cache[sample] = (genes, coords)
    return genes, coords


def main():
    strains = list(csv.DictReader(open(STRAINS), delimiter="\t"))
    out_path = OUTDIR / "rm_loci.tsv"
    n_rows = 0
    missing = []
    with open(out_path, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["cluster", "cohort", "subject", "species", "timepoint",
                    "sample", "contig", "operon", "system_types",
                    "n_MT", "n_RE", "n_SP", "operon_start", "operon_end",
                    "motifs", "gene_ids"])
        for st in strains:
            cluster = st["cluster"]
            for member in st["members"].split(";"):
                tp, contig = member.split(":", 1)
                # contig looks like infant_14_31_C ; sample is the prefix minus _<num>_<C|L>
                mm = re.match(r"^(?P<sample>.+)_\d+_[CL]$", contig)
                if not mm:
                    continue
                sample = mm.group("sample")
                genes, coords = sample_data(sample)
                if genes is None:
                    missing.append(sample)
                    continue
                # collect this contig's R-M genes grouped by operon
                ops = defaultdict(list)
                for g in genes:
                    if g["Gene"].startswith(contig + "_") or g["Gene"].rsplit("_", 1)[0] == contig:
                        # ensure the gene truly belongs to this contig (gene id = <contig>_<idx>)
                        if re.match(rf"^{re.escape(contig)}_\d+$", g["Gene"]):
                            ops[g["Operon"]].append(g)
                for operon, glist in ops.items():
                    types = sorted({g["System Type"] for g in glist})
                    n_mt = sum(1 for g in glist if g["Gene type"] == "MT")
                    n_re = sum(1 for g in glist if g["Gene type"] == "RE")
                    n_sp = sum(1 for g in glist if g["Gene type"] == "SP")
                    starts, ends = [], []
                    for g in glist:
                        c = coords.get(g["Gene"])
                        if c:
                            starts.append(c[0]); ends.append(c[1])
                    op_start = min(starts) if starts else ""
                    op_end = max(ends) if ends else ""
                    motifs = sorted({g["Homolog motif"] for g in glist
                                     if g.get("Homolog motif", "").strip()})
                    gene_ids = ";".join(g["Gene"] for g in glist)
                    w.writerow([cluster, st["cohort"], st["subject"], st["species"],
                                tp, sample, contig, operon, ",".join(types),
                                n_mt, n_re, n_sp, op_start, op_end,
                                ",".join(motifs), gene_ids])
                    n_rows += 1

    print(f"Wrote {out_path}: {n_rows} operon rows")
    if missing:
        print(f"  samples with no R-M table ({len(set(missing))}): {sorted(set(missing))}")


if __name__ == "__main__":
    main()
