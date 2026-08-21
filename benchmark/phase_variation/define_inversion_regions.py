#!/usr/bin/env python3
"""Stage 2a: for each multi-SP Type I candidate, define the invertible segment.

For every candidate strain (a dRep cluster carrying a Type I R-M operon with
>=2 specificity (SP) subunits), at each time point:
  - pull the operon's SP gene coordinates from the MicrobeMod .faa headers,
  - extract the operon window (+/- flank) from the contig fasta,
  - run EMBOSS einverted to find inverted repeats bracketing the TRD (the
    structural signature of an invertible Type I specificity locus),
  - report a candidate inversion region (SP-span, refined by the nearest
    flanking inverted-repeat pair when found).

Validation: for E. faecalis (161_4 / infant_14_31_C) the region should land on
the published inversion locus ~2323857-2327116.

Login-node-safe: einverted + seqkit subseq on ~30 small (<12 kb) windows.

Output: /home/shuaiw/borg/revision/phase_variation/candidate_inversion_regions.tsv
"""

import csv
import re
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

RUN2 = Path("/groups/banfield/projects/multienv/methylation/data/borg/paper/run2")
OUTDIR = Path("/home/shuaiw/borg/revision/phase_variation")
RM_LOCI = OUTDIR / "rm_loci.tsv"
FLANK = 1500  # bp added around the SP span before searching for inverted repeats

FAA_HDR = re.compile(r"^>(?P<gid>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>-?1)")


def rm_dir(sample):
    for m in ("methylation4", "methylation3"):
        d = RUN2 / sample / f"{sample}_{m}" / "RM_systems"
        if (d / "all_ctgs_RM.rm.genes.tsv").exists():
            return d, m
    return None, None


def contig_fasta(sample, m, contig):
    return RUN2 / sample / f"{sample}_{m}" / "contigs" / f"{contig}.fa"


def load_faa_coords(faa):
    coords = {}
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                mm = FAA_HDR.match(line)
                if mm:
                    coords[mm.group("gid")] = (int(mm.group("start")),
                                               int(mm.group("end")),
                                               int(mm.group("strand")))
    return coords


def gene_types(tsv):
    """gene_id -> (Gene type, System Type, Homolog motif)."""
    out = {}
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[r["Gene"]] = (r["Gene type"], r["System Type"],
                              r.get("Homolog motif", "").strip())
    return out


def run_einverted(region_fa):
    """Return list of (loop_start, loop_end) inverted-repeat spans (region coords)."""
    with tempfile.TemporaryDirectory() as td:
        outfile = Path(td) / "e.out"
        outseq = Path(td) / "e.fa"
        cmd = ["einverted", "-sequence", str(region_fa), "-gap", "12",
               "-threshold", "30", "-match", "3", "-mismatch", "-4",
               "-outfile", str(outfile), "-outseq", str(outseq)]
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        except Exception:
            return []
        spans = []
        if outfile.exists():
            nums = []
            for line in open(outfile):
                m = re.findall(r"(\d+)\s+[acgtnACGTN]+\s+(\d+)", line)
                for a, b in m:
                    nums.append((int(a), int(b)))
            # einverted prints the two arms; pair them into (min,max)
            for i in range(0, len(nums) - 1, 2):
                s = min(nums[i][0], nums[i][1], nums[i + 1][0], nums[i + 1][1])
                e = max(nums[i][0], nums[i][1], nums[i + 1][0], nums[i + 1][1])
                spans.append((s, e))
        return spans


def main():
    rows = list(csv.DictReader(open(RM_LOCI), delimiter="\t"))
    # candidates = (cluster) with any Type I operon n_SP>=2
    cand_clusters = sorted({r["cluster"] for r in rows
                            if "Type_I" in r["system_types"] and int(r["n_SP"]) >= 2})

    out = OUTDIR / "candidate_inversion_regions.tsv"
    with open(out, "w", newline="") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow(["cluster", "cohort", "subject", "species", "timepoint",
                    "sample", "contig", "operon", "sp_genes", "sp_span_start",
                    "sp_span_end", "mtase_motif", "n_inverted_repeats",
                    "region_start", "region_end", "contig_fasta"])
        for r in rows:
            if r["cluster"] not in cand_clusters:
                continue
            if "Type_I" not in r["system_types"] or int(r["n_SP"]) < 2:
                continue
            sample, contig, operon = r["sample"], r["contig"], r["operon"]
            d, m = rm_dir(sample)
            if d is None:
                continue
            coords = load_faa_coords(d / "all_ctgs_RM.rm.genes.faa")
            gtypes = gene_types(d / "all_ctgs_RM.rm.genes.tsv")
            gids = r["gene_ids"].split(";")
            sp_coords, mt_motifs = [], []
            for g in gids:
                gt = gtypes.get(g, ("", "", ""))
                if gt[0] == "SP" and g in coords:
                    sp_coords.append(coords[g])
                if gt[0] == "MT" and gt[2]:
                    mt_motifs.append(gt[2])
            if len(sp_coords) < 2:
                continue
            sp_start = min(c[0] for c in sp_coords)
            sp_end = max(c[1] for c in sp_coords)
            reg_s = max(1, sp_start - FLANK)
            reg_e = sp_end + FLANK

            fa = contig_fasta(sample, m, contig)
            n_ir = 0
            region_start, region_end = sp_start, sp_end
            if fa.exists():
                with tempfile.TemporaryDirectory() as td:
                    sub = Path(td) / "region.fa"
                    with open(sub, "w") as so:
                        subprocess.run(
                            ["seqkit", "subseq", "-r", f"{reg_s}:{reg_e}", str(fa)],
                            check=True, stdout=so, stderr=subprocess.DEVNULL)
                    spans = run_einverted(sub)
                    n_ir = len(spans)
                    if spans:
                        # refine region to the widest IR-bracketed span (in contig coords)
                        offs = reg_s - 1
                        s = min(sp[0] for sp in spans) + offs
                        e = max(sp[1] for sp in spans) + offs
                        region_start = min(s, sp_start)
                        region_end = max(e, sp_end)
            w.writerow([r["cluster"], r["cohort"], r["subject"], r["species"],
                        r["timepoint"], sample, contig, operon,
                        ";".join(g for g in gids if gtypes.get(g, ("",))[0] == "SP"),
                        sp_start, sp_end,
                        ",".join(sorted(set(mt_motifs))), n_ir,
                        region_start, region_end, str(fa)])
    print(f"Wrote {out}")
    print(f"Candidate clusters: {cand_clusters}")


if __name__ == "__main__":
    main()
