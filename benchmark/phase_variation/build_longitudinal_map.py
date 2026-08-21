#!/usr/bin/env python3
"""Stage 0 of the phase-variation screen.

Build the within-subject longitudinal design and the list of strains tracked
across >=2 time points, from the prefix tables + the dRep 99% cross-sample
cluster table.

Outputs (to /home/shuaiw/borg/revision/phase_variation/):
  - sample_timepoint_map.tsv : sample -> cohort, subject, timepoint, raw_bam
  - longitudinal_strains.tsv : (dRep cluster, subject) present at >=2 timepoints,
                               with species and the member contig at each timepoint
"""

import re
import csv
from pathlib import Path
from collections import defaultdict

PREFIX_TABLE = Path("/home/shuaiw/MODIFI/assembly_pipe/prefix_table.tab")
CDB = Path("/groups/banfield/projects/multienv/methylation/data/borg/paper/"
           "specificity/dRep_99_out/data_tables/Cdb.csv")
CLADE = Path("/home/shuaiw/MODIFI/tmp/figures/strain_diff/drep_99/clade_data.csv")
OUTDIR = Path("/home/shuaiw/borg/revision/phase_variation")

# prefix patterns
#   infant: NANO_<batch>_INF<subject>_<week>G1_pacbio  -> subject INF<subject>, week <week>
#           NANO_<batch>_MAT<subject>_<week>G1_pacbio  -> maternal (kept, flagged)
#   asthma: TIPS_<subject>_<N>Mo_LR                    -> subject <subject>, month <N>
RE_INFANT = re.compile(r"NANO_\d+_(INF|MAT)(\d+)_(\d+)G1")
RE_ASTHMA = re.compile(r"TIPS_(\d+)_(\d+)Mo")

# genome contig name: <sample>_<contig>_<C|L>[.fa]
RE_GENOME = re.compile(r"^(?P<sample>.+)_(?P<contig>\d+)_(?P<topo>[CL])$")


def parse_prefix_table(path):
    """Return sample -> dict(cohort, subject, timepoint, timepoint_unit, raw_bam)."""
    out = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            prefix, sample, bam = fields[0], fields[1], fields[2]

            m = RE_INFANT.search(prefix)
            if m:
                kind, subj, week = m.group(1), m.group(2), int(m.group(3))
                out[sample] = dict(
                    cohort="infant" if kind == "INF" else "maternal",
                    subject=f"{kind}{subj}",
                    timepoint=week,
                    timepoint_unit="week",
                    raw_bam=bam,
                )
                continue

            m = RE_ASTHMA.search(prefix)
            if m:
                subj, month = m.group(1), int(m.group(2))
                out[sample] = dict(
                    cohort="asthma",
                    subject=subj,
                    timepoint=month,
                    timepoint_unit="month",
                    raw_bam=bam,
                )
                continue
            # non-longitudinal cohorts (soil/cow/ocean/...) are skipped
    return out


def parse_sample_from_genome(genome):
    g = genome[:-3] if genome.endswith(".fa") else genome
    m = RE_GENOME.match(g)
    if not m:
        return None, None
    return m.group("sample"), g  # sample, contig_id (without .fa)


def load_species(path):
    """cluster (secondary_cluster) -> (phylum, species) from clade_data.csv."""
    sp = {}
    if not path.exists():
        return sp
    with open(path) as f:
        for row in csv.DictReader(f):
            sp[row["cluster"]] = (row.get("phylum", ""), row.get("species", ""))
    return sp


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    samples = parse_prefix_table(PREFIX_TABLE)
    species = load_species(CLADE)

    # write sample_timepoint_map.tsv
    map_path = OUTDIR / "sample_timepoint_map.tsv"
    with open(map_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "cohort", "subject", "timepoint", "timepoint_unit", "raw_bam"])
        for s in sorted(samples):
            d = samples[s]
            w.writerow([s, d["cohort"], d["subject"], d["timepoint"],
                        d["timepoint_unit"], d["raw_bam"]])

    # subjects with >=2 timepoints (longitudinal), infant + asthma only
    tp_by_subject = defaultdict(set)
    for s, d in samples.items():
        if d["cohort"] in ("infant", "asthma"):
            tp_by_subject[(d["cohort"], d["subject"])].add(d["timepoint"])
    longitudinal_subjects = {k for k, tps in tp_by_subject.items() if len(tps) >= 2}

    # walk Cdb: cluster -> subject -> list of (timepoint, sample, contig)
    cluster_subj = defaultdict(lambda: defaultdict(list))
    with open(CDB) as f:
        for row in csv.DictReader(f):
            genome = row["genome"]
            cluster = row["secondary_cluster"]
            sample, contig = parse_sample_from_genome(genome)
            if sample is None or sample not in samples:
                continue
            d = samples[sample]
            if d["cohort"] not in ("infant", "asthma"):
                continue
            key = (d["cohort"], d["subject"])
            if key not in longitudinal_subjects:
                continue
            cluster_subj[cluster][key].append((d["timepoint"], sample, contig))

    # keep (cluster, subject) present at >=2 distinct timepoints
    strains_path = OUTDIR / "longitudinal_strains.tsv"
    n_strains = 0
    with open(strains_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["cluster", "cohort", "subject", "phylum", "species",
                    "n_timepoints", "timepoints", "members"])
        for cluster in sorted(cluster_subj):
            phy, sp = species.get(cluster, ("", ""))
            for (cohort, subject), recs in sorted(cluster_subj[cluster].items()):
                tps = sorted({r[0] for r in recs})
                if len(tps) < 2:
                    continue
                n_strains += 1
                members = ";".join(f"{tp}:{cont}" for tp, samp, cont in sorted(recs))
                w.writerow([cluster, cohort, subject, phy, sp,
                            len(tps), ",".join(map(str, tps)), members])

    # summary to screen
    print(f"Wrote {map_path}")
    print(f"  longitudinal samples: "
          f"{sum(1 for d in samples.values() if d['cohort'] in ('infant','asthma'))}")
    print(f"  longitudinal subjects (>=2 timepoints): {len(longitudinal_subjects)}")
    for (cohort, subj), tps in sorted(tp_by_subject.items()):
        if (cohort, subj) in longitudinal_subjects:
            print(f"    {cohort:7s} {subj:12s} timepoints={sorted(tps)}")
    print(f"Wrote {strains_path}")
    print(f"  longitudinal strains (cluster x subject, >=2 timepoints): {n_strains}")


if __name__ == "__main__":
    main()
