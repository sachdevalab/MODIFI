#!/usr/bin/env python3
"""
Extract a small subset of viral structural-marker HMM models from a large Pfam-A.hmm
by accession, without needing an hmmfetch SSI index (streams the file record by record).

Writes:
  <out_hmm>              -- concatenated HMMER3 records for the requested accessions
  <out_dir>/marker_map.tsv -- accession<TAB>class<TAB>pfam_name<TAB>has_ga

Usage:
  python build_viral_hmm.py <pfam_a_hmm> <out_hmm>
"""
import os
import sys

# Curated Pfam accessions for phage structural markers -> marker class.
# (verified present in Pfam-A 37 with matching DESC lines)
VIRAL_CLASSES = {
    "PF03354": "terminase_large",   # Terminase large subunit, ATPase domain
    "PF04466": "terminase_large",   # Phage terminase large subunit
    "PF05876": "terminase_large",   # Phage terminase large subunit gpA, ATPase domain
    "PF03237": "terminase_large",   # Terminase large subunit, T4likevirus-type, N-terminal
    "PF03592": "terminase_small",   # Terminase small subunit
    "PF05119": "terminase_small",   # Phage terminase, small subunit
    "PF05065": "major_capsid",      # Phage capsid family
    "PF04233": "major_capsid",      # Phage Mu protein F like protein
    "PF03864": "major_capsid",      # Phage major capsid protein E
    "PF04860": "portal",            # Phage portal protein
    "PF05136": "portal",            # Phage portal protein, lambda family
}

# Plasmid hallmark markers: replication-initiation + partition (ParA/ParB) proteins.
PLASMID_CLASSES = {
    "PF01446": "rep_initiation",    # Replication protein (Rep_1)
    "PF01719": "rep_initiation",    # Plasmid replication protein, origin binding (Rep_2)
    "PF01051": "rep_initiation",    # Initiator Replication protein, WH1 (Rep_3)
    "PF02486": "rep_initiation",    # Replication initiation factor (Rep_trans)
    "PF18008": "rep_initiation",    # Replication initiator protein A, C-terminal
    "PF06970": "rep_initiation",    # Replication initiator protein A (RepA), N-terminus
    "PF02195": "partition",         # ParB/Sulfiredoxin (ParB)
    "PF17762": "partition",         # HTH domain found in ParB protein
    "PF18607": "partition",         # ParA helix-turn-helix domain
    "PF01656": "partition",         # CobQ/CobB/MinD/ParA nucleotide-binding (ParA ATPase)
}

MARKER_SETS = {"viral": VIRAL_CLASSES, "plasmid": PLASMID_CLASSES}


def base_acc(acc: str) -> str:
    """Strip the Pfam version suffix, e.g. PF03354.14 -> PF03354."""
    return acc.split(".")[0]


def main():
    pfam_hmm, out_hmm = sys.argv[1], sys.argv[2]
    marker_set = sys.argv[3] if len(sys.argv) > 3 else "viral"
    classes = MARKER_SETS[marker_set]
    out_dir = os.path.dirname(os.path.abspath(out_hmm))
    os.makedirs(out_dir, exist_ok=True)

    wanted = set(classes)
    found = {}  # base_acc -> (pfam_name, has_ga)

    with open(pfam_hmm) as fin, open(out_hmm, "w") as fout:
        record = []
        rec_acc = None
        rec_name = None
        rec_has_ga = False
        in_record = False
        for line in fin:
            if line.startswith("HMMER3"):
                # start of a new record
                record = [line]
                rec_acc = rec_name = None
                rec_has_ga = False
                in_record = True
                continue
            if not in_record:
                continue
            record.append(line)
            if line.startswith("NAME "):
                rec_name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC "):
                rec_acc = base_acc(line.split(None, 1)[1].strip())
            elif line.startswith("GA "):
                rec_has_ga = True
            elif line.startswith("//"):
                if rec_acc in wanted and rec_acc not in found:
                    fout.writelines(record)
                    found[rec_acc] = (rec_name, rec_has_ga)
                in_record = False

    # marker_map.tsv
    map_path = os.path.join(out_dir, "marker_map.tsv")
    with open(map_path, "w") as m:
        m.write("accession\tclass\tpfam_name\thas_ga\n")
        for acc in classes:
            name, has_ga = found.get(acc, (None, None))
            status = "MISSING" if acc not in found else ""
            m.write(f"{acc}\t{classes[acc]}\t{name or status}\t{int(bool(has_ga))}\n")

    missing = wanted - set(found)
    print(f"[build_viral_hmm] extracted {len(found)}/{len(wanted)} models -> {out_hmm}")
    print(f"[build_viral_hmm] all have GA cutoffs: {all(v[1] for v in found.values())}")
    if missing:
        print(f"[build_viral_hmm] WARNING missing accessions: {sorted(missing)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
