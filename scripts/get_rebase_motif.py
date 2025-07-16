import re

output = []
unique_raw_motifs = {}

with open("rebase/withrefm.txt") as f:
    for line in f:
        if line.startswith("<3>"):
            raw = line.strip()[3:]  # Remove "<3>"
            if not raw:
                continue
            if raw in unique_raw_motifs:
                continue
            unique_raw_motifs[raw] = True
            # Clean motif (remove ^)
            cleaned = raw.replace("^", "")
            # Find methylation site (position of ^)
            meth_pos = raw.find("^")
            output.append((raw, cleaned, meth_pos + 1))  # +1 for 1-based

# Save as TSV
with open("rebase/motifs_cleaned.tsv", "w") as out:
    out.write("RawMotif\tCleanedMotif\tMethylationPosition\n")
    for row in output:
        out.write("\t".join(map(str, row)) + "\n")

print(f"Extracted {len(output)} motifs with positions.")

