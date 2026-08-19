#!/bin/bash
# Download the 8 JF8 reference genomes (from ref_ID.table) and map Mock_JF8.fa contigs onto them
# to assign each JF8 contig to a species. Downloads need internet -> run on the login node.
set -euo pipefail

REFDIR=/home/shuaiw/borg/revision/motif_class/jf8_refs
OUT=/home/shuaiw/borg/revision/motif_class
MOCK=/home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa
mkdir -p "$REFDIR"
cd "$REFDIR"

# species<TAB>accession  (from ref_ID.table). Underscores in species -> use a safe tag.
map_accession() {
  # $1 = accession ; echo the species tag
  :
}

# species tag | accession | type(gca/nuc)
cat > refs.list <<'EOF'
Bacteroides_caccae	GCA_000169015	gca
Bacteroides_ovatus	NZ_CP012938	nuc
Bacteroides_thetaiotaomicron	NC_004663	nuc
Bacteroides_vulgatus	NC_009614.1	nuc
Collinsella_aerofaciens	GCA_000169035	gca
Clostridium_bolteae	GCA_000154365	gca
Escherichia_coli	NC_000913	nuc
Ruminococcus_gnavus	GCA_000169475	gca
EOF

: > combined_ref.fa
while IFS=$'\t' read -r sp acc typ; do
  [ -z "$sp" ] && continue
  echo ">>> $sp  $acc  ($typ)"
  if [ "$typ" = "gca" ]; then
    pre=$(echo "$acc" | sed -E 's/GCA_([0-9]{3})([0-9]{3})([0-9]{3}).*/GCA\/\1\/\2\/\3/')
    base="https://ftp.ncbi.nlm.nih.gov/genomes/all/${pre}/"
    sub=$(curl -s -L "$base" | grep -oE "${acc}[^/\"]*" | head -1)
    curl -s -L "${base}${sub}/${sub}_genomic.fna.gz" -o "${sp}.fna.gz"
    zcat "${sp}.fna.gz" > "${sp}.fna"
  else
    efetch -db nuccore -id "$acc" -format fasta > "${sp}.fna"
  fi
  # relabel every contig header with the species tag so we can read species off the target name
  seqkit replace -p '^' -r "${sp}|" "${sp}.fna" >> combined_ref.fa
  echo "    $(grep -c '>' ${sp}.fna) sequences"
done < refs.list

echo ">>> mapping Mock_JF8.fa onto combined reference"
minimap2 -x asm20 -t 8 --secondary=no combined_ref.fa "$MOCK" > jf8_vs_refs.paf 2>/dev/null

# assign each JF8 contig to the species with the most aligned bases
python3 - "$MOCK" jf8_vs_refs.paf "$OUT/jf8_contig_species.tsv" <<'PY'
import sys, csv
from collections import defaultdict
mock, paf, out = sys.argv[1], sys.argv[2], sys.argv[3]
SPECIES = {
 "Bacteroides_caccae":"Bacteroides caccae","Bacteroides_ovatus":"Bacteroides ovatus",
 "Bacteroides_thetaiotaomicron":"Bacteroides thetaiotaomicron","Bacteroides_vulgatus":"Bacteroides vulgatus",
 "Collinsella_aerofaciens":"Collinsella aerofaciens","Clostridium_bolteae":"Clostridium bolteae",
 "Escherichia_coli":"Escherichia coli","Ruminococcus_gnavus":"Ruminococcus gnavus",
}
aln = defaultdict(lambda: defaultdict(int))   # query -> species_tag -> matching bases
for line in open(paf):
    f = line.rstrip("\n").split("\t")
    if len(f) < 12: continue
    q = f[0]; tgt = f[5]; nmatch = int(f[9])   # residue matches
    sp_tag = tgt.split("|", 1)[0]
    aln[q][sp_tag] += nmatch
# all query names (include unmapped contigs)
names = []
for rec in open(mock):
    if rec.startswith(">"):
        names.append(rec[1:].split()[0])
rows = []
for q in names:
    if q in aln and aln[q]:
        best = max(aln[q].items(), key=lambda kv: kv[1])
        rows.append((q, SPECIES.get(best[0], best[0]), best[1]))
    else:
        rows.append((q, "unassigned", 0))
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    for r in rows: w.writerow(r)
from collections import Counter
c = Counter(r[1] for r in rows)
sys.stderr.write(f"[map] {len(rows)} contigs -> {out}\n")
for k,v in sorted(c.items(), key=lambda x:-x[1]):
    sys.stderr.write(f"  {v:4} {k}\n")
PY

rm -f test_gca.zip gca_test.fna.gz
echo ">>> done: $OUT/jf8_contig_species.tsv"
