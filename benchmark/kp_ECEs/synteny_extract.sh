#!/bin/bash
# Extract the 18 members of cluster infant_15_35_C into per-member FASTAs + a combined FASTA,
# for the LoVis4u synteny plot. Sequences come from the aggregated ECE fasta (no run2 re-extraction).
set -euo pipefail
SEQKIT="conda run -n tldr seqkit"
SRC=/home/shuaiw/borg/revision/kp_eces/seqs/all_metagenome_eces.fna
OUT=/home/shuaiw/borg/revision/kp_eces/synteny
mkdir -p "$OUT/contigs"

MEMBERS="infant_15_35_C asthma_20_19_C infant_4_322_L infant_28_15_C infant_1_8_C infant_4_48_L \
infant_23_17_C asthma_20_24_C infant_1_15_C infant_12_5_C infant_25_33_L infant_9_12_L \
infant_8_26_C infant_20_17_C infant_27_17_C asthma_6_14_C infant_25_149_L infant_3_3_C"

printf "%s\n" $MEMBERS > "$OUT/members.txt"
$SEQKIT grep -f "$OUT/members.txt" "$SRC" > "$OUT/members.fna"

# some members are absent from the aggregated fasta -> pull them from their run2 assembly by name
RUN2=/groups/banfield/projects/multienv/methylation/data/borg/paper/run2
found=$($SEQKIT seq -n "$OUT/members.fna")
for m in $MEMBERS; do
  if ! grep -qx "$m" <<< "$found"; then
    samp=$(echo "$m" | sed -E 's/_[0-9]+_[CL]$//')
    asm="$RUN2/$samp/$samp.hifiasm.p_ctg.rename.fa"
    echo "recovering $m from $asm"
    $SEQKIT grep -p "$m" "$asm" >> "$OUT/members.fna"
  fi
done

# per-member fasta
for m in $MEMBERS; do
  $SEQKIT grep -p "$m" "$OUT/members.fna" > "$OUT/contigs/$m.fna"
done

$SEQKIT fx2tab -nl "$OUT/members.fna" > "$OUT/members_len.tsv"

n=$($SEQKIT stats -T "$OUT/members.fna" | awk 'NR==2{print $4}')
echo "extracted sequences: $n (expected 18)"
echo "--- lengths ---"; cat "$OUT/members_len.tsv"
