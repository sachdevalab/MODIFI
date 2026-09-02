#!/bin/bash
# Add fibertools to the compute-resource scaling benchmark (Figure S2).
# fibertools is HiFi-only, so we REUSE MODIFI's existing HiFi alignment
# (modifi_hifi_full/align_bam/aligned.bam = soil CCS aligned to the full SR-VP
# assembly) and subset it to each test_N reference's contigs (index-based, no
# re-alignment), then run ft predict-m6a + ft pileup, timed like MODIFI/ipdSummary.
# Runs in the terminal with 64 threads. Usage: batch_fibertools.sh [test_100 ...]
set -uo pipefail

THREADS=64
SAMTOOLS=/shared/software/bin/samtools
PBINDEX=/home/shuaiw/smrtlink/pbindex
FT=/home/shuaiw/miniconda3/envs/fibertools/bin/ft
REF_DIR=/home/shuaiw/borg/contigs
ALIGNED=/home/shuaiw/borg/paper/ipdsummary/compare_all_meta/modifi_hifi_full/align_bam/aligned.bam
OUTDIR=/home/shuaiw/borg/paper/ipdsummary/soil_1
FT_OUT=$OUTDIR/fibertools.out
mkdir -p "$FT_OUT"

REFS=("$@"); [ ${#REFS[@]} -eq 0 ] && REFS=(test_100 test_200 test_300 test_500)

for base in "${REFS[@]}"; do
    ref=$REF_DIR/$base.fa
    ccs_bam=$OUTDIR/$base.align.ccs.bam
    m6a_bam=$OUTDIR/$base.m6a.bam
    bed=$OUTDIR/$base.ccs.bed
    log=$FT_OUT/$base.log
    echo "=== [$base] start $(date) ===" | tee "$log"

    # 1) subset the existing HiFi alignment to this reference's contigs
    if [ ! -s "$ccs_bam" ]; then
        [ -s "$ref.fai" ] || $SAMTOOLS faidx "$ref"
        awk 'BEGIN{OFS="\t"}{print $1,0,$2}' "$ref.fai" > "$bed"
        echo "[$base] contigs=$(wc -l < "$bed"); subsetting aligned.bam ..." | tee -a "$log"
        $SAMTOOLS view -@ $THREADS -b -M -L "$bed" "$ALIGNED" -o "$ccs_bam" >>"$log" 2>&1
        $SAMTOOLS index -@ $THREADS "$ccs_bam" >>"$log" 2>&1
        $PBINDEX "$ccs_bam" >>"$log" 2>&1
    fi
    echo "[$base] ccs reads=$($SAMTOOLS view -c -@ $THREADS "$ccs_bam" 2>/dev/null)" | tee -a "$log"

    # 2) fibertools predict-m6a (timed) -> pileup (timed); alignment/subset excluded from timing
    /usr/bin/time -v -o "$FT_OUT/$base.ft_predict.time" \
        $FT predict-m6a -t $THREADS "$ccs_bam" "$m6a_bam" >>"$log" 2>&1
    echo "[$base] predict-m6a done $(date)" | tee -a "$log"
    $SAMTOOLS index -@ $THREADS "$m6a_bam" >>"$log" 2>&1
    /usr/bin/time -v -o "$FT_OUT/$base.ft_pileup.time" \
        $FT pileup -m --per-base -o "$FT_OUT/$base.pileup.bed" "$m6a_bam" >>"$log" 2>&1
    echo "[$base] pileup done $(date); pileup lines=$(wc -l < "$FT_OUT/$base.pileup.bed" 2>/dev/null)" | tee -a "$log"

    pcpu=$(grep -E 'User time|System time' "$FT_OUT/$base.ft_predict.time" 2>/dev/null | awk '{s+=$NF}END{print s}')
    echo "=== [$base] DONE $(date) ===" | tee -a "$log"
done
