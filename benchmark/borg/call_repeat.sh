# python3 tandem_repeats.py \
#  -i /home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.fasta \
#  --min_repeat_count 3 -t 16 -o /home/shuaiw/borg/paper/borg_data/batch_export2/repeat \
#  -g /home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.gff


# python3 tandem_repeats.py \
#  -i /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.contigs.fa \
#  --min_repeat_count 3 -t 16 -o /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.repeat \
#  -g /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.prodigal.gff

# python3 tandem_repeats.py \
#  -i /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.contigs.fa \
#  --min_repeat_count 3 -t 16 -o /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.repeat \
#  -g /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.prodigal.gff


# python3 tandem_repeats.py \
#  -i /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.contigs.fa \
#  --min_repeat_count 3 -t 16 -o /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.repeat \
#  -g /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.prodigal.gff


# python3 tandem_repeats.py \
#  -i /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.fa \
#  --min_repeat_count 3 -t 16 -o /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.repeat \
#  -g /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.prodigal.gff


# ctg_name=SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L
# ref=/home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L.fa
# outfir=/home/shuaiw/borg/paper/borg_data/batch_export2/green_borgs
# prodigal -i $ref -f gff -o $outfir/$ctg_name.prodigal.gff
# python3 tandem_repeats.py \
#  -i $ref \
#  --min_repeat_count 3 -t 16 -o $outfir/$ctg_name.repeat \
#  -g $outfir/$ctg_name.prodigal.gff

# ctg_name=SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_METAMDBG_723659_L
# ref=/home/shuaiw/borg/paper/gg_run3/soil_100/soil_100_methylation4/contigs/SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_METAMDBG_723659_L.fa
# outfir=/home/shuaiw/borg/paper/borg_data/batch_export2/green_borgs
# prodigal -i $ref -f gff -o $outfir/$ctg_name.prodigal.gff
# python3 tandem_repeats.py \
#  -i $ref \
#  --min_repeat_count 3 -t 16 -o $outfir/$ctg_name.repeat \
#  -g $outfir/$ctg_name.prodigal.gff

#  ctg_name=SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_6165_L
# ref=/home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_6165_L.fa
# outfir=/home/shuaiw/borg/paper/borg_data/batch_export2/green_borgs
# prodigal -i $ref -f gff -o $outfir/$ctg_name.prodigal.gff
# python3 tandem_repeats.py \
#  -i $ref \
#  --min_repeat_count 3 -t 16 -o $outfir/$ctg_name.repeat \
#  -g $outfir/$ctg_name.prodigal.gff

for file in /home/shuaiw/borg/paper/borg_data/all_borgs/*.fa; do
    ctg_name=$(basename "$file" .fa)
    ref="$file"
    outdir=/home/shuaiw/borg/paper/borg_data/batch_export2/all_borgs
    
    prodigal -i "$ref" -f gff -o "$outdir/$ctg_name.prodigal.gff"
    python3 tandem_repeats.py \
        -i "$ref" \
        --min_repeat_count 3 -t 16 -o "$outdir/$ctg_name.repeat" \
        -g "$outdir/$ctg_name.prodigal.gff"
done
