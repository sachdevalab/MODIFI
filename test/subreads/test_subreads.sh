## remove output directory if exists
rm -rf output/

python ../../main.py --aligned_bam SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_100_L.bam \
-o output/ \
-r SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_100_L.fa \
--read_type subreads \
--kmer_mean_db ../../control_db/control_db.up7.down3.mean.dat \
--kmer_num_db ../../control_db/control_db.up7.down3.num.dat \
--threads 5