## remove output directory if exists
rm -rf output/

python ../main.py --whole_bam ERR10042285_2_L.bam \
--work_dir output/ \
--whole_ref ERR10042285_2_L.fa \
--read_type hifi \
--kmer_mean_db ../control_db/control_db.up7.down3.mean.dat \
--kmer_num_db ../control_db/control_db.up7.down3.num.dat \
--threads 5