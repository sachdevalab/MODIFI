## remove output directory if exists
rm -rf output/

python ../main.py --whole_bam infant_1_29_C.bam \
--work_dir output/ \
--whole_ref infant_1_29_C.fa \
--read_type hifi 