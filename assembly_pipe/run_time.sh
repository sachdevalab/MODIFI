# sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/96plex/time.txt python /home/shuaiw/mGlu/main.py \
#                 --work_dir /home/shuaiw/borg/paper/run2/96plex/96plex_methylation_time \
#                 --whole_bam /home/shuaiw/borg/paper/run2/96plex/96plex.align.bam \
#                 --whole_ref /home/shuaiw/borg/paper/run2/96plex/96plex.hifiasm.p_ctg.rename.fa \
#                 --read_type hifi \
#                 --clean \
#                 --threads 64" \
#                 --job-name=96plex
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/cow_bioreactor_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/cow_bioreactor_1/cow_bioreactor_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/cow_bioreactor_1/cow_bioreactor_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/cow_bioreactor_1/cow_bioreactor_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=cow_bioreactor_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/cow_bioreactor_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/cow_bioreactor_2/cow_bioreactor_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/cow_bioreactor_2/cow_bioreactor_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/cow_bioreactor_2/cow_bioreactor_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=cow_bioreactor_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/cow_bioreactor_3/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/cow_bioreactor_3/cow_bioreactor_3_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/cow_bioreactor_3/cow_bioreactor_3.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/cow_bioreactor_3/cow_bioreactor_3.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=cow_bioreactor_3
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/cow_bioreactor_4/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/cow_bioreactor_4/cow_bioreactor_4_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/cow_bioreactor_4/cow_bioreactor_4.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/cow_bioreactor_4/cow_bioreactor_4.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=cow_bioreactor_4
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/cow_bioreactor_5/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/cow_bioreactor_5/cow_bioreactor_5_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/cow_bioreactor_5/cow_bioreactor_5.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/cow_bioreactor_5/cow_bioreactor_5.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=cow_bioreactor_5
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/cow_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/cow_1/cow_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/cow_1/cow_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/cow_1/cow_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=cow_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/ERR12723528_mice/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=ERR12723528_mice
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/ERR12723529_mice/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/ERR12723529_mice/ERR12723529_mice_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/ERR12723529_mice/ERR12723529_mice.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/ERR12723529_mice/ERR12723529_mice.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=ERR12723529_mice
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/ERR5621427_sludge/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/ERR5621427_sludge/ERR5621427_sludge_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/ERR5621427_sludge/ERR5621427_sludge.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/ERR5621427_sludge/ERR5621427_sludge.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=ERR5621427_sludge
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/ERR5621429_sludge/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/ERR5621429_sludge/ERR5621429_sludge_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/ERR5621429_sludge/ERR5621429_sludge.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/ERR5621429_sludge/ERR5621429_sludge.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=ERR5621429_sludge
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/ERR5621430_sludge/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/ERR5621430_sludge/ERR5621430_sludge_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/ERR5621430_sludge/ERR5621430_sludge.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/ERR5621430_sludge/ERR5621430_sludge.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=ERR5621430_sludge
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_1/infant_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_1/infant_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_1/infant_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_2/infant_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_2/infant_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_2/infant_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_3/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_3/infant_3_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_3/infant_3.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_3/infant_3.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_3
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_4/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_4/infant_4_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_4/infant_4.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_4/infant_4.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_4
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_5/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_5/infant_5_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_5/infant_5.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_5/infant_5.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_5
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_6/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_6/infant_6_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_6/infant_6.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_6/infant_6.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_6
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_7/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_7/infant_7_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_7/infant_7.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_7/infant_7.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_7
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_8/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_8/infant_8_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_8/infant_8.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_8/infant_8.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_8
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_9/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_9/infant_9_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_9/infant_9.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_9/infant_9.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_9
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_10/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_10/infant_10_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_10/infant_10.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_10/infant_10.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_10
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_11/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_11/infant_11_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_11/infant_11.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_11/infant_11.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_11
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_12/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_12/infant_12_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_12/infant_12.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_12/infant_12.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_12
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_13/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_13/infant_13_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_13/infant_13.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_13/infant_13.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_13
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_14/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_14/infant_14_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_14/infant_14.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_14/infant_14.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_14
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_15/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_15/infant_15_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_15/infant_15.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_15/infant_15.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_15
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_16/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_16/infant_16_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_16/infant_16.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_16/infant_16.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_16
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_17/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_17/infant_17_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_17/infant_17.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_17/infant_17.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_17
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_18/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_18/infant_18_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_18/infant_18.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_18/infant_18.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_18
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_19/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_19/infant_19_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_19/infant_19.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_19/infant_19.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_19
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_20/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_20/infant_20_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_20/infant_20.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_20/infant_20.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_20
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_21/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_21/infant_21_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_21/infant_21.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_21/infant_21.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_21
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_22/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_22/infant_22_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_22/infant_22.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_22/infant_22.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_22
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_23/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_23/infant_23_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_23/infant_23.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_23/infant_23.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_23
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_24/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_24/infant_24_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_24/infant_24.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_24/infant_24.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_24
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_25/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_25/infant_25_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_25/infant_25.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_25/infant_25.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_25
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_26/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_26/infant_26_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_26/infant_26.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_26/infant_26.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_26
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_27/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_27/infant_27_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_27/infant_27.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_27/infant_27.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_27
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/infant_28/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/infant_28/infant_28_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/infant_28/infant_28.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/infant_28/infant_28.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=infant_28
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/ocean_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/ocean_1/ocean_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/ocean_1/ocean_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/ocean_1/ocean_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=ocean_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_s3_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_s3_1/soil_s3_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_s3_1/soil_s3_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_s3_1/soil_s3_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_s3_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_s3_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_s3_2/soil_s3_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_s3_2/soil_s3_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_s3_2/soil_s3_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_s3_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_s4_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_s4_1/soil_s4_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_s4_1/soil_s4_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_s4_1/soil_s4_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_s4_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_s4_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_s4_2/soil_s4_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_s4_2/soil_s4_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_s4_2/soil_s4_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_s4_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_s1_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_s1_1/soil_s1_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_s1_1/soil_s1_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_s1_1/soil_s1_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_s1_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_s1_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_s1_2/soil_s1_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_s1_2/soil_s1_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_s1_2/soil_s1_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_s1_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_1/soil_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_1/soil_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/soil_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/soil_2/soil_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/soil_2/soil_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/soil_2/soil_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=soil_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/SRR14074352_human/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/SRR14074352_human/SRR14074352_human_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/SRR14074352_human/SRR14074352_human.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/SRR14074352_human/SRR14074352_human.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=SRR14074352_human
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/SRR23446539_sugarcane/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/SRR23446539_sugarcane/SRR23446539_sugarcane_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/SRR23446539_sugarcane/SRR23446539_sugarcane.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/SRR23446539_sugarcane/SRR23446539_sugarcane.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=SRR23446539_sugarcane
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/SRR23446540_sugarcane/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/SRR23446540_sugarcane/SRR23446540_sugarcane_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/SRR23446540_sugarcane/SRR23446540_sugarcane.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/SRR23446540_sugarcane/SRR23446540_sugarcane.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=SRR23446540_sugarcane
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_1/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_1/asthma_1_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_1/asthma_1.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_1/asthma_1.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_1
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_2/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_2/asthma_2_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_2/asthma_2.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_2/asthma_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_2
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_3/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_3/asthma_3_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_3/asthma_3.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_3/asthma_3.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_3
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_4/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_4/asthma_4_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_4/asthma_4.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_4/asthma_4.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_4
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_5/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_5/asthma_5_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_5/asthma_5.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_5/asthma_5.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_5
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_6/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_6/asthma_6_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_6/asthma_6.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_6/asthma_6.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_6
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_7/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_7/asthma_7_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_7/asthma_7.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_7/asthma_7.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_7
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_8/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_8/asthma_8_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_8/asthma_8.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_8/asthma_8.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_8
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_9/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_9/asthma_9_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_9/asthma_9.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_9/asthma_9.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_9
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_10/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_10/asthma_10_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_10/asthma_10.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_10/asthma_10.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_10
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_11/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_11/asthma_11_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_11/asthma_11.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_11/asthma_11.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_11
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_12/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_12/asthma_12_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_12/asthma_12.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_12/asthma_12.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_12
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_13/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_13/asthma_13_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_13/asthma_13.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_13/asthma_13.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_13
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_14/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_14/asthma_14_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_14/asthma_14.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_14/asthma_14.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_14
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_15/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_15/asthma_15_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_15/asthma_15.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_15/asthma_15.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_15
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_16/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_16/asthma_16_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_16/asthma_16.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_16/asthma_16.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_16
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_17/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_17/asthma_17_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_17/asthma_17.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_17/asthma_17.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_17
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_18/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_18/asthma_18_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_18/asthma_18.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_18/asthma_18.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_18
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_19/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_19/asthma_19_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_19/asthma_19.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_19/asthma_19.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_19
sbatch  --partition standard --wrap "/usr/bin/time -v -o  /home/shuaiw/borg/paper/run2/asthma_20/time.txt python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2/asthma_20/asthma_20_methylation_time \
                --whole_bam /home/shuaiw/borg/paper/run2/asthma_20/asthma_20.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2/asthma_20/asthma_20.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --clean \
                --threads 64" \
                --job-name=asthma_20
