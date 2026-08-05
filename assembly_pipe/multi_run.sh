#!/bin/bash

REF=/home/shuaiw/borg/paper/run2/infant_8/infant_8.hifiasm.p_ctg.rename.fa
WORK_BASE=/home/shuaiw/methylation/binning/cross_sample/infant_8

# sbatch --partition standard --job-name=infant_1 --wrap "python /home/shuaiw/MODIFI/main.py \
#         -o ${WORK_BASE}/infant_1 \
#         -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240037_3G1_pacbio.bam \
#         -r ${REF} \
#         --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_2 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_2 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_2G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_3 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_3 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_3G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_4 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_4 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_4G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_5 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_5 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_2G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_6 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_6 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_3G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_7 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_7 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330001_4G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_8 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_8 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330004_3G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_9 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_9 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330004_4G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_10 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_10 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_2G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_11 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_11 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_3G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_12 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_12 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_4G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_13 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_13 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_2_MAT1240037_2G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_14 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_14 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_5G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_15 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_15 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_6G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_16 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_16 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1310001_8G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_17 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_17 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1310007_7G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_18 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_18 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_5G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_19 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_19 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_7G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_20 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_20 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_8G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_21 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_21 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340008_8G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_22 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_22 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_6G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_23 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_23 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_7G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_24 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_24 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240043_5G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_25 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_25 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240045_6G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_26 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_26 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240119_7G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_27 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_27 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_4_INF1330101_8G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"

sbatch --partition standard --job-name=infant_28 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o ${WORK_BASE}/infant_28 \
        -b /home/shuaiw/borg/paper/aws/infant/NANO_4_INF1340021_5G1_pacbio.bam \
        -r ${REF} \
        --read_type hifi --threads 64"
