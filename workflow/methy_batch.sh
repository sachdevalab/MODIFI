
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_3_INF1330004_5PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_3_INF1330004_5PB/NANO_3_INF1330004_5PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_3_INF1330004_5PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_3_INF1330004_5PB
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_3_INF1330004_7PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_3_INF1330004_7PB/NANO_3_INF1330004_7PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_3_INF1330004_7PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_3_INF1330004_7PB
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_3_INF1330004_8PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_3_INF1330004_8PB/NANO_3_INF1330004_8PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_3_INF1330004_8PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_3_INF1330004_8PB
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_2_INF1340011_2PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_2_INF1340011_2PB/NANO_2_INF1340011_2PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_2_INF1340011_2PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_2_INF1340011_2PB
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_2_INF1340011_3PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_2_INF1340011_3PB/NANO_2_INF1340011_3PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_2_INF1340011_3PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_2_INF1340011_3PB
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_3_INF1340011_6PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_3_INF1340011_6PB/NANO_3_INF1340011_6PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_3_INF1340011_6PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_3_INF1340011_6PB
        

        sbatch --partition standard --wrap "snakemake \
    --config whole_bam=/home/shuaiw/borg/allison/batch/NANO_3_INF1340011_7PB.align.ccs.bam \
    whole_ref=/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/NANO_3_INF1340011_7PB/NANO_3_INF1340011_7PB_HR_HIFIASM_META_scaffold_min1000.fa \
    work_dir=/home/shuaiw/borg/allison/batch/NANO_3_INF1340011_7PB_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
    --job-name=NANO_3_INF1340011_7PB
        
