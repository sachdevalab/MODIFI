#!/bin/bash
#SBATCH --job-name=plex_infant 
 #SBATCH --partition=standard

            ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.raw.p100.infant.merge.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.raw.bam
            samtools sort -T /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100 -@ $SLURM_CPUS_ON_NODE -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.raw.bam
            rm /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.raw.bam
            samtools index /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.bam
            /home/shuaiw//smrtlink/pbindex /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.bam
        

            ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.raw.p10.infant.merge.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.raw.bam
            samtools sort -T /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10 -@ $SLURM_CPUS_ON_NODE -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.raw.bam
            rm /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.raw.bam
            samtools index /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.bam
            /home/shuaiw//smrtlink/pbindex /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.bam
        

            ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.raw.p20.infant.merge.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.raw.bam
            samtools sort -T /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20 -@ $SLURM_CPUS_ON_NODE -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.raw.bam
            rm /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.raw.bam
            samtools index /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.bam
            /home/shuaiw//smrtlink/pbindex /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.bam
        

            ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.raw.p30.infant.merge.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.raw.bam
            samtools sort -T /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30 -@ $SLURM_CPUS_ON_NODE -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.raw.bam
            rm /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.raw.bam
            samtools index /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.bam
            /home/shuaiw//smrtlink/pbindex /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.bam
        

            ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.raw.p50.infant.merge.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.raw.bam
            samtools sort -T /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50 -@ $SLURM_CPUS_ON_NODE -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.raw.bam
            rm /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.raw.bam
            samtools index /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.bam
            /home/shuaiw//smrtlink/pbindex /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.bam
        

            ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.raw.p05.infant.merge.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.raw.bam
            samtools sort -T /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05 -@ $SLURM_CPUS_ON_NODE -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.raw.bam
            rm /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.raw.bam
            samtools index /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.bam
            /home/shuaiw//smrtlink/pbindex /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.bam
        
