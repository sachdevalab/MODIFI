import os


def subsample(out_dir, full_bam, run_cmd):
    f = open(run_cmd, "w")
    for p in ["100", "10", "20", "30", "50", "05"]:
        prefix = f"m64004_210929_143746.p{p}"
        out_bam = os.path.join(f"/home/shuaiw/borg/paper/linkage/m64004_210929_143746.p{p}.bam")
        cmd = f"""
            samtools view -s 42.{p} -b {full_bam} > {out_bam}
            /home/shuaiw//smrtlink/pbindex {out_bam}
            samtools index {out_bam}
        """
        # os.system(cmd)
        run = f"""
            sbatch  --partition standard --wrap " /usr/bin/time -v -o {out_dir}/{prefix}.time python /home/shuaiw/Methy/main.py \\
            --work_dir {out_dir}/{prefix} \\
            --whole_bam {out_bam} \\
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \\
            --read_type hifi \\
            --min_len 1000 \\
            --min_cov 1 \\
            --min_frac 0.3 \\
            --min_score 30 \\
            --min_sites 100 \\
            --min_iden 0.97 \\
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \\
            --threads 64 --run_steps host" \\
            --job-name=p{p}_pure
        """
        print (run, file = f)
    f.close()

def pure_main():
    full_bam = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.align.ccs.bam"
    out_dir = "/home/shuaiw/borg/paper/linkage/"
    out_dir = "/home/shuaiw/borg/paper/linkage/pure"
    run_cmd = "run_dp.sh"
    subsample(out_dir, full_bam, run_cmd)

def meta_subsample(out_dir, raw_96, soil):
    f = open("run_meta_dp.sh", "w")
    print ("#!/bin/bash\n#SBATCH --job-name=plex_soil \n #SBATCH --partition=standard", file = f)
    for p in ["100", "10", "20", "30", "50", "05"]:
        prefix = f"m64004_210929_143746.raw.p{p}"
        out_bam = os.path.join(out_dir, f"{prefix}.p{p}.bam")
        merge_bam = os.path.join(out_dir, f"{prefix}.soil.merge.bam")
        cmd = f"""
            samtools view -s 42.{p} -b {raw_96} > {out_bam}
            /home/shuaiw//smrtlink/pbindex {out_bam}
            samtools index {out_bam}
            samtools merge -f {merge_bam} {out_bam} {soil}
            /home/shuaiw//smrtlink/pbindex {merge_bam}
            samtools index {merge_bam}
        """
        # os.system(cmd)
        
        # alignment = f"""
        #     ~/smrtlink/pbmm2 align --preset CCS -j $SLURM_CPUS_ON_NODE /home/shuaiw/borg/contigs/soil_zymo.fa {merge_bam} {prefix}.align.raw.bam
        #     samtools sort -T {prefix} -@ $SLURM_CPUS_ON_NODE -o {prefix}.align.bam {prefix}.align.raw.bam
        #     rm {prefix}.align.raw.bam
        #     samtools index {prefix}.align.bam
        #     /home/shuaiw//smrtlink/pbindex {prefix}.align.bam
        # """
        # print (alignment, file = f)
        align_bam = f"/home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p{p}.align.bam"
        prefix = f"/home/shuaiw/borg/paper/linkage/meta2/m64004_210929_143746.p{p}"
        run = f"""
            sbatch  --partition standard --wrap "/usr/bin/time -v -o {prefix}.run.time python /home/shuaiw/Methy/main.py \\
            --work_dir {prefix}/ \\
            --whole_bam {align_bam} \\
            --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa \\
            --read_type hifi \\
            --min_len 1000 \\
            --min_cov 1 \\
            --min_frac 0.4 \\
            --min_score 30 \\
            --min_sites 100 \\
            --min_iden 0.97 \\
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \\
            --threads 64 --run_steps host" \\
            --job-name=p{p}
        """
        print (run, file = f)
    f.close()


soil = "/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam"
raw_96 = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.hifi_reads.bam"
out_dir = "/home/shuaiw/borg/paper/linkage/meta"
meta_subsample(out_dir, raw_96, soil)

# pure_main()

