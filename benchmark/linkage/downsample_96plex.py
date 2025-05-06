import os


def subsample(out_dir, full_bam, run_cmd):
    f = open(run_cmd, "w")
    for p in ["10", "20", "30", "50", "05"]:
        prefix = f"m64004_210929_143746.p{p}"
        out_bam = os.path.join(out_dir, f"m64004_210929_143746.p{p}.bam")
        cmd = f"""
            samtools view -s 42.{p} -b {full_bam} > {out_bam}
            /home/shuaiw//smrtlink/pbindex {out_bam}
            samtools index {out_bam}
        """
        # os.system(cmd)
        run = f"""
            /usr/bin/time -v -o {out_dir}/{prefix}.time python /home/shuaiw/Methy/main.py \
            --work_dir {out_dir}/{prefix} \
            --whole_bam {out_bam} \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --max_NM 3000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 30 \
            --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list
        """
        print (run, file = f)
    f.close()


full_bam = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.align.ccs.bam"
out_dir = "/home/shuaiw/borg/paper/linkage"
run_cmd = "run_dp.sh"
subsample(out_dir, full_bam, run_cmd)