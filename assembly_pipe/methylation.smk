# Snakefile

# Configuration
configfile: "config.yaml"

rule all:
    input:
        methy_finish = f"{config['work_dir']}/methylation.finish"


rule map_reads_to_assembly:
    input:
        fa = f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        bam = config["hifi_bam"],
        # checkm_finish = f"{config['work_dir']}/{config['prefix']}.checkm.finish"
    output:
        align_bam = temp(f"{config['work_dir']}/{config['prefix']}.align.bam"),
        bai=f"{config['work_dir']}/{config['prefix']}.align.bam.bai"
    threads: config["threads"]
    shell:
        """
        samtools faidx {input.fa}
        ~/smrtlink/pbmm2 align --preset CCS -j {threads} {input.fa} {input.bam} {config[work_dir]}/{config[prefix]}.raw.bam
        samtools sort -T {config[work_dir]}/{config[prefix]} -@ {threads} -o {output.align_bam} {config[work_dir]}/{config[prefix]}.raw.bam
        rm {config[work_dir]}/{config[prefix]}.raw.bam

        samtools index {output.align_bam}
        /home/shuaiw/smrtlink/pbindex {output.align_bam}
        """

rule call_methylation:
    input:
        bam=f"{config['work_dir']}/{config['prefix']}.align.bam",
        fa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        mge_file = f"{config['work_dir']}/all_mge.tsv",
    output:
        time=f"{config['work_dir']}/{config['prefix']}.methyl.time",
        summary=f"{config['work_dir']}/{config['prefix']}_methylation/summary.csv",
        methy_finish = f"{config['work_dir']}/methylation.finish"
    threads: config["threads"]
    shell:
        """
        /usr/bin/time -v -o {output.time} python /home/shuaiw/Methy/main.py \
          --work_dir {config[work_dir]}/{config[prefix]}_methylation \
          --whole_bam {input.bam} \
          --whole_ref {input.fa} \
          --read_type hifi \
          --min_len 1000 \
          --max_NM 3 \
          --min_cov 1 \
          --min_frac 0.4 \
          --min_score 30 \
          --min_sites 30 \
          --mge_file {input.mge_file} \
          --threads {threads} 
        touch {output.methy_finish}

        """