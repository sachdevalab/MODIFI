# Snakefile

# Configuration
configfile: "config.yaml"

rule all:
    input:
        mge_file = f"{config['work_dir']}/all_mge2.tsv",
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation4/host_summary.csv",




rule get_ctg_mge:
    input:
    output:
        mge_finish = f"{config['work_dir']}/ctg_mge.done",
        mge_file = f"{config['work_dir']}/all_mge2.tsv",
        host_file = f"{config['work_dir']}/all_host_ctgs.tsv"
    params:
        output_dir = config['work_dir'],
        prefix = config["prefix"]
    shell:
        """
        python merge_MGEs.py {params.output_dir} {params.prefix}
        touch {output.mge_finish}
        """

rule call_host:
    input:
        bam=f"{config['work_dir']}/{config['prefix']}.align.bam",
        fa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        mge_file = f"{config['work_dir']}/all_mge2.tsv",
    output:
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation4/host_summary.csv"
    threads: config["threads"]
    shell:
        """
        /home/shuaiw/miniconda3/envs/methy3/bin/python /home/shuaiw/mGlu/main.py \
          --work_dir {config[work_dir]}/{config[prefix]}_methylation4 \
          --whole_bam {input.bam} \
          --whole_ref {input.fa} \
          --read_type hifi \
          --min_len 1000 \
          --min_cov 3 \
          --min_frac 0.3 \
          --min_score 30 \
          --min_sites 100 \
          --mge_file {input.mge_file} \
          --threads {threads} 
        """