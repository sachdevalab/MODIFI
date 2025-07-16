# Snakefile

# Configuration
configfile: "config.yaml"

host_ref = config["host_ref"]

rule all_host:
    input:
        mge_file = f"{config['work_dir']}/all_mge.tsv",
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation/host_summary.csv",
        methyl_time = f"{config['work_dir']}/{config['prefix']}.methyl.time",

rule call_host:
    input:
        bam=f"{config['work_dir']}/{config['prefix']}.align.bam",
        fa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        mge_file = f"{config['work_dir']}/all_mge.tsv",
    output:
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation/host_summary.csv"
    threads: config["threads"]
    shell:
        """
        python /home/shuaiw/Methy/main.py \
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
          --run_steps host \
          --mge_file {input.mge_file} \
          --threads {threads} 
        """