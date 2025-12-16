# Snakefile

# Configuration
configfile: "config_soil.yaml"

rule all_annotation:
    input:
        ctg_mge = f"{config['work_dir']}/ctg_mge.done",
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation4/host_summary.csv",


rule get_ctg_mge:
    input:
        genomad_finish=f"{config['work_dir']}/genomad.done"
    output:
        mge_finish = f"{config['work_dir']}/ctg_mge.done",
        mge_file = f"{config['work_dir']}/all_mge2.tsv",
        host_file = f"{config['work_dir']}/all_host_ctgs.tsv"
    params:
        output_dir = config['work_dir'],
        prefix = config["prefix"],
        ref=config['ref']
    shell:
        """
        python merge_MGEs.py {params.output_dir} {params.prefix} {params.ref}
        touch {output.mge_finish}
        """

rule call_host:
    input:
        fa=f"{config['ref']}",
        mge_file = f"{config['work_dir']}/all_mge2.tsv",
    output:
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation4/host_summary.csv"
    threads: config["threads"]
    shell:
        """
        /home/shuaiw/miniconda3/envs/methy3/bin/python /home/shuaiw/mGlu/main.py \
          --work_dir {config[work_dir]}/{config[prefix]}_methylation4 \
          --unaligned_bam {config[hifi_bam]} \
          --whole_ref {input.fa} \
          --read_type hifi \
            --min_len 1000 \
            --min_cov 3 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --run_steps host \
          --mge_file {input.mge_file} \
          --threads {threads} 
        """