# Snakefile

# Configuration
configfile: "config_soil.yaml"

rule all_annotation:
    input:
        checkm=f"{config['work_dir']}/checkM2/quality_report.tsv",
        gtdb_finish=f"{config['work_dir']}/GTDB/gtdbtk.done",
        genomad_finish=f"{config['work_dir']}/genomad.done",
        ctg_mge = f"{config['work_dir']}/ctg_mge.done",
        checkv_finish=f"{config['work_dir']}/checkV.done",
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation3/host_summary.csv",


rule checkM:
    input:
        fasta=f"{config['ref']}",
    output:
        checkm=f"{config['work_dir']}/checkM2/quality_report.tsv",
        finish=f"{config['work_dir']}/{config['prefix']}.checkm.finish"
    threads: config["threads"]
    shell:
        """ 
        python split_ctgs.py {input.fasta} {config[work_dir]}/bins/
        checkm2 predict --input {config[work_dir]}/bins/ --output-directory  {config[work_dir]}/checkM2 --force -x .fasta --threads {threads}
        touch {output.finish}
        """

rule checkV:
    input:
        checkm=f"{config['work_dir']}/checkM2/quality_report.tsv",
        fasta=f"{config['ref']}",
    output:
        checkv_finish=f"{config['work_dir']}/checkV.done"
    threads: config["threads"]
    shell:
        """
        checkv end_to_end \
            {input.fasta} \
            {config[work_dir]}/checkV \
            -t {threads} \
            --remove_tmp
        touch {output.checkv_finish}
        """

rule GTDB:
    input:
        checkv_finish=f"{config['work_dir']}/checkV.done",
    params:
        bin_dir=f"{config['work_dir']}/bins/",
    output:
        gtdb_finish=f"{config['work_dir']}/GTDB/gtdbtk.done"
    threads: config["threads"]
    shell:
        """
        gtdbtk classify_wf \
          --genome_dir {params.bin_dir} \
          --out_dir {config[work_dir]}/GTDB \
          --skip_ani_screen \
          --cpus {threads} \
          -x fasta
        touch {output.gtdb_finish}
        """

rule genomad:
    input:
        fasta=f"{config['ref']}",
        gtdb_finish=f"{config['work_dir']}/GTDB/gtdbtk.done"
    output:
        genomad_finish=f"{config['work_dir']}/genomad.done"
    threads: config["threads"]
    shell:
        """
        genomad end-to-end --relaxed --cleanup --enable-score-calibration \
            --threads {threads} --sensitivity 7.0 --force-auto \
            {input.fasta} \
            {config[work_dir]}/Genomad/ \
            /groups/diamond/databases/genomad/v1.7/
        touch {output.genomad_finish}
        """

rule get_ctg_mge:
    input:
        # prokka_finish = f"{config['work_dir']}/prokka/prokka.finish",
        # vibrantr=f"{config['work_dir']}/vibrant.done",
        genomad_finish=f"{config['work_dir']}/genomad.done"
    output:
        mge_finish = f"{config['work_dir']}/ctg_mge.done",
        mge_file = f"{config['work_dir']}/all_mge.tsv",
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
        fa=f"{config['ref']}",
        mge_file = f"{config['work_dir']}/all_mge.tsv",
    output:
        host_summary = f"{config['work_dir']}/{config['prefix']}_methylation3/host_summary.csv"
    threads: config["threads"]
    shell:
        """
        /home/shuaiw/miniconda3/envs/methy3/bin/python /home/shuaiw/mGlu/main.py \
          --work_dir {config[work_dir]}/{config[prefix]}_methylation3 \
          --unaligned_bam {config[hifi_bam]} \
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