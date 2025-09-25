# Snakefile

# Configuration
configfile: "config.yaml"

rule all:
    input:
        checkm=f"{config['work_dir']}/checkM2_2/quality_report.tsv",
        gtdb_finish=f"{config['work_dir']}/GTDB_2/gtdbtk.done",

rule checkM:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        checkm=f"{config['work_dir']}/checkM2_2/quality_report.tsv",
        finish=f"{config['work_dir']}/{config['prefix']}.checkm.finish"
    threads: config["threads"]
    shell:
        """ 
        mkdir {config[work_dir]}/bins2/
        python split_ctgs.py {input.fasta} {config[work_dir]}/bins2/
        checkm2 predict --input {config[work_dir]}/bins2/ --output-directory  {config[work_dir]}/checkM2_2 --force -x .fasta --threads {threads}
        touch {output.finish}
        """

rule GTDB:
    input:
        checkm=f"{config['work_dir']}/checkM2_2/quality_report.tsv",
    params:
        bin_dir=f"{config['work_dir']}/bins2/",
    output:
        gtdb_finish=f"{config['work_dir']}/GTDB_2/gtdbtk.done"
    threads: config["threads"]
    shell:
        """
        gtdbtk classify_wf \
          --genome_dir {params.bin_dir} \
          --out_dir {config[work_dir]}/GTDB_2 \
          --skip_ani_screen \
          --cpus {threads} \
          -x fasta
        touch {output.gtdb_finish}
        """

