# Snakefile

# Configuration
configfile: "config.yaml"

rule all:
    input:
        drep_97_finish = f"{config['work_dir']}/dRep97.finish",

rule dRep_97:
    input: f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa"
    output:
        drep_clu_file = f"{config['work_dir']}/dRep_out_97/data_tables/Cdb.csv",
        drep_97_finish = f"{config['work_dir']}/dRep_97.finish",
    params:
        genome_list = f"{config['work_dir']}/genome.list",
    threads: config["threads"]
    shell:
        """
        find {config[work_dir]}/bins/ -name "*.fasta" > {params.genome_list}
        dRep dereplicate \
        -p {threads} \
        -g {params.genome_list} \
        -comp 50 \
        -con 10 \
        --S_algorithm skani \
        -ms 10000 \
        -sa 0.97 \
        -nc 0.7 {config[work_dir]}/dRep_out_97
        touch {output.drep_97_finish}
        """