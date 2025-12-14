# Snakefile

# Configuration
configfile: "config_soil.yaml"

rule all_annotation:
    input:
        gtdb_finish=f"{config['work_dir']}/GTDB/gtdbtk.done",

rule GTDB:
    input:
        
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

