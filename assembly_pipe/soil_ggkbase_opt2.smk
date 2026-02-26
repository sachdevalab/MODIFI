# Snakefile

# Configuration
configfile: "config_soil.yaml"

rule all_annotation:
    input:
        methy_finish=f"{config['work_dir']}/{config['prefix']}.methy.finish",


rule call_methy:
    input:
        fa=f"{config['ref']}",
    output:
        methy_finish=f"{config['work_dir']}/{config['prefix']}.methy.finish",
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
            --run_steps anno \
          --threads {threads} 
        touch {output.methy_finish}
        """



# rule GTDB:
#     input:
#         fasta=f"{config['ref']}",
#         checkv_finish=f"{config['work_dir']}/checkV.done",
#         methy_finish=f"{config['work_dir']}/{config['prefix']}.methy.finish",
#     params:
#         bin_dir=f"{config['work_dir']}/bins/",
#     output:
#         gtdb_finish=f"{config['work_dir']}/GTDB/gtdbtk.done"
#     threads: config["threads"]
#     shell:
#         """
#         python split_ctgs.py {input.fasta} {params.bin_dir}
#         gtdbtk classify_wf \
#           --genome_dir {params.bin_dir} \
#           --out_dir {config[work_dir]}/GTDB \
#           --skip_ani_screen \
#           --cpus {threads} \
#           -x fasta
#         touch {output.gtdb_finish}
#         """


