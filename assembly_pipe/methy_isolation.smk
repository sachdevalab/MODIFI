# Snakefile

# Configuration
configfile: "config.yaml"

rule all:
    input:
        # methy_finish = f"{config['work_dir']}/{config['prefix']}_methylation2/methylation.finish"
        rm_finish = f"{config['work_dir']}/{config['prefix']}_methylation2/RM_systems/all_ctgs_RM.rm.genes.tsv"


rule call_methylation:
    input:
        bam=f"{config['work_dir']}/{config['prefix']}.align.bam",
        fa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        mge_file = f"{config['work_dir']}/all_mge.tsv",
    output:
        # time = f"{config['work_dir']}/{config['prefix']}_methylation2/methylation.time",
        # methy_finish = f"{config['work_dir']}/{config['prefix']}_methylation2/methylation.finish",
        rm_finish = f"{config['work_dir']}/{config['prefix']}_methylation2/RM_systems/all_ctgs_RM.rm.genes.tsv"

    threads: config["threads"]
    # conda: "methy3"
    shell:
        """
        /home/shuaiw/miniconda3/envs/methy3/bin/python /home/shuaiw/mGlu/main.py \
          --work_dir {config[work_dir]}/{config[prefix]}_methylation2 \
          --whole_bam {input.bam} \
          --whole_ref {input.fa} \
          --read_type hifi \
          --min_len 1000 \
          --min_cov 10 \
          --min_frac 0.1 \
          --min_score 30 \
          --min_sites 10 \
          --annotate_rm \
          --run_steps anno \
          --mge_file {input.mge_file} \
          --threads 64 \
          --kmer_mean_db /home/shuaiw/mGlu/control_db/control_db.up7.down3.mean.dat \
          --kmer_num_db /home/shuaiw/mGlu/control_db/control_db.up7.down3.num.dat
        """

