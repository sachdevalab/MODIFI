# Snakefile

# Configuration
configfile: "config.yaml"

host_ref = config["host_ref"]

rule all_assembly:
    input:
        final_fa = f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        bam_index = f"{config['work_dir']}/{config['prefix']}.align.bam.bai",
        methyl_time = f"{config['work_dir']}/{config['prefix']}.methyl.time",
        methy_finish = f"{config['work_dir']}/methylation.finish",

rule bam_to_fastq_gz:
    input:
        bam=config["hifi_bam"]
    output:
        fastq_gz=f"{config['work_dir']}/{config['prefix']}.hifi.fastq.gz"
    params:
        prefix  = config["prefix"],
        work_dir = config["work_dir"]
    threads: config["threads"]
    shell:
        """
        mkdir -p {config[work_dir]}
        /home/shuaiw/smrtlink/pbindex {input.bam}
        /home/shuaiw/smrtlink/bam2fastq {input.bam} -j {threads} -o {params.work_dir}/{params.prefix}
        mv {params.work_dir}/{params.prefix}.fastq.gz {output.fastq_gz}
        """
        ### transform bam to fastq, bam2fastq  /home/shuaiw/smrtlink/bam2fastq

rule trim_reads_bbduk:
    input:
        fastq_gz=f"{config['work_dir']}/{config['prefix']}.hifi.fastq.gz",
    output:
        trimmed=f"{config['work_dir']}/{config['prefix']}.hifi.qc.fq.gz"
    threads: config["threads"]
    shell:
        """
        /shared/software/bbmap/39.01/bbduk.sh \
            in={input.fastq_gz} \
            out={output.trimmed} \
            minavgquality=20 \
            qtrim=rl \
            trimq=20 \
            threads={threads}
        """

rule assemble_hifiasm:
    input:
        fq=f"{config['work_dir']}/{config['prefix']}.hifi.qc.fq.gz"
    output:
        gfa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.gfa",
    log: f"{config['work_dir']}/{config['prefix']}.hifiasm.log"
    threads: config["threads"]
    shell:
        """
        hifiasm --primary -l0 -o {config[work_dir]}/{config[prefix]}.hifiasm \
            -t {threads} {input.fq} > {log} 
        """

rule gfa_to_fasta:
    input:
        gfa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.gfa"
    output:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa"
    shell:
        # Rohan's version:
        """
        awk '$1==\"S\" {{printf \">%s\\n%s\\n\", $2, $3}}' {input.gfa} > {config[work_dir]}/{config[prefix]}.p_ctg.fa
        python /home/shuaiw/Methy/assembly_pipe/gg_rename_assembly_iso.py \
            -i {config[work_dir]}/{config[prefix]}.p_ctg.fa \
            -o {output.fasta} -s {config[prefix]}
        """

rule map_reads_to_assembly:
    input:
        fa = f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        bam = config["hifi_bam"],
        # checkm_finish = f"{config['work_dir']}/{config['prefix']}.checkm.finish"
    output:
        align_bam = f"{config['work_dir']}/{config['prefix']}.align.bam",
        bai=f"{config['work_dir']}/{config['prefix']}.align.bam.bai",
        read_count = f"{config['work_dir']}/{config['prefix']}.align.count.csv",
    threads: config["threads"]
    shell:
        """
        samtools faidx {input.fa}
        ~/smrtlink/pbmm2 align --preset CCS -j {threads} {input.fa} {input.bam} {config[work_dir]}/{config[prefix]}.raw.bam
        samtools sort -T {config[work_dir]}/{config[prefix]} -@ {threads} -o {output.align_bam} {config[work_dir]}/{config[prefix]}.raw.bam
        rm {config[work_dir]}/{config[prefix]}.raw.bam

        samtools index {output.align_bam}
        /home/shuaiw/smrtlink/pbindex {output.align_bam}
        python /home/shuaiw/Methy/assembly_pipe/count_assembly.py {config[prefix]} {config[work_dir]} {input.bam}
        """

rule checkM:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        checkm=f"{config['work_dir']}/checkM2/quality_report.tsv",
        finish=f"{config['work_dir']}/{config['prefix']}.checkm.finish"
    threads: config["threads"]
    shell:
        """ 
        mkdir {config[work_dir]}/bins/
        cp {input.fasta} {config[work_dir]}/bins/
        checkm2 predict --input {config[work_dir]}/bins/ --output-directory  {config[work_dir]}/checkM2 --force -x .fa --threads {threads}
        touch {output.finish}
        """

rule GTDB:
    input:
        checkm=f"{config['work_dir']}/checkM2/quality_report.tsv",
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
          -x fa
        touch {output.gtdb_finish}
        """

rule genomad:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
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

rule call_methylation:
    input:
        genomad_finish=f"{config['work_dir']}/genomad.done",
        mge_finish = f"{config['work_dir']}/ctg_mge.done",
        bam=f"{config['work_dir']}/{config['prefix']}.align.bam",
        fa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        mge_file = f"{config['work_dir']}/all_mge.tsv",
    output:
        time=f"{config['work_dir']}/{config['prefix']}.methyl.time",
        summary=f"{config['work_dir']}/{config['prefix']}_methylation/summary.csv",
        methy_finish = f"{config['work_dir']}/methylation.finish"
    threads: config["threads"]
    # conda: "methy3"
    shell:
        """
        which python
        /usr/bin/time -v -o {output.time} /home/shuaiw/miniconda3/envs/methy3/bin/python /home/shuaiw/Methy/main.py \
          --work_dir {config[work_dir]}/{config[prefix]}_methylation \
          --whole_bam {input.bam} \
          --whole_ref {input.fa} \
          --read_type hifi \
          --min_len 1000 \
          --min_cov 3 \
          --min_iden 0.97 \
          --min_frac 0.4 \
          --min_score 30 \
          --min_sites 100 \
          --annotate_rm \
          --mge_file {input.mge_file} \
          --threads {threads} \
          --kmer_mean_db /home/shuaiw/borg/paper/run2/96plex/96plex_methylation/control/control_db.up7.down3.mean.dat \
          --kmer_num_db /home/shuaiw/borg/paper/run2/96plex/96plex_methylation/control/control_db.up7.down3.num.dat
        touch {output.methy_finish}
        """

