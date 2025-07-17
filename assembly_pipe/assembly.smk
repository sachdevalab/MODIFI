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

rule remove_host:
    input:
        fastq=f"{config['work_dir']}/{config['prefix']}.hifi.qc.fq.gz"
    output:
        bam=temp(f"{config['work_dir']}/{config['prefix']}.minimap2.human.sort.bam"),
        human_reads=f"{config['work_dir']}/{config['prefix']}-human-mapped-reads.fastq",
        # stats=f"{config['work_dir']}/{config['prefix']}-human-stats.tsv",
        filtered_fastq=f"{config['work_dir']}/{config['prefix']}.hifi.qc.HR.fq.gz"
    params:
        ref=host_ref
    threads: config["threads"]
    shell:
        """
        minimap2 -t {threads} -ax map-hifi {params.ref} {input.fastq} | \
            samtools view -b -T {params.ref} -F 4 | \
            sambam > {output.bam}
        samtools view {output.bam} | cut -f 1 | sort | uniq > {output.human_reads}
        filterbyname.sh in={input.fastq} out={output.filtered_fastq} names={output.human_reads} include=f
        """

rule assemble_hifiasm:
    input:
        fq=f"{config['work_dir']}/{config['prefix']}.hifi.qc.HR.fq.gz"
    output:
        gfa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.gfa",
    log: f"{config['work_dir']}/{config['prefix']}.hifiasm.log"
    threads: config["threads"]
    shell:
        """
        hifiasm_meta --no-binning -o {config[work_dir]}/{config[prefix]}.hifiasm \
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
        /home/rohan/dev/pipeline/workflow/scripts/gg_rename_assembly.py \
            -i {config[work_dir]}/{config[prefix]}.p_ctg.fa \
            -o {output.fasta} -s {config[prefix]}
        """

rule map_reads_to_assembly:
    input:
        fa = f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        bam = config["hifi_bam"],
        # checkm_finish = f"{config['work_dir']}/{config['prefix']}.checkm.finish"
    output:
        align_bam = temp(f"{config['work_dir']}/{config['prefix']}.align.bam"),
        bai=f"{config['work_dir']}/{config['prefix']}.align.bam.bai"
    threads: config["threads"]
    shell:
        """
        samtools faidx {input.fa}
        ~/smrtlink/pbmm2 align --preset CCS -j {threads} {input.fa} {input.bam} {config[work_dir]}/{config[prefix]}.raw.bam
        samtools sort -T {config[work_dir]}/{config[prefix]} -@ {threads} -o {output.align_bam} {config[work_dir]}/{config[prefix]}.raw.bam
        rm {config[work_dir]}/{config[prefix]}.raw.bam

        samtools index {output.align_bam}
        /home/shuaiw/smrtlink/pbindex {output.align_bam}
        """

rule call_methylation:
    input:
        bam=f"{config['work_dir']}/{config['prefix']}.align.bam",
        fa=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        time=f"{config['work_dir']}/{config['prefix']}.methyl.time",
        summary=f"{config['work_dir']}/{config['prefix']}_methylation/summary.csv",
        methy_finish = f"{config['work_dir']}/methylation.finish"
    threads: config["threads"]
    shell:
        """
        /usr/bin/time -v -o {output.time} python /home/shuaiw/Methy/main.py \
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
          --threads {threads} 
        touch {output.methy_finish}

        """