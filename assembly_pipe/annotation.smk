# Snakefile

# Configuration
configfile: "config.yaml"

rule all_annotation:
    input:
        checkm=f"{config['work_dir']}/checkM2/quality_report.tsv",
        gtdb_finish=f"{config['work_dir']}/GTDB/gtdbtk.done",
        genomad_finish=f"{config['work_dir']}/genomad.done",
        virsorter2=f"{config['work_dir']}/virsorter2.done",
        vibrantr=f"{config['work_dir']}/vibrant.done",
        faa = f"{config['work_dir']}/prodigal/{config['prefix']}.faa",
        dram_finish = f"{config['work_dir']}/dram2/dram.finish",
        prokka_finish = f"{config['work_dir']}/prokka/prokka.finish",
        ctg_mge = f"{config['work_dir']}/ctg_mge.done",
        finish=f"{config['work_dir']}/prodigal/{config['prefix']}.prodigal.finish",
        checkv_finish=f"{config['work_dir']}/checkV.done",
        anvi_done=f"{config['work_dir']}/anvi.done",
        drep_clu_file = f"{config['work_dir']}/dRep_out/data_tables/Cdb.csv",
        drep_finish = f"{config['work_dir']}/dRep.finish",

rule checkM:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
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
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
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

rule virsorter2:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        genomad_finish=f"{config['work_dir']}/genomad.done"

    output:
        virsorter2=f"{config['work_dir']}/virsorter2.done"
    threads: config["threads"]
    shell:
        """
        virsorter run all \
            --seqfile {input.fasta} \
            -w  {config[work_dir]}/virsorter2/ \
            -j {threads} \
            --tmpdir {config[work_dir]}/tmp/ 
        touch {output.virsorter2}
        """

rule vibrant:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        virsorter2=f"{config['work_dir']}/virsorter2.done"
    output:
        vibrantr=f"{config['work_dir']}/vibrant.done"
    threads: config["threads"]
    shell:
        """
        VIBRANT_run.py -i {input.fasta} -folder {config[work_dir]}/vibrant/ -t {threads} 
        touch {output.vibrantr}
        """

## plasX to call plasmids
rule anvi:
    input:
        vibrantr=f"{config['work_dir']}/vibrant.done",
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        anvi_done=f"{config['work_dir']}/anvi.done"
    params:
        prefix=f"{config['work_dir']}/{config['prefix']}",
    threads: config["threads"]
    shell:
        """
        anvi-gen-contigs-database -L 0 -T {threads} --project-name {params.prefix} -f {input.fasta} -o {params.prefix}.db

        anvi-export-gene-calls --gene-caller prodigal -c {params.prefix}.db -o {params.prefix}-gene-calls.txt

        anvi-run-ncbi-cogs -T {threads} --cog-version COG14 --cog-data-dir /home/shuaiw/borg/paper/anvio_db/COG_2014 -c {params.prefix}.db

        anvi-run-pfams -T {threads} --pfam-data-dir /home/shuaiw/borg/paper/anvio_db/Pfam_v32 -c {params.prefix}.db

        anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c {params.prefix}.db -o {params.prefix}-cogs-and-pfams.txt

        touch {output.anvi_done}
        """
## use pyrodigal-gv set threads
rule prodigal_gv:
    input:
        anvi_done=f"{config['work_dir']}/anvi.done",
        vibrantr=f"{config['work_dir']}/vibrant.done",
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        faa=f"{config['work_dir']}/prodigal/{config['prefix']}.faa",
        finish=f"{config['work_dir']}/prodigal/{config['prefix']}.prodigal.finish"
    threads: config["threads"]
    shell:
        """
        pyrodigal-gv -p meta -m -i {input.fasta} -a {output.faa} -j {threads}
        touch {output.finish}
        """

rule dram:
    input:
        faa = f"{config['work_dir']}/prodigal/{config['prefix']}.faa",

    output:
        f"{config['work_dir']}/dram2/annotations.tsv",
        dram_finish = f"{config['work_dir']}/dram2/dram.finish"

    params:
        output_dir =  f"{config['work_dir']}/dram2",

    threads:
        config["threads"]

    shell:
        """
        rm -r {params.output_dir}
        DRAM.py annotate_genes \
        -i {input.faa} \
        -o {params.output_dir} \
        --threads {threads} \

        touch {output.dram_finish}
        """

#### how to do rRNA and tRNA from rohan
rule prokka:
    input:
        faa = f"{config['work_dir']}/prodigal/{config['prefix']}.faa",
        fasta = f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        dram_finish = f"{config['work_dir']}/dram2/dram.finish",
    output:
        faa = f"{config['work_dir']}/prokka/{config['prefix']}.faa",
        tsv = f"{config['work_dir']}/prokka/{config['prefix']}.tsv",
        prokka_finish = f"{config['work_dir']}/prokka/prokka.finish"

    params:
        output_dir = f"{config['work_dir']}/prokka",
        prefix = config["prefix"],

    threads:
        config["threads"]

    shell:
        """
        prokka --outdir {params.output_dir} --prefix {params.prefix} \
            --cpus {threads} --force \
            --proteins {input.faa} {input.fasta}
        touch {output.prokka_finish}
        """

rule get_ctg_mge:
    input:
        prokka_finish = f"{config['work_dir']}/prokka/prokka.finish"
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


rule drep:
    input:
        mge_finish = f"{config['work_dir']}/ctg_mge.done",
    output:
        drep_clu_file = f"{config['work_dir']}/dRep_out/data_tables/Cdb.csv",
        drep_finish = f"{config['work_dir']}/dRep.finish",
    params:
        genome_list = f"{config['work_dir']}/genome.list",
    threads: config["threads"]
    shell:
        """
        ls {config[work_dir]}/bins/*.fasta > {params.genome_list}
        dRep dereplicate \
        -p {threads} \
        -g {params.genome_list} \
        -comp 50 \
        -con 10 \
        --S_algorithm skani \
        -ms 10000 \
        -sa 0.95 \
        -nc 0.7 {config[work_dir]}/dRep_out
        touch {output.drep_finish}
        """


### motif enrichment analysis
# /home/shuaiw/Methy/benchmark/orphan/motif_enrichment.py
