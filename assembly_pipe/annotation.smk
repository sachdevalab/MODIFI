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
        touch {output.finish}
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
            -j {threads}
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
# rule plasx:
#     input:
#         fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
#     output:
#         plasx_done=f"{config['work_dir']}/plasx.done"
#     threads: config["threads"]
#     shell:
#         """
#         plasx predict -i {input.fasta} -o {config[work_dir]}/plasx/ -t {threads}
#         touch {output.plasx_done}
#         """
## use pyrodigal-gv set threads
rule prodigal_gv:
    input:
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
        faa = f"{config['work_dir']}/tr_spacer_hit_filtered/prokka/tr.faa",
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
        finish = f"{config['work_dir']}/ctg_mge.done"
    params:
        output_dir = config['work_dir'],
        prefix = config["prefix"]
    shell:
        """
        python merge_MGEs.py {params.output_dir} {params.prefix}
        touch {output.finish}
        """






## genomad PlasFlow  VirSorter / VirSorter2 VIBRANT

