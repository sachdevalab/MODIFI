# Snakefile

# Configuration
configfile: "config.yaml"

rule all_mge:
    input:
        virsorter2=f"{config['work_dir']}/virsorter2.done",
        vibrantr=f"{config['work_dir']}/vibrant.done",

        anvi_done=f"{config['work_dir']}/anvi.done",
        dram_finish = f"{config['work_dir']}/dram2/dram.finish",
        prokka_finish = f"{config['work_dir']}/prokka/prokka.finish",
        prodigal_gv_finish=f"{config['work_dir']}/prodigal/{config['prefix']}.prodigal.finish",
        plasX_done=f"{config['work_dir']}/plasX.done",

        # enrichment_finish = f"{config['work_dir']}/motif_enrichment.done",



# plasX to call plasmids
# rule anvi:
#     input:
#         fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
#     output:
#         anvi_done=f"{config['work_dir']}/anvi.done"
#     params:
#         prefix=f"{config['work_dir']}/{config['prefix']}",
#     threads: config["threads"]
#     shell:
#         """
#         anvi-gen-contigs-database -L 0 -T {threads} --project-name {params.prefix} -f {input.fasta} -o {params.prefix}.db

#         anvi-export-gene-calls --gene-caller prodigal -c {params.prefix}.db -o {params.prefix}-gene-calls.txt

#         anvi-run-ncbi-cogs -T {threads} --cog-version COG14 --cog-data-dir /home/shuaiw/borg/paper/anvio_db/COG_2014 -c {params.prefix}.db

#         anvi-run-pfams -T {threads} --pfam-data-dir /home/shuaiw/borg/paper/anvio_db/Pfam_v32 -c {params.prefix}.db

#         anvi-export-functions --annotation-sources COG14_FUNCTION,Pfam -c {params.prefix}.db -o {params.prefix}-cogs-and-pfams.txt

#         touch {output.anvi_done}
#         """

rule plasX:
    input:
        anvi_done=f"{config['work_dir']}/anvi.done",
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        plasX_done=f"{config['work_dir']}/plasX.done"
    params:
        prefix=f"{config['work_dir']}/{config['prefix']}",
    threads: config["threads"]
    conda:"/home/sdiamond/miniconda3/envs/plasx"
    shell:
        """
        /home/sdiamond/miniconda3/envs/plasx/bin/plasx search_de_novo_families \
            -g {params.prefix}-gene-calls.txt \
            -o {params.prefix}-de-novo-families.txt \
            --threads {threads} \
            --overwrite

        /home/sdiamond/miniconda3/envs/plasx/bin/plasx predict -a {params.prefix}-cogs-and-pfams.txt {params.prefix}-de-novo-families.txt \
            -g {params.prefix}-gene-calls.txt \
            -o {params.prefix}-plasmid_nr-scores.txt \
            --overwrite
        touch {output.plasX_done}
        """



rule virsorter2:
    input:
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
        plasX_done=f"{config['work_dir']}/plasX.done",

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


## use pyrodigal-gv set threads
rule prodigal_gv:
    input:
        anvi_done=f"{config['work_dir']}/anvi.done",
        vibrantr=f"{config['work_dir']}/vibrant.done",
        fasta=f"{config['work_dir']}/{config['prefix']}.hifiasm.p_ctg.rename.fa",
    output:
        faa=f"{config['work_dir']}/prodigal/{config['prefix']}.faa",
        prodigal_gv_finish=f"{config['work_dir']}/prodigal/{config['prefix']}.prodigal.finish"
    threads: config["threads"]
    shell:
        """
        pyrodigal-gv -p meta -m -i {input.fasta} -a {output.faa} -j {threads}
        touch {output.prodigal_gv_finish}
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