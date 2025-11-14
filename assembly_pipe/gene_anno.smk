import os
import glob

# Get all FASTA files from the meta_genomes directory
META_GENOMES_DIR = "/home/shuaiw/borg/paper/isolation/GTDB_tree/meta_genomes"
RESULT_DIR = "/home/shuaiw/borg/paper/gene_anno/meta"
FASTA_FILES = glob.glob(os.path.join(META_GENOMES_DIR, "*.fa"))
GENOME_NAMES = [os.path.basename(f).replace('.fa', '') for f in FASTA_FILES]

rule all:
    input:
        expand(f"{RESULT_DIR}/{{genome}}/prokka/prokka.finish", genome=GENOME_NAMES)



#### how to do rRNA and tRNA from rohan
rule prokka:
    input:
        fasta = f"{META_GENOMES_DIR}/{{genome}}.fa"
    output:
        faa = f"{RESULT_DIR}/{{genome}}/prokka/{{genome}}.faa",
        tsv = f"{RESULT_DIR}/{{genome}}/prokka/{{genome}}.tsv",
        prokka_finish = f"{RESULT_DIR}/{{genome}}/prokka/prokka.finish"

    params:
        output_dir = f"{RESULT_DIR}/{{genome}}/prokka",
        prefix = "{genome}"
    threads: 64
    shell:
        """
        mkdir -p {params.output_dir}
        prokka --outdir {params.output_dir} --prefix {params.prefix} \
            --cpus {threads} --force {input.fasta}
        touch {output.prokka_finish}
        """

