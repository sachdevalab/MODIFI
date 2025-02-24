# Snakefile

import glob
import os

configfile: "config.yaml"

whole_bam = config["whole_bam"]
whole_ref = config["whole_ref"]
work_dir = config["work_dir"]

print (whole_bam)
print (whole_ref)
print (work_dir)

## create work_dir
os.makedirs(work_dir, exist_ok=True)

rule all:
    input:
        os.path.join(work_dir, "split_done.txt"),
        os.path.join(work_dir, "load_done.txt"),
        os.path.join(work_dir, "control_done.txt"),
        os.path.join(work_dir, "compare_done.txt"),
        os.path.join(work_dir, "motif_done.txt")

rule split:
    input:
    output:
        # touch(os.path.join(work_dir, "split_done.txt"))
        os.path.join(work_dir, "split_done.txt")
    shell:
        """
        cd split 
        /usr/bin/time -v -o {output} snakemake --config whole_bam={whole_bam} whole_ref={whole_ref} work_dir={work_dir}
        cd ..
        """

rule load:
    input:
        os.path.join(work_dir, "split_done.txt")
    output:
        os.path.join(work_dir, "load_done.txt")
    params:
        maxAlignments=config.get("maxAlignments", 10000)
    shell:
        """
        cd load 
        /usr/bin/time -v -o {output} snakemake --config whole_bam={whole_bam} whole_ref={whole_ref} work_dir={work_dir} maxAlignments={params.maxAlignments}
        cd ..
        """

rule control:
    input:
        os.path.join(work_dir, "load_done.txt")
    output:
        os.path.join(work_dir, "control_done.txt")
    shell:
        """
        cd ipdtools 
        /usr/bin/time -v -o {output} snakemake --config whole_bam={whole_bam} whole_ref={whole_ref} work_dir={work_dir}
        cd ..
        """ 

rule compare:
    input:
        os.path.join(work_dir, "control_done.txt")
    output:
        os.path.join(work_dir, "compare_done.txt")
    shell:
        """
        cd compare 
        /usr/bin/time -v -o {output} snakemake --config whole_bam={whole_bam} whole_ref={whole_ref} work_dir={work_dir}
        cd ..
        """

rule motif:
    input:
        os.path.join(work_dir, "compare_done.txt")
    output:
        os.path.join(work_dir, "motif_done.txt")
    shell:
        """
        cd motif 
        /usr/bin/time -v -o {output} snakemake --config whole_bam={whole_bam} whole_ref={whole_ref} work_dir={work_dir}
        cd ..
        """
