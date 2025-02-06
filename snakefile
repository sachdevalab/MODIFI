# Snakefile

import glob

# Define directories
# bam_dir = "/home/shuaiw/methylation/data/borg/b_contigs"
# fa_dir = "/home/shuaiw/methylation/data/borg/b_contigs/contigs"
# output_dir = "/home/shuaiw/methylation/data/borg/b_contigs/test1"

# bam_dir = "/home/shuaiw/methylation/data/borg/large_contigs/bams"
# fa_dir = "/home/shuaiw/methylation/data/borg/large_contigs/contigs"
# output_dir = "/home/shuaiw/methylation/data/borg/large_contigs/ipd1"

whole_bam = "/home/shuaiw/borg/all_borg/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
whole_ref = "/home/shuaiw/borg/all_borg/all_borg.fasta"
work_dir = "/home/shuaiw/methylation/data/borg/new_test"

# bam_dir = "/home/shuaiw/borg/all_borg/bams"
# fa_dir = "/home/shuaiw/borg/all_borg/contigs"
# output_dir = "/home/shuaiw/borg/all_borg/ipd1"



## check if output_dir exists, if not, create it
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
if not os.path.exists(output_dir+"/ipd"):
    os.makedirs(output_dir+"/ipd")
if not os.path.exists(output_dir+"/control"):
    os.makedirs(output_dir+"/control")
if not os.path.exists(output_dir+"/ipd_ratio"):
    os.makedirs(output_dir+"/ipd_ratio")

# Get unique sample names by removing ".align.bam" or ".fa"
samples = list(set(f.split("/")[-1].replace(".bam", "").replace(".fa", "") for f in glob.glob(f"{bam_dir}/*.bam")))

# rule all:
#     input:
#         expand(f"{output_dir}/{{sample}}.ipd1.csv", sample=samples)

rule all:
    input:
        expand(f"{output_dir}/ipd_ratio/{{sample}}.ipd3.csv", sample=samples)

rule split_bam:
    input:
        whole_bam = f"{whole_bam}",
        whole_ref = f"{whole_ref}",
        work_dir = f"{work_dir}"
    output:
        f"{work_dir}/bams/{{sample}}.bam"
    shell:
        """
        set -euo pipefail
        python split_bam.py --bam {input.whole_bam} --whole_ref {input.whole_ref} --workdir {input.work_dir} --threads 100
        """

rule load_ipd:
    input:
        bam=f"{work_dir}/bams/{{sample}}.bam",
        ref=f"{work_dir}/contigs/{{sample}}.fa"
    output:
        f"{work_dir}/ipd/{{sample}}.ipd1.csv"
    shell:
        """
        set -euo pipefail
        python standard_load5.py {input.bam} {input.ref} {output} 
        """

rule control:
    input:
        ipd=f"{work_dir}/ipd/{{sample}}.ipd1.csv"
    output:
        f"{work_dir}/control/{{sample}}.ipd2.csv"
    shell:
        """
        ls {work_dir}/contigs/*.fa > {work_dir}/contigs_list.txt
        set -euo pipefail
        ./src/test {input.ipd} {output} {work_dir}/contigs_list.txt 100
        """

rule ipd_ratio:
    input:
        control=f"{work_dir}/control/{{sample}}.ipd2.csv"
    output:
        f"{work_dir}/ipd_ratio/{{sample}}.ipd3.csv"
    shell:
        """
        set -euo pipefail
        python comp_ipd_ratio.py {input.control} {output}
        """
