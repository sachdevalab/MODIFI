import os
import re

def read_list(bam_list, cmd_file, prefix_table):
    """
    Read a list of BAM files from a given file.
    """
    prefix_dict = read_prefix_table(prefix_table)
    w = open(cmd_file, 'w')
    m = open("run_methy2.sh", 'w')
    anno = open("run_anno.sh", 'w')
    borg = open("run_borg.sh", 'w')
    orphan = open("run_orphan.sh", 'w')
    drep = open("run_drep.sh", 'w')
    spacer = open("run_spacer.sh", 'w')
    # h = open(prefix_table, 'w')
    i = 1
    with open(prefix_table, 'r') as f:
        for line in f:
            items = line.strip().split()
            raw_prefix = items[0]
            prefix = items[1]
            hifi_bam = items[2]
            # print (items)
            environment = items[3]
            work_dir = os.path.join("/home/shuaiw/borg/paper/run2", prefix)

            prokka = f"{work_dir}/prokka/{prefix}.gff"

            cmd = f"""
            #### number {i}
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \\
                hifi_bam={hifi_bam} \\
                prefix={prefix} \\
                work_dir={work_dir} -j 64" \\
                --job-name={prefix}
            """
            print (cmd.strip())
            print (cmd, file=w)

            cmd = f"""
            #### number {i}
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \\
                hifi_bam={hifi_bam} \\
                prefix={prefix} \\
                work_dir={work_dir} -j 64 --use-conda" \\
                --job-name={prefix} 
            """
            if not os.path.exists(prokka):
                print (cmd, file=anno)

            spacer_cmd = f"""
                python /home/shuaiw/Methy/benchmark/spacer/spacer_match.py \\
                {prefix} \\
                 /home/shuaiw/borg/paper/run2/{prefix}/
            """
            # spacer_cmd = f"""
            # python merge_MGEs.py \\
            # /home/shuaiw/borg/paper/run2/{prefix}/ \\
            # {prefix} 
            # """
            print (spacer_cmd, file=spacer)

            borg_cmd = f"""
            python find_borg.py  /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
                /home/shuaiw/borg/paper/run2/{prefix}/borg/ 
            """
            print (borg_cmd, file=borg)

            methy2_dir = os.path.join(work_dir, f"{prefix}_methylation2")
            # if os.path.exists(methy2_dir):
            #     continue
            # methy_cmd = f"""
            #     sbatch  --partition standard --wrap "python /home/shuaiw/Methy/main.py \\
            #     --work_dir /home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation_kmer3 \\
            #     --whole_bam /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.bam \\
            #     --whole_ref /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
            #     --read_type hifi \\
            #     --min_len 2000 \\
            #     --min_cov 3 \\
            #     --min_iden 0.97 \\
            #     --min_frac 0.3 \\
            #     --min_score 30 \\
            #     --min_sites 30 \\
            #     --mge_file /home/shuaiw/borg/paper/run2/{prefix}/all_mge.tsv \\
            #     --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \\
            #     --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \\
            #     --threads 64 --clean" \\
            #     --job-name=clip_{prefix}
            # """
            methy_cmd = f"""
                sbatch  --partition standard --wrap "python /home/shuaiw/Methy/main.py \\
                --work_dir /home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3 \\
                --whole_bam /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.bam \\
                --whole_ref /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
                --read_type hifi \\
                --min_len 2000 \\
                --min_cov 3 \\
                --min_iden 0.97 \\
                --min_frac 0.3 \\
                --min_score 30 \\
                --min_sites 100 \\
                --mge_file /home/shuaiw/borg/paper/run2/{prefix}/all_mge.tsv \\
                --run_steps host \\
                --threads 64" \\
                --job-name={prefix}
            """
            print (methy_cmd.strip(), file=m)

            orphan_cmd = f"""
            /home/shuaiw/miniconda3/envs/methy3/bin/python \\
                /home/shuaiw/Methy/benchmark/orphan/motif_enrichment.py \\
                    /home/shuaiw/borg/paper/run2/{prefix}/prokka/{prefix}.gff \\
                    /home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation2
            """
            if not os.path.exists(prokka):
                print (orphan_cmd.strip(), file=orphan)

            cmd = f"""
            #### number {i}
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \\
                hifi_bam={hifi_bam} \\
                prefix={prefix} \\
                work_dir={work_dir} -j 64" \\
                --job-name={prefix}
            """
            print (cmd, file=drep)

            # print (f"{prefix}\t{prefix}\t{hifi_bam}", file=h)
            i += 1
    w.close()
    m.close()
    borg.close()
    orphan.close()
    drep.close()
    spacer.close()
    # h.close()

def read_prefix_table(prefix_table):
    """
    Read a prefix table file.
    """
    prefix_dict = {}
    with open(prefix_table, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            items = line.strip().split()
            raw_prefix = items[0]
            prefix = items[1]
            hifi_bam = items[2]
            prefix_dict[raw_prefix] = prefix
    return prefix_dict

def batch_asthma(cmd_file, prefix_table, outdir):
    """
    Read a list of BAM files from a given file.
    """
    prefix_dict = read_prefix_table(prefix_table)
    w = open(cmd_file, 'w')
    borg = open("run_borg.sh", 'w')
    i = 1
    with open(prefix_table, 'r') as f:
        for line in f:
            items = line.strip().split()
            raw_prefix = items[0]
            prefix = items[1]
            hifi_bam = items[2]
            # print (items)
            environment = items[3]

            work_dir = os.path.join(outdir, prefix)


            cmd = f"""
            #### number {i}
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \\
                hifi_bam={hifi_bam} \\
                prefix={prefix} \\
                work_dir={work_dir} -j 64" \\
                --job-name={prefix}
            """
            # cmd = f"""
            # # #### number {i}
            # snakemake -s annotation.smk --config \\
            #     hifi_bam={hifi_bam} \\
            #     prefix={prefix} \\
            #     work_dir={work_dir} -j 64 
            # """
            # cmd = f"""
            # # #### number {i}
            # sbatch --partition standard --wrap "snakemake -s annotation.smk --config \\
            #     hifi_bam={hifi_bam} \\
            #     prefix={prefix} \\
            #     work_dir={work_dir} -j 64" \\
            #     --job-name={prefix}
            # """

            cmd = f"""
                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \\
                --work_dir /home/shuaiw/borg/paper/borg_data/methy/{prefix}/{prefix}_methylation3 \\
                --whole_bam {outdir}/{prefix}/{prefix}.align.bam \\
                --whole_ref {outdir}/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
                --read_type hifi \\
                --min_len 2000 \\
                --min_cov 3 \\
                --min_iden 0.97 \\
                --min_frac 0.3 \\
                --min_score 30 \\
                --min_sites 100 \\
                --run_steps anno \\
                --mge_file {outdir}/{prefix}/all_mge.tsv \\
                --threads 64" \\
                --job-name=sue_{i}
            """

            borg_cmd_for = f"""
                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \\
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/{prefix}/{prefix}_methylation3 \\
                --unaligned_bam {hifi_bam} \\
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \\
                --read_type hifi \\
                --min_len 1000 \\
                --min_cov 2 \\
                --min_ctg_cov 2 \\
                --min_iden 0.95 \\
                --min_frac 0.1 \\
                --min_score 30 \\
                --min_sites 100 \\
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \\
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \\
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \\
                --threads 64 --visu_ipd --detect_misassembly" \\
                --job-name=borg_{i}
            """

            borg_cmd_rev = f"""
                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \\
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_rev2/{prefix}/{prefix}_methylation3 \\
                --unaligned_bam {hifi_bam} \\
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.revcomp.fa \\
                --read_type hifi \\
                --min_len 1000 \\
                --min_cov 2 \\
                --min_ctg_cov 2 \\
                --min_iden 0.95 \\
                --min_frac 0.1 \\
                --min_score 30 \\
                --min_sites 100 \\
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \\
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \\
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \\
                --threads 64 --visu_ipd --detect_misassembly" \\
                --job-name=borg_{i}
            """

            # borg_cmd = f"""
            #     sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \\
            #     --work_dir /home/shuaiw/borg/paper/borg_data/methy3/{prefix}/{prefix}_methylation3 \\
            #     --unaligned_bam {hifi_bam} \\
            #     --whole_ref /home/shuaiw/borg/paper/borg_data/jill_borgs/borg_set2.fa \\
            #     --read_type hifi \\
            #     --min_len 1000 \\
            #     --min_cov 3 \\
            #     --min_ctg_cov 3 \\
            #     --min_iden 0.90 \\
            #     --min_frac 0.2 \\
            #     --min_score 30 \\
            #     --min_sites 100 \\
            #     --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \\
            #     --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \\
            #     --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \\
            #     --threads 64" \\
            #     --job-name=borg_{i}
            # """

            # borg_cmd = f"""
            # nohup python /home/shuaiw/mGlu/assembly_pipe/../benchmark/borg/find_borg.py  /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
            #     /home/shuaiw/borg/paper/run2/{prefix}/borg/ --prefix {prefix} &
            # """

            if environment == "soil":
                print (borg_cmd_for, file=borg)
                print (borg_cmd_rev, file=borg)


            # if i   in [9]:

            # print (cmd.strip())
            # print (cmd, file=w)

            i += 1
    w.close()
    borg.close()


def batch_soil():
    """
    Read a list of BAM files from a given file.
    """
    outdir = "/home/shuaiw/borg/paper/gg_run/"
    prefix_table = "prefix_table_soil.tab"
    w = open("run_soil.sh", 'w')
    i = 1
    with open(prefix_table, 'r') as f:
        for line in f:
            items = line.strip().split()
            ref = items[0]
            prefix = items[1]
            hifi_bam = items[2]
            # print (items)
            environment = items[3]

            work_dir = os.path.join(outdir, prefix)
            cmd = f"""
            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_opt.smk --config \\
                hifi_bam={hifi_bam} \\
                prefix={prefix} \\
                work_dir={work_dir} \\
                ref={ref} -j 64 "  --job-name={prefix}
            """
            print(cmd, file=w)
    w.close()


if __name__ == "__main__":
    # bam_list = "/home/shuaiw/borg/paper/aws/bam.list"
    # cmd_file = "run.sh"
    # prefix_table = "prefix_table.tab"
    # read_list(bam_list, cmd_file, prefix_table)

    outdir = "/home/shuaiw/borg/paper/run2/"
    cmd_file = "run_asthma.sh"
    prefix_table = "prefix_table.tab"
    batch_asthma(cmd_file, prefix_table, outdir)

    # batch_soil()
