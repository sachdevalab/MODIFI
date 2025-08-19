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
    # h = open(prefix_table, 'w')
    i = 1
    with open(bam_list, 'r') as f:
        for line in f:
            hifi_bam = line.strip()
            field = hifi_bam.split("/")[-1].split(".")
            prefix = field[0]
            if re.search("soil", hifi_bam) and prefix != "XRSBK_20221007_S64018_PL100268287-1_C01":
                prefix = f"{field[0]}_{field[-2]}"
            ## replace prefix with the one in the prefix table
            if prefix in prefix_dict:
                prefix = prefix_dict[prefix]
            else:
                print(f"Warning: {prefix} not found in prefix table, using original prefix.")
            work_dir = os.path.join("/home/shuaiw/borg/paper/run2", prefix)
            cmd = f"""
            #### number {i}
            sbatch --partition standard --wrap "snakemake --config \\
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
            print (cmd, file=anno)

            borg_cmd = f"""
            python find_borg.py  /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
                /home/shuaiw/borg/paper/run2/{prefix}/borg/ 
            """
            print (borg_cmd, file=borg)

            methy2_dir = os.path.join(work_dir, f"{prefix}_methylation2")
            # if os.path.exists(methy2_dir):
            #     continue
            methy_cmd = f"""
                sbatch  --partition standard --wrap "python /home/shuaiw/Methy/main.py \\
                --work_dir /home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation2 \\
                --whole_bam /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.bam \\
                --whole_ref /home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa \\
                --read_type hifi \\
                --min_len 1000 \\
                --max_NM 30000000 \\
                --min_cov 1 \\
                --min_frac 0.4 \\
                --min_score 30 \\
                --min_sites 30 \\
                --run_steps host \\
                --mge_file /home/shuaiw/borg/paper/run2/{prefix}/all_mge.tsv \\
                --threads 64" \\
                --job-name={prefix}_test
            """
            print (methy_cmd.strip(), file=m)
            

            # print (f"{prefix}\t{prefix}\t{hifi_bam}", file=h)
            i += 1
    w.close()
    m.close()
    borg.close()
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
            items = line.strip().split("\t")
            raw_prefix = items[0]
            prefix = items[1]
            hifi_bam = items[2]
            prefix_dict[raw_prefix] = prefix
    return prefix_dict


if __name__ == "__main__":
    bam_list = "/home/shuaiw/borg/paper/aws/bam.list"
    cmd_file = "run.sh"
    prefix_table = "prefix_table.tab"
    read_list(bam_list, cmd_file, prefix_table)