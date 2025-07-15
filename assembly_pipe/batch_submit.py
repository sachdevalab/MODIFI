import os

def read_list(bam_list, cmd_file):
    """
    Read a list of BAM files from a given file.
    """
    w = open(cmd_file, 'w')
    i = 1
    with open(bam_list, 'r') as f:
        for line in f:
            hifi_bam = line.strip()
            prefix = hifi_bam.split("/")[-1].split(".")[0]
            work_dir = os.path.join("/home/shuaiw/borg/paper", prefix)
            cmd = f"""
            #### number {i}
            sbatch --partition standard --wrap "snakemake --config \\
                hifi_bam={hifi_bam} \\
                prefix={prefix} \\
                work_dir={work_dir}" \\
                --job-name={prefix}
            """
            print (cmd.strip())
            print (cmd, file=w)
            i += 1
    w.close()




if __name__ == "__main__":
    bam_list = "/home/shuaiw/borg/paper/aws/bam.list"
    cmd_file = "run.sh"
    read_list(bam_list, cmd_file)