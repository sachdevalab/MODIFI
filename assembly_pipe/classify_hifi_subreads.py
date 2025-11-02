import pysam
import subprocess



def classify_bam(bam):
    ## read the first read of the bam using biopython
    ## see if the read has the tag
    
    # Use grep to search for the tag in the BAM file (as text)
    # Note: BAM is a binary format, so we need to convert it to SAM first
    cmd = f"samtools view {bam} | head -n 10"
    result = subprocess.check_output(cmd, shell=True, text=True)
    if "\tfi:" in result and "\tri:" in result:
        return "hifi"
    elif "\tip:" in result and "\tpw:" in result:
        return "subreads"
    else:
        return "unknown"
    

# bam = "/home/shuaiw/borg/paper/aws/isolate/bacteria/set1/ERR10087071/demultiplex.bc1053T__bc1053T.bam"
folder = "/home/shuaiw/borg/paper/aws/isolate/bacteria/set1/"
ccs_bam_dir = "/home/shuaiw/borg/paper/aws/isolate/bacteria/ccs_bams/"

folder = "/groups/banfield/projects/multienv/methylation_temp/aws_methylation/"
ccs_bam_dir = "/groups/banfield/projects/multienv/methylation_temp/batch2_ccs_bam/"


## enumerate subfolders in folder
import os
cmd_list = []
for subfolder in os.listdir(folder):
    ## skip if it is not a directory
    if not os.path.isdir(os.path.join(folder, subfolder)):
        continue
    sra = subfolder.split(".")[0]
    subfolder_path = os.path.join(folder, subfolder)

    ## list all bam files in subfolder_path
    bam_num = 0
    for file in os.listdir(subfolder_path):
        if file.endswith(".bam") and not file.endswith(".ccs.bam"):
            bam = os.path.join(subfolder_path, file)
            bam_type = classify_bam(bam)
            print (f"{sra}\t{file}\t{bam_type}")
            bam_num += 1
    if bam_num > 1:
        print (f"############### {sra}\t has over 1 bams")
        continue
    if bam_num == 0:
        print (f"################{sra}\t has no bams")
        continue
    # print (f"Total BAM files in {subfolder}: {bam_num}")

    ccs_bam = os.path.join(ccs_bam_dir, f"{sra}.ccs.bam")
    if bam_type == "unknown":
        print (f"############### {sra}\t has unknown bam type")
        continue
    elif bam_type == "subreads":
        cmd = f'sbatch  --partition standard --job-name={sra} --wrap "/home/shuaiw/smrtlink/ccs {bam} {ccs_bam} --hifi-kinetics --min-rq 0.99 --min-passes 3  --num-threads 64" '
        cmd = f'/home/shuaiw/smrtlink/ccs {bam} {ccs_bam} --hifi-kinetics --min-rq 0.99 --min-passes 3  --num-threads 64'
        cmd_list.append(cmd)
    else:
        print (f"############### {sra}\t is already hifi")
        ## link the bam to ccs_bam
        cmd = f'ln -s {bam} {ccs_bam}'
        os.system(cmd)
    # cmd =  f"""sbatch  --partition standard --wrap "snakemake -s isolation.smk \\
    #             --config hifi_bam={ccs_bam} \\
    #             prefix={sra} \\
    #             work_dir=/home/shuaiw/borg/paper/isolation/bacteria/{sra} \\
    #             -j 64"  --job-name={sra[3:]}"""
    # cmd = f"""
    #     sbatch  --partition standard --wrap "python /home/shuaiw/Methy/main.py \\
    #     --work_dir /home/shuaiw/borg/paper/isolation/bacteria/{sra}/{sra}_methylation \\
    #     --whole_bam {ccs_bam} \\
    #     --whole_ref /home/shuaiw/borg/paper/isolation/bacteria/{sra}/{sra}.hifiasm.p_ctg.rename.fa \\
    #     --read_type hifi \\
    #     --min_len 1000 \\
    #     --min_cov 10 \\
    #     --min_frac 0.1 \\
    #     --min_score 30 \\
    #     --min_sites 10 \\
    #     --mge_file /home/shuaiw/borg/paper/isolation/bacteria/{sra}/all_mge.tsv \\
    #     --threads 30 --run_steps motif merge profile host" \\
    #     --job-name={sra[3:]}
    # """
    # cmd =  f"""sbatch  --partition standard --wrap "snakemake -s gtdb_isolation.smk \\
    #             --config hifi_bam={ccs_bam} \\
    #             prefix={sra} \\
    #             work_dir=/home/shuaiw/borg/paper/isolation/bacteria/{sra} \\
    #             -j 64"  --job-name={sra[3:]}"""
    


convert = "get_hifi.sh"
f = open(convert, "w")
## separately run the commands using ten scripts
num_scripts = 10
for i in range(num_scripts):
    script_file = f"batch/run_ccs_part_{i}.sh"
    with open(script_file, "w") as sf:
        for j in range(i, len(cmd_list), num_scripts):
            sf.write(cmd_list[j] + "\n")
    print (f"sbatch  --partition standard --wrap 'bash {script_file}' --job-name={i}", file = f)
f.close()
