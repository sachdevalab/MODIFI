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
    

bam = "/home/shuaiw/borg/paper/aws/isolate/bacteria/set1/ERR10087071/demultiplex.bc1053T__bc1053T.bam"
folder = "/home/shuaiw/borg/paper/aws/isolate/bacteria/set1/"
ccs_bam_dir = "/home/shuaiw/borg/paper/aws/isolate/bacteria/ccs_bams/"
convert = "get_hifi.sh"
f = open(convert, "w")
## enumerate subfolders in folder
import os
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
            print (f"{sra}\t{file}\t{classify_bam(bam)}")
            bam_num += 1
    if bam_num > 1:
        print (f"{sra}\t has over 1 bams")
    if bam_num == 0:
        print (f"{sra}\t has no bams")
    # print (f"Total BAM files in {subfolder}: {bam_num}")

    ccs_bam = os.path.join(ccs_bam_dir, f"{sra}.ccs.bam")
    # cmd = f'sbatch  --partition standard --job-name={sra} --wrap "/home/shuaiw/smrtlink/ccs {bam} {ccs_bam} --hifi-kinetics --min-rq 0.99 --min-passes 3  --num-threads 64" '
    cmd =  f"""sbatch  --partition standard --wrap "snakemake -s isolation.smk \\
                --config hifi_bam={ccs_bam} \\
                prefix={sra} \\
                work_dir=/home/shuaiw/borg/paper/isolation/bacteria/{sra} \\
                -j 64"  --job-name={sra[3:]}"""
    print (cmd, file=f)
f.close()
