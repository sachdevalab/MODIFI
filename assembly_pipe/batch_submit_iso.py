import pysam
import subprocess
import os


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
    


def convert_hifi(folder, ccs_bam_dir):
    cmd_list = []
    num = 0
    for subfolder in os.listdir(folder):
        ## skip if it is not a directory
        if not os.path.isdir(os.path.join(folder, subfolder)):
            continue
        sra = subfolder.split(".")[0]
        subfolder_path = os.path.join(folder, subfolder)
        ccs_bam = os.path.join(ccs_bam_dir, f"{sra}.ccs.bam")

        ## skip if ccs_bam already exists
        if os.path.exists(ccs_bam):
            continue
        # print (sra)
        num += 1
        ## list all bam files in subfolder_path
        bam_num = 0
        bam_dict = {}
        for file in os.listdir(subfolder_path):
            
            if file.endswith(".bam") and not file.endswith(".merge.bam"):
                bam = os.path.join(subfolder_path, file)
                ## record the bam size
                bam_size = os.path.getsize(bam)
                bam_dict[bam] = bam_size
                bam_type = classify_bam(bam)
                # print (f"{sra}\t{file}\t{bam_type}")
                bam_num += 1
        if bam_num > 1:
            print (f"############### {sra}\t has over 1 bams")
            ## select the largest bam
            # sorted_bam = sorted(bam_dict.items(), key=lambda x: x[1], reverse=True)
            # bam = sorted_bam[0][0]
            bam = os.path.join(subfolder_path, f"{sra}.merge.bam")
            # use samtools to merge the bam files in bam_list into bam
            merge_cmd = f"samtools merge {bam} " + " ".join(list(bam_dict.keys()))
            # if not os.path.exists(bam):
            #     os.system(merge_cmd)    
            print(merge_cmd)
        if bam_num == 0:
            print (f"################{sra}\t has no bams")
            continue
        # print (f"Total BAM files in {subfolder}: {bam_num}")

        
        if bam_type == "unknown":
            print (f"############### {sra}\t has unknown bam type")
            continue
        elif bam_type == "subreads":
            cmd = f'sbatch  --partition standard --job-name={sra} --wrap "/home/shuaiw/smrtlink/ccs {bam} {ccs_bam} --hifi-kinetics --min-rq 0.99 --min-passes 3  --num-threads 64" '
            cmd = f'/home/shuaiw/smrtlink/ccs {bam} {ccs_bam} --hifi-kinetics --min-rq 0.99 --min-passes 3  --num-threads 64'
            if not os.path.exists(ccs_bam):
                # os.system(cmd)
                cmd_list.append(cmd)
        else:
            print (f"############### {sra}\t is already hifi")
            ## link the bam to ccs_bam (force overwrite if exists)
            cmd = f'ln -sf {bam} {ccs_bam}'
            os.system(cmd)
        
    convert = "get_hifi.sh"
    f = open(convert, "w")
    ## separately run the commands using ten scripts
    num_scripts = 5
    for i in range(num_scripts):
        script_file = f"batch/run_ccs_part_{i}.sh"
        with open(script_file, "w") as sf:
            for j in range(i, len(cmd_list), num_scripts):
                sf.write(cmd_list[j] + "\n")
        # print (f"sbatch  --partition standard --wrap 'bash {script_file}' --job-name={i}", file = f)
        print (f"nohup bash {script_file} &", file = f)
    f.close()
    print (num)
    print (len(cmd_list))


def batch_run(ccs_bam_dir, work_dir):
    ## get all ccs bam files
    skip_sra = ["ERR13342925", "ERR13342926","ERR6536208","SRR13008124"]
    cmd_list = []
    for file in os.listdir(ccs_bam_dir):
        if file.endswith(".ccs.bam"):
            ccs_bam = os.path.join(ccs_bam_dir, file)
            sra = file.split(".")[0]
            if sra in skip_sra:
                continue
            cmd =  f"""snakemake -s isolation.smk \\
                        --config hifi_bam={ccs_bam} \\
                        prefix={sra} \\
                        work_dir={work_dir}/{sra} \\
                        -j 64  """
            finish_file = os.path.join(work_dir, sra, f"GTDB_2/gtdbtk.done")
            if not os.path.exists(finish_file):
                cmd_list.append(cmd)
    batch_file = "run_all_isolation.sh"
    num_scripts = 1
    f = open(batch_file, "w")
    for i in range(num_scripts):
        script_file = f"batch/run_isolation_part_{i}.sh"
        with open(script_file, "w") as sf:
            for j in range(i, len(cmd_list), num_scripts):
                sf.write(cmd_list[j] + "\n")
        print (f"sbatch  --partition standard --wrap 'bash {script_file}' --job-name=iso_{i}", file = f)
    print (f"Total isolation jobs: {len(cmd_list)}")
    f.close()

def methy_run(results_dir):
    ## get all ccs bam files
    cmd_list = []
    for folder in os.listdir(results_dir):
        prefix = folder
        aligned_bam = os.path.join(results_dir, folder, f"{prefix}.aligned.bam")
        finish_file = os.path.join(results_dir, folder, f"{prefix}_methylation4/methylation.finish")
        rm_file = os.path.join(results_dir, folder, f"{prefix}_methylation4//RM_systems/all_ctgs_RM.rm.genes.tsv")
        if not os.path.exists(finish_file):

            cmd =  f"""snakemake -s methy_isolation.smk --config prefix={prefix} \
                work_dir={results_dir}/{prefix} -j 64"""
            cmd_list.append(cmd)
    batch_file = "run_methy_isolation.sh"
    num_scripts = 9
    f = open(batch_file, "w")
    for i in range(num_scripts):
        script_file = f"batch/methy_isolation_part_{i}.sh"
        with open(script_file, "w") as sf:
            ## add header
            sf.write("#!/bin/bash\n")
            for j in range(i, len(cmd_list), num_scripts):
                sf.write(cmd_list[j] + "\n")
        # print (f"nohup bash {script_file} &", file = f)
        print (f"sbatch  --partition standard --wrap 'bash {script_file}'  --job-name=iso_{i} ", file = f)
    f.close()
    print (len(cmd_list))

def RM_run(results_dir):
    ## get all ccs bam files
    cmd_list = []
    for folder in os.listdir(results_dir):
        prefix = folder
        aligned_bam = os.path.join(results_dir, folder, f"{prefix}.aligned.bam")
        finish_file = os.path.join(results_dir, folder, f"{prefix}_methylation2/methylation.finish")
        rm_file = os.path.join(results_dir, folder, f"{prefix}_methylation2/RM_systems/all_ctgs_RM.rm.genes.tsv")
        if not os.path.exists(rm_file):

            cmd =  f"""snakemake -s methy_isolation.smk --config prefix={prefix} \
                work_dir={results_dir}/{prefix} -j 64"""
            cmd_list.append(cmd)
    batch_file = "run_methy_isolation.sh"
    num_scripts = 9
    f = open(batch_file, "w")
    for i in range(num_scripts):
        script_file = f"batch/methy_isolation_part_{i}.sh"
        with open(script_file, "w") as sf:
            for j in range(i, len(cmd_list), num_scripts):
                sf.write(cmd_list[j] + "\n")
        print (f"sbatch  --partition standard --wrap 'bash {script_file}'  --job-name=iso_{i}", file = f)
    f.close()
    print (len(cmd_list))

# bam = "/home/shuaiw/borg/paper/aws/isolate/bacteria/set1/ERR10087071/demultiplex.bc1053T__bc1053T.bam"
# folder = "/home/shuaiw/borg/paper/aws/isolate/bacteria/set1/"
# ccs_bam_dir = "/home/shuaiw/borg/paper/aws/isolate/bacteria/ccs_bams/"

folder = "/home/shuaiw/borg/paper/isolation/aws_methylation2/"
ccs_bam_dir = "/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/"
results_dir = "/home/shuaiw/borg/paper/isolation/batch2_results/"
# batch_run(ccs_bam_dir, results_dir)
# convert_hifi(folder, ccs_bam_dir)
methy_run(results_dir)
# RM_run(results_dir)

