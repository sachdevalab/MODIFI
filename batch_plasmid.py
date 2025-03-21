

all_file = "/home/shuaiw/borg/allison/AG-PacBio-files-03-20-2025.txt" 

def map_cmd():
    f = open(all_file)
    f.readline()


    cmd_file = "/home/shuaiw/borg/allison/map_batch.sh"
    h = open(cmd_file, "w")

    header = """#!/bin/bash
    #SBATCH --job-name=allison      # Job name
    #SBATCH --partition=standard # Partition name

    ## construct output directory if not exists
    """
    print (header, file=h)

    for line in f:
        line = line.strip()
        field = line.split("\t")
        # print (field)

        cmd = f"""
        outdir=/home/shuaiw/borg/allison/batch
        ref={field[2]}

        if [ ! -d $outdir ]; then
            mkdir $outdir
        fi


        ccs_bam={field[1]}
        prefix=$outdir/{field[0]}
        align_bam=$prefix.align.ccs.bam


        /usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $ccs_bam $align_bam \
        --preset CCS -j $SLURM_CPUS_ON_NODE --sort -J $SLURM_CPUS_ON_NODE
        samtools index $align_bam
        /home/shuaiw//smrtlink/pbindex $align_bam
        """
        
        print (cmd, file=h)

    f.close()
    h.close()


def mod_cmd():
    f = open(all_file)
    f.readline()


    cmd_file = "/home/shuaiw/borg/allison/mt_batch.sh"
    h = open(cmd_file, "w")

    header = """#!/bin/bash
    #SBATCH --job-name=mod      # Job name
    #SBATCH --partition=standard # Partition name

    ## construct output directory if not exists
    """
    print (header, file=h)

    for line in f:
        line = line.strip()
        field = line.split("\t")
        # print (field)

        cmd = f"""
        outdir=/home/shuaiw/borg/allison/MTase
        prefix=$outdir/{field[0]}_RM

        MicrobeMod annotate_rm  -f {field[2]} -o $prefix -t 10
        """
        
        print (cmd, file=h)

    f.close()
    h.close()

def methy_cmd():
    f = open(all_file)
    f.readline()


    cmd_file = "/home/shuaiw/Methy/workflow/methy_batch.sh"
    h = open(cmd_file, "w")

    header = """
    """
    print (header, file=h)

    for line in f:
        line = line.strip()
        field = line.split("\t")
        # print (field)
        if field[0] == "NANO_3_INF1340011_8PB":
            continue

        cmd = f"""
        sbatch --partition standard --wrap "snakemake \\
    --config whole_bam=/home/shuaiw/borg/allison/batch/{field[0]}.align.ccs.bam \\
    whole_ref={field[2]} \\
    work_dir=/home/shuaiw/borg/allison/batch/{field[0]}_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\\
    --job-name={field[0]}
        """
        cmd2 = f"""
        sbatch --partition standard --wrap "snakemake \\
    --config whole_bam=/home/shuaiw/borg/allison/batch/{field[0]}.align.ccs.bam \\
    whole_ref={field[2]} \\
    work_dir=/home/shuaiw/borg/allison/batch/{field[0]} read_type=ccs min_len=5000 max_NM=30 min_cov=5"\\
    --job-name={field[0]}
        """
        
        print (cmd, file=h)

    f.close()
    h.close()

methy_cmd()


    
