import os, re
import requests
from subprocess import call

# Define the input file and output directory
input_file = '/home/shuaiw/Methy/benchmark.data.source.txt'
output_dir = '/home/shuaiw/methylation/data/published_data/fanggang/bax'
bam_dir = '/home/shuaiw/methylation/data/published_data/fanggang/bam'
fasta_dir = "/home/shuaiw/methylation/data/published_data/fanggang/fasta/"
align_dir = "/home/shuaiw/methylation/data/published_data/fanggang/align/"

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Function to download a file from a URL
def download_file(url, output_path):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
        print(f"Downloaded {output_path}")
    else:
        print(f"Failed to download {url}")

def download_all():
    # Read the input file and process each URL
    annotation = None
    i = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('###'):
                annotation = line.strip().split('###')[1].strip()
                strain = annotation.split(',')[0].strip().replace(' ', '_').replace(';', '_').replace('(', '_').replace(')', '_')
                print (strain)
            elif line.startswith('http'):
                url = line.strip()
                file_name = os.path.basename(url)
                output_path = os.path.join(output_dir, file_name)
                url2 = url.replace('_p0.bas', '_p0.1.bax')
                file_name2 = os.path.basename(url2)
                output_path2 = os.path.join(output_dir, file_name2)
                url3 = url.replace('_p0.bas', '_p0.2.bax')
                file_name3 = os.path.basename(url3)
                output_path3 = os.path.join(output_dir, file_name3)
                url4 = url.replace('_p0.bas', '_p0.3.bax')
                file_name4 = os.path.basename(url4)
                output_path4 = os.path.join(output_dir, file_name4)


                download_file(url, output_path)
                download_file(url2, output_path2)
                download_file(url3, output_path3)
                download_file(url4, output_path4)

                cmd = f"""
                bax2bam -o {bam_dir}/{strain} {output_dir}/{file_name2} {output_dir}/{file_name3} {output_dir}/{file_name4}
                """ 
                print(cmd)
                os.system(cmd)
            i += 1

    print("All files downloaded and converted.")

def get_fasta():
    file = fasta_dir + '/align_name.tab'
    name_dict = {}
    for line in open(file):
        field = line.strip().split('\t')
        strain = '_'.join(field[0].split())
        name_dict[strain] = fasta_dir + field[1]
    return name_dict

def command(bam, ref, outdir, align_bam, prefix):
    return f"""#!/bin/bash
#SBATCH --job-name=methy     # Job name
#SBATCH --partition=standard # Partition name

    subreads_bam={bam}
    ref={ref}
    outdir={outdir}
    align_bam={align_bam}
    prefix={prefix}

    ## construct output directory if not exists
    if [ ! -d $outdir ]; then
        mkdir $outdir
    fi

    /usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align --preset SUBREAD -j $SLURM_CPUS_ON_NODE $ref $subreads_bam $prefix.raw.bam 
    samtools sort -T $prefix -@ $SLURM_CPUS_ON_NODE -o $align_bam $prefix.raw.bam

    samtools index $align_bam
    /home/shuaiw//smrtlink/pbindex $align_bam

    #### remove the reads aligned with NM > 3

    /usr/bin/time -v -o $prefix.ipdSummary.time \
    ~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers $SLURM_CPUS_ON_NODE \
    --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix


    /usr/bin/time -v -o $prefix.motifmaker.time ~/smrtlink/motifMaker find -f $ref -g $prefix.gff -j $SLURM_CPUS_ON_NODE -o $prefix.motif.csv -m 30

    /usr/bin/time -v -o $prefix.reprocess.time ~/smrtlink/motifMaker reprocess -m $prefix.motif.csv \
    -f $ref -g $prefix.gff -c $prefix.csv -o $prefix.reprocess.gff -j $SLURM_CPUS_ON_NODE
    """

def align_all():
    # Read the input file and process each URL
    annotation = None
    i = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('###'):
                annotation = line.strip().split('###')[1].strip()
                strain = annotation.split(',')[0].strip().replace(' ', '_').replace(';', '_').replace('(', '_').replace(')', '_')
                bam = bam_dir + '/' + strain + '.subreads.bam'
                for name in name_dict:
                    if re.search(name, strain):
                        print (name, strain)
                        fasta = name_dict[name]
                        # cmd = f"samtools faidx {fasta}"
                        # os.system(cmd)
                        align_bam = align_dir + strain + '.align.bam'
                        outdir = align_dir 
                        ref = fasta
                        prefix = align_dir + strain 
                        # script = '/home/shuaiw/borg/methy_workflow8.sh'

                        cmd = command(bam, ref, outdir, align_bam, prefix)
                        with open(f'{prefix}.sh', 'w') as f:
                            f.write(cmd)
                        break
                break


            i += 1





name_dict = get_fasta()
print (name_dict)
align_all()