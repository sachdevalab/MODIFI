### given a bam, extract the read alignment for each contig, and store it in a separate bam file


import pysam
import os
from multiprocessing import Pool
import sys


def split_bam(bam, split_bam_dir):
    # Ensure the output directory exists
    os.makedirs(split_bam_dir, exist_ok=True)
    samfile = pysam.AlignmentFile(bam, "rb")
    i = 1
    args = []
    for contig in samfile.references:
        ## only consider the contigs that are longer than 1M
        # if samfile.lengths[samfile.get_tid(contig)] < 500000:
        #     continue
        contig_bam = split_bam_dir + contig + ".bam"
        ref = split_bam_dir + contig + ".fasta"
        ### each contig is a task, run them with multiprocessing
        # handle_each_contig(contig,ref,contig_bam,samfile)
        args.append((contig,ref,contig_bam, bam))
        i += 1
    samfile.close()
    # Use multiprocessing to handle each contig in parallel
    with Pool(processes=10) as pool:
        pool.starmap(handle_each_contig, args)

def handle_each_contig(contig,ref,contig_bam,bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    os.system(f"samtools faidx {whole_ref} {contig} > {ref}")
    os.system(f"samtools faidx {ref}")

    # Create a new header that only includes the specific contig
    new_header = samfile.header.to_dict().copy()
    
    new_header['SQ'] = [sq for sq in new_header['SQ'] if sq['SN'] == contig]
    # print (new_header['SQ'])

    contig_samfile = pysam.AlignmentFile(contig_bam, "wb", header=new_header)
    for read in samfile.fetch(contig):
        read.reference_id = 0
        contig_samfile.write(read)
    contig_samfile.close()
    samfile.close()
    detect_methy(contig_bam, ref, contig)


def detect_methy(bam, ref, contig, threads=1):
    """
    Given each bam, index it and detect the methylation
    """
    cmd = f"""
    align_bam={bam}
    ref={ref}
    prefix={split_bam_dir}/{contig}
    threads={threads}

    /home/shuaiw//smrtlink/pbindex $align_bam
    ~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers $threads \
    --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix


    ~/smrtlink/motifMaker find -f $ref -g $prefix.gff -j $threads -o $prefix.motif.csv -m 30

    #  ~/smrtlink/motifMaker reprocess -m $prefix.motif.csv \
    # -f $ref -g $prefix.gff -c $prefix.csv -o $prefix.reprocess.gff 

    # rm $align_bam
    # rm $align_bam.pbi
    rm $prefix.pickle
    # rm $prefix.fasta*
    # rm $prefix.csv
    # rm $prefix.m5C.gff
    """
    os.system(cmd)


# if len(sys.argv) < 3:
#     print("Usage: python split_bam_single.py <arg1> <arg2>", sys.argv, len(sys.argv))
#     sys.exit(1)
### 36608 contigs   25575 are detected
bam = "/home/shuaiw/methylation/data/borg/all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
whole_ref="/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
split_bam_dir = "/home/shuaiw/methylation/data/borg/split_bam_dir3/"
# split_bam(bam, split_bam_dir)
# bam = sys.argv[1]
# ref = sys.argv[2]
# contig = sys.argv[3]
all_params = sys.argv[1]
para_list = all_params.split()
bam = para_list[0]
ref = para_list[1]
contig = para_list[2]
# print (bam, ref, contig)
detect_methy(bam, ref, contig)