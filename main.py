import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import os
import subprocess
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import time
import psutil
import logging
import shutil
import yaml

## add the path to scripts
sys.path.append(os.path.join(sys.path[0], "scripts"))

from split_bam import split_bam
from standard_load7 import IPD_load_worker
from generate_contig_list import get_contig_list
from comp_ipd_ratio import get_ipd_ratio
from segment_genome import process_depth_and_gff
from collect_motifs import collect_motifs
from motif_profile import motif_profile_worker
from merge_profile import merge_profile_worker
from estimate_linkage import batch_MGE_invade
from analyze_RM import RM_main

from load_cfg import load_binaries

# -------------------------------
# Version
# -------------------------------
__version__ = "0.0.0"

# -------------------------------
# Parameters
# -------------------------------
DEPTH_THRESHOLD = 5
MAX_GAP = 10

# motif_maker_bin = "/home/shuaiw/smrtlink/motifMaker"
# motif_maker_bin = "/home/shuaiw/smrtlink/pbmotifmaker"
# pbmm2_bin = "/home/shuaiw/smrtlink/pbmm2"
# pbindex_bin = "/home/shuaiw/smrtlink/pbindex"
## pbindex

motif_maker_bin = None
pbmm2_bin = None
pbindex_bin = None

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Set level for the logger


# config_file = os.path.join(sys.path[0], "config.yaml")
# if os.path.exists(config_file):
#     with open(config_file, "r") as cf:
#         config = yaml.safe_load(cf)
#         motif_maker_bin = os.path.join(config.get("smrtlink_bin", ""), "pbmotifmaker")
#         pbmm2_bin = os.path.join(config.get("smrtlink_bin", ""), "pbmm2")
#         pbindex_bin = os.path.join(config.get("smrtlink_bin", ""), "pbindex")
# print (f"Loaded binaries from config.yaml: motif_maker_bin={motif_maker_bin}, pbmm2_bin={pbmm2_bin}, pbindex_bin={pbindex_bin}")

def record_resource_usage(step_name, func, *args, **kwargs):
    """
    Record CPU time and peak memory usage for a given function.

    Parameters:
        step_name (str): Name of the step being executed.
        func (callable): The function to execute.
        *args: Positional arguments for the function.
        **kwargs: Keyword arguments for the function.

    Returns:
        The return value of the function.
    """
    logger.info(f"🔄 Starting step: {step_name}")
    start_time = time.time()  # Start wall-clock time
    process = psutil.Process()

    # Run the function
    result = func(*args, **kwargs)

    # Record resource usage
    end_time = time.time()  # End wall-clock time
    elapsed_time = (end_time - start_time) / 3600  # Convert to hours
    peak_memory = process.memory_info().rss / (1024 * 1024)  # Convert to MB

    logger.info(f"✅ Step '{step_name}' completed.")
    logger.info(f"   ⏱️ Wall-Clock Time: {elapsed_time:.2f} hours")
    logger.info(f"   📈 Peak Memory: {peak_memory:.2f} MB\n")

    return result

def parse_arguments():
    # Custom formatter to combine default values and raw description formatting
    class CustomFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    
    parser = argparse.ArgumentParser(
        description="""DNA modification detection and MGE-host linkage inference.

Example: python main.py \\
        --work_dir /path-to/output \\
        --unaligned_bam raw.hifi_reads.bam \\
        --whole_ref metagenomic_assembly.fasta \\
        --read_type hifi
""",
        formatter_class=CustomFormatter
    )

    parser.add_argument("-v", "--version", action="version", version=f"MODIFI v{__version__}")
    parser.add_argument("--whole_bam", type=str, required=False,
                        help="Input aligned BAM file with kinetic data (HiFi or subreads). Use this for pre-aligned BAM files.\
                              Must be aligned by pbmm2 to ensure kinetic tags are present.")
    parser.add_argument("--unaligned_bam", type=str, required=False,
                        help="Input unaligned BAM file with kinetic data (HiFi or subreads). Will be aligned using pbmm2.")
    parser.add_argument("--whole_ref", type=str, required=True,
                        help="Reference FASTA file for contigs.")
    parser.add_argument("--work_dir", type=str, required=True,
                        help="Working directory for all output files.")
    # parser.add_argument("--maxAlignments", type=int, default=100000,
    #                     help="Maximum number of alignments to process.")
    parser.add_argument("--read_type", choices=["subreads", "hifi"], default="hifi",
                        help="Type of reads in BAM file.")
    parser.add_argument("--min_iden", type=float, default=None,
                        help="Minimum identity allowed for read alignment (default: 0.97 for hifi, 0.85 for subreads).")
    # parser.add_argument("--max_NM", type=int, default=None,
    #                     help="Maximum number of mismatches allowed (None to disable).")
    parser.add_argument("--min_len", type=int, default=1000,
                        help="Minimum contig length to process.")
    parser.add_argument("--min_cov", type=int, default=1,
                        help="Minimum read coverage required to retain a base.")
    parser.add_argument("--min_ctg_cov", type=int, default=5,
                        help="Minimum read coverage required to retain a contig.")
    parser.add_argument("--kmer_mean_db", type=str, default=None,
                        help="Path to optional k-mer mean IPD database.")
    parser.add_argument("--kmer_num_db", type=str, default=None,
                        help="Path to optional k-mer count database.")
    parser.add_argument("--min_frac", type=float, default=0.4,
                        help="Minimum methylation fraction to retain a motif.")
    parser.add_argument("--min_sites", type=int, default=30,
                        help="Minimum number of methylated sites per motif.")
    parser.add_argument("--min_score", type=int, default=30,
                        help="Minimum modification score for motif calling.")
    parser.add_argument("--mge_file", type=str, default="NA",
                        help="MGE table file (sep by tab), can be output of geNomad, with at least one column with header: seq_name.")
    parser.add_argument("--threads", type=int, default=64,
                        help="Number of threads to use for processing.")
    parser.add_argument("--up", type=int, default=7,
                        help="Number of upstream bases to consider for k-mer analysis.")
    parser.add_argument("--down", type=int, default=3,
                        help="Number of downstream bases to consider for k-mer analysis.")
    parser.add_argument("--clean", action="store_true", dest="clean", help="Enable cleaning step")
    parser.add_argument("--bin_file", type=str, required=False, default=None,
                        help="Path to the binning file containing contig-to-bin mappings. Sep by tab, with two columns: contig_name and bin_name, no headers. ")
    parser.add_argument("--segment", action="store_true", dest="segment", 
                        help="Enable segmentation of the contigs by depth, increase recall for low-depth contigs, but cost more time.")
    parser.add_argument("--visu_ipd", action="store_true",
                        help="Enable visulization of IPD distribution.")
    parser.add_argument("--annotate_rm", action="store_true",
                        help="Enable RM system annotation, requiring MicrobeMod installed in system path.")
    parser.add_argument("--rm_gene_file", type=str,
                        help="RM gene annotation file by MicrobeMod, with suffix .rm.genes.tsv (only for testing)")
    ## add a new parameter to control the whether to use detect misassembly
    parser.add_argument("--detect_misassembly", action="store_true",
                        help="Enable manual detection of misassembly by IPD ratio line continuity.")
    parser.add_argument("--binning", action="store_true",
                        help="Enable binning based on methylation (in testing).")
    parser.add_argument(
        "--run_steps",
        nargs="+",
        choices=[
            "split", "load", "control", "compare", "motif", "profile", "merge", "host", "anno"
        ],
        default=["split", "load", "control", "compare", "motif", "profile", "merge", "host"],
        help="Steps to run in the pipeline (default: all), sep by space. Only for easy testing."
    )

    args = parser.parse_args()
    
    # Set min_iden default based on read_type if not explicitly provided
    if args.min_iden is None:
        if args.read_type.lower() == "hifi":
            args.min_iden = 0.97
        elif args.read_type.lower() == "subreads":
            args.min_iden = 0.85
        else:
            raise ValueError(f"Unknown read type: {args.read_type}. Must be 'hifi' or 'subreads'")
    
    # Validate that either whole_bam or unaligned_bam is provided, but not both
    if args.whole_bam and args.unaligned_bam:
        parser.error("Cannot specify both --whole_bam and --unaligned_bam. Use one or the other.")
    if not args.whole_bam and not args.unaligned_bam:
        parser.error("Must specify either --whole_bam (for aligned BAM) or --unaligned_bam (for unaligned BAM).")
    
    return args

def align_bam_with_pbmm2(input_bam, reference_fasta, output_bam, read_type, threads=1):
    """
    Align unaligned BAM file to reference using pbmm2.
    
    Args:
        input_bam: Path to input unaligned BAM file
        reference_fasta: Path to reference FASTA file
        output_bam: Path to output aligned BAM file
        read_type: 'hifi' or 'subreads'
        threads: Number of threads to use
    """
    logger.info(f"Aligning {input_bam} to {reference_fasta} using pbmm2...")
    
    # Set preset based on read type
    if read_type.lower() == "hifi":
        preset = "CCS"
    elif read_type.lower() == "subreads":
        preset = "SUBREAD"
    else:
        raise ValueError(f"Unknown read type: {read_type}. Must be 'hifi' or 'subreads'")
    
    # Create temporary raw BAM file
    temp_raw_bam = output_bam.replace(".bam", ".raw.bam")
    
    try:
        # Run pbmm2 alignment
        cmd = [
            pbmm2_bin,
            "align",
            "--preset", preset,
            "-j", str(threads),
            reference_fasta,
            input_bam,
            temp_raw_bam
        ]
        
        logger.info(f"Running pbmm2 command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        
        # Sort the aligned BAM file
        logger.info(f"Sorting aligned BAM file...")
        sort_cmd = [
            "samtools", "sort",
            "-T", output_bam.replace(".bam", ""),
            "-@", str(threads),
            "-o", output_bam,
            temp_raw_bam
        ]
        
        subprocess.run(sort_cmd, check=True)
        
        # Index the sorted BAM file
        logger.info(f"Indexing aligned BAM file...")
        subprocess.run(["samtools", "index", output_bam], check=True)
        subprocess.run([pbindex_bin, output_bam], check=True)
        
        # Clean up temporary file
        if os.path.exists(temp_raw_bam):
            os.remove(temp_raw_bam)
            
        logger.info(f"Successfully aligned and sorted BAM file: {output_bam}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during alignment: {e}")
        # Clean up temporary files on error
        for temp_file in [temp_raw_bam]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        raise
    except Exception as e:
        logger.error(f"Unexpected error during alignment: {e}")
        raise

def load_ipd_parallel(args, paras):
    logger.info ("Loading IPD data in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for bam in os.listdir(paras["bam_dir"]):
            if not bam.endswith(".bam"):
                continue
            if bam.endswith(".filtered.bam"):
                continue
            bam_path = os.path.join(paras["bam_dir"], bam)
            ctg_name = bam.replace(".bam", "")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            # outputfile = os.path.join(paras["ipd_dir"], f"{ctg_name}.count")
            ipd_file = os.path.join(paras["ipd_dir"], f"{ctg_name}.ipd1.csv")

            ## if the output file already exists and not empty, skip it
            if os.path.exists(ipd_file) and os.path.getsize(ipd_file) > 0:
                # logger.info(f"Output file {outputfile} already exists. Skipping...")
                continue

            future = executor.submit(
                IPD_load_worker,
                fasta=fasta,
                subread_bam=bam_path,
                max_mismatch=100000000,
                read_type=args.read_type,
                ipd_file = ipd_file,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()

            yield finish_code

def get_control_parallele(args, paras):

    ctg_num = get_contig_list(paras["ipd_dir"], paras["ctg_dir"], paras["ctg_list_file"])
    ## if ctg_num is 0, report the ctg num is 0 and exit
    if ctg_num == 0:
        logger.warning("No contigs found with corresponding IPD files. Exiting.")
        sys.exit(0)
    cmd = [
        paras["kmer_bin"],
        "--ipd_dir", paras["ipd_dir"],
        "--control_dir", paras["control"],
        "--fasta_list", paras["ctg_list_file"],
        "--thread_num", str(args.threads),
        "--kmer_mean_file", str(args.kmer_mean_db),
        "--kmer_num_file", str(args.kmer_num_db),
        "--up", str(args.up),
        "--down", str(args.down),
    ]

    logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def safe_get_ipd_ratio(**kwargs):
    try:
        return get_ipd_ratio(**kwargs)
    except Exception as e:
        logger.error(f"Failed on {kwargs.get('csv')}: {e}")
        return None


def compare_ipd_parallel(args, paras):
    logger.info ("Detect modified bases in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for bam in os.listdir(paras["ipd_dir"]):
            if not bam.endswith(".ipd1.csv"):
                continue
            ctg_name = bam.replace(".ipd1.csv", "")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            control = os.path.join(paras["control"], f"{ctg_name}.ipd2.csv")
            ipd_ratio_file = os.path.join(paras["ipd_ratio"], f"{ctg_name}.ipd3.csv")
            gff = os.path.join(paras["gffs"], f"{ctg_name}.gff")
            figure_file = os.path.join(paras["figs"], f"{ctg_name}.png")

            future = executor.submit(
                safe_get_ipd_ratio,
                csv = control,
                output = ipd_ratio_file,
                gff = gff,
                ref = fasta,
                figure_file = figure_file,
                min_cov = args.min_cov,
                visu_flag = args.visu_ipd,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()

            yield finish_code

def depth_worker(bam, fasta, gff, seg_ref, seg_gff,):
    ## run samtools depth
    depth_file = f"{bam}.depth"
    cmd = [
        "samtools", "depth", "-aa", bam
    ]

    with open(depth_file, "w") as depth_output:
        # logger.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, stdout=depth_output, stderr=subprocess.DEVNULL, check=True)

    mean_depth = process_depth_and_gff(depth_file, fasta, gff, seg_ref, seg_gff,
                        depth_threshold=DEPTH_THRESHOLD, max_gap=MAX_GAP)
    return mean_depth

def depth_parallel(args, paras):
    logger.info ("segment by depth in parallel...")

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for gff in os.listdir(paras["gffs"]):
            if not gff.endswith(".gff"):
                continue
            if gff.endswith(".reprocess.gff"):
                continue
            # logger.info ("####", paras["gffs"], gff)
            ctg_name = gff.replace(".gff", "")
            bam = os.path.join(paras["bam_dir"], f"{ctg_name}.bam")
            gff = os.path.join(paras["gffs"], gff)
            seg_ref = os.path.join(paras["segs"], f"{ctg_name}.seg.fa")
            seg_gff = os.path.join(paras["segs"], f"{ctg_name}.seg.gff")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            motifs = os.path.join(paras["motifs"], f"{ctg_name}.motifs.csv")

            future = executor.submit(
                depth_worker,
                bam=bam,
                fasta=fasta,
                gff=gff,
                seg_ref=seg_ref,
                seg_gff=seg_gff,
            )
            futures.append(future)

        for future in tqdm(as_completed(futures), total=len(futures)):
            mean_depth = future.result()


def motif_worker(ctg_name, seg_ref, seg_gff, motif, min_score, each_thread=1):
    ## run samtools depth
    # depth_file = f"{bam}.depth"
    # cmd = [
    #     "samtools", "depth", "-aa", bam
    # ]

    # with open(depth_file, "w") as depth_output:
    #     # logger.info(f"Running command: {' '.join(cmd)}")
    #     subprocess.run(cmd, stdout=depth_output, stderr=subprocess.DEVNULL, check=True)

    # mean_depth = process_depth_and_gff(depth_file, fasta, gff, seg_ref, seg_gff,
    #                     depth_threshold=DEPTH_THRESHOLD, max_gap=MAX_GAP)
    cmd = [   ## old version
        motif_maker_bin,
        "find",
        "-f", seg_ref,
        "-g", seg_gff,
        "-o", motif,
        "-j", str(each_thread),
        "-m", str(min_score),
    ]

    cmd = [   ## new version 1.2.0
        motif_maker_bin,
        "find",
        seg_ref,
        seg_gff,
        motif,
        "-j", str(each_thread),
        "--min-score", str(min_score),
    ]
    if motif_maker_bin.endswith(".jar"):
        print ("Using MultiMotifMaker for motif finding")
        cmd = ["java", "-jar", motif_maker_bin, 
               "find", 
               "-f", seg_ref, 
               "-g", seg_gff,
               "-o", motif, 
                "-t", str(each_thread), 
                "-m", str(min_score),
        ]

    # logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return ctg_name, 0

def motif_parallel(args, paras):
    logger.info(f"Detect motif in parallel (memory-aware)... Using {args.threads} threads.")

    MAX_MEMORY_USAGE_RATIO = 0.9  # Don't exceed 90% of total memory
    CHECK_INTERVAL = 1  # seconds between memory checks

    def memory_safe():
        mem = psutil.virtual_memory()
        return mem.available / mem.total > (1 - MAX_MEMORY_USAGE_RATIO)
    
    ctg_depth_dict = {}
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for gff in os.listdir(paras["gffs"]):
            if not gff.endswith(".gff"):
                continue
            if gff.endswith(".reprocess.gff"):
                continue
            # logger.info ("####", paras["gffs"], gff)
            ctg_name = gff.replace(".gff", "")
            bam = os.path.join(paras["bam_dir"], f"{ctg_name}.bam")
            gff = os.path.join(paras["gffs"], gff)
            if args.segment: ## enable segment
                seg_ref = os.path.join(paras["segs"], f"{ctg_name}.seg.fa")
                seg_gff = os.path.join(paras["segs"], f"{ctg_name}.seg.gff")
            else:
                seg_ref = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
                seg_gff = os.path.join(paras["gffs"], f"{ctg_name}.gff")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            motifs = os.path.join(paras["motifs"], f"{ctg_name}.motifs.csv")

            # Memory-aware submission
            while not memory_safe():
                logger.debug("Waiting for memory to free up...")
                time.sleep(CHECK_INTERVAL)

            future = executor.submit(
                motif_worker,
                ctg_name = ctg_name,
                seg_ref=seg_ref,
                seg_gff=seg_gff,
                motif = motifs,
                min_score = args.min_score,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            ctg_name, mean_depth = future.result()
            # logger.info ("*****", ctg_name, mean_depth)
            ctg_depth_dict[ctg_name] = mean_depth
            # yield finish_code
    # logger.info (ctg_depth_dict)
    return ctg_depth_dict
    
def collect_motifs_worker(args, paras):
    collect_motifs(
        folder=paras["motifs"],
        all_motif=paras["all_motifs"],
        MIN_FRAC=args.min_frac,
        MIN_detect=args.min_sites,)
    logger.info ("Motif collection done.")

def profile_parallel(args, paras):
    logger.info ("Profile motif in parallel...")
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        all_motifs = pd.read_csv(paras["all_motifs"]) 
        for gff in os.listdir(paras["gffs"]):
            if not gff.endswith(".gff"):
                continue
            if gff.endswith(".reprocess.gff"):
                continue
            # logger.info ("####", paras["gffs"], gff)
            ctg_name = gff.replace(".gff", "")
            bam = os.path.join(paras["bam_dir"], f"{ctg_name}.bam")
            gff = os.path.join(paras["gffs"], gff)
            seg_ref = os.path.join(paras["segs"], f"{ctg_name}.seg.fa")
            seg_gff = os.path.join(paras["segs"], f"{ctg_name}.seg.gff")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            motifs = os.path.join(paras["motifs"], f"{ctg_name}.motifs.csv")
            profile = os.path.join(paras["profiles"], f"{ctg_name}.motifs.profile.csv")
            ipd_ratio = os.path.join(paras["ipd_ratio"], f"{ctg_name}.ipd3.csv")

            ## skip if the profile file already exists
            # if os.path.exists(profile):
            #     logger.info(f"Profile file {profile} already exists. Skipping...")
            #     continue

            
            future = executor.submit(
                motif_profile_worker,
                my_ref = fasta,
                gff = gff,
                all_motifs = all_motifs,
                profile = profile,
                ipd_ratio_file = ipd_ratio,
                min_frac = args.min_frac,
                min_sites = args.min_sites,
                score_cutoff = args.min_score,
                min_cov = args.min_cov,
                misassembly = args.detect_misassembly,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()
            yield finish_code

            # print (finish_code)

def run_merge_profile(args, paras):
        merge_profile_worker(
            work_dir = args.work_dir, 
            heat_map = paras["heatmap"],
            profile_list = paras["profiles"], 
            total_profile = paras["total_profile"], 
            min_frac = args.min_frac, 
            whole_ref = args.whole_ref, 
            bin_file = args.bin_file,
            threads = args.threads,
            bin_flag = args.binning,
        )

def depth_analysis(paras, ctg_depth_dict):
    ## change the dict to df
    depth_df = (
        pd.DataFrame.from_dict(ctg_depth_dict, orient="index", columns=["depth", "length"])
        .reset_index()
        .rename(columns={"index": "contig"})
    )
    logger.info(f"{len(depth_df)} contigs in the assembly.")
    over0_depth = depth_df[depth_df["depth"] > 0]
    logger.info(f"{len(over0_depth)} contigs with depth > 0.")
    ### report mean and median depth
    mean_depth = over0_depth["depth"].mean()
    median_depth = over0_depth["depth"].median()
    logger.info(f"Mean depth: {mean_depth:.2f}, Median depth: {median_depth:.2f} in {len(over0_depth)} contigs with depth > 0.")
    ## sort depth_df by depth in reverse
    depth_df = depth_df.sort_values(by="depth", ascending=False)
    depth_df.to_csv(paras["depth_file"], index=False)
    ## plot the depth distribution
    sns.set(style="whitegrid")
    # Remove the lowest and highest 10% contigs by depth for plotting
    lower_quantile = over0_depth["depth"].quantile(0.01)
    upper_quantile = over0_depth["depth"].quantile(0.99)
    filtered_depth = over0_depth[(over0_depth["depth"] >= lower_quantile) & (over0_depth["depth"] <= upper_quantile)]
    sns.histplot(filtered_depth, x="depth")
    plt.savefig(paras["depth_plot"])
    ## log the number of contigs with depth larger than min_ctg_cov
    logger.info(f"{len(over0_depth[over0_depth['depth'] >= args.min_ctg_cov])} contigs with depth >= {args.min_ctg_cov}.")

def predict_host_worker(args, paras):
    os.makedirs(paras["hosts"], exist_ok = True)
    if args.mge_file != 'NA' and os.path.exists(args.mge_file):
        batch_MGE_invade(args.mge_file, paras["profiles"], paras["hosts"], args.whole_ref, bin_file=args.bin_file, min_frac=args.min_frac, threads=args.threads, min_ctg_cov = args.min_ctg_cov, min_detect = args.min_sites)

def get_paras(args):
    """
    Get parameters for the pipeline.
    """

    paras = {}
    paras["bam_dir"] = os.path.join(args.work_dir, "bams")
    paras["ctg_dir"] = os.path.join(args.work_dir, "contigs")
    paras["ipd_dir"] = os.path.join(args.work_dir, "ipd")
    paras["control"] = os.path.join(args.work_dir, "control")
    paras["ipd_ratio"] = os.path.join(args.work_dir, "ipd_ratio")
    paras["gffs"] = os.path.join(args.work_dir, "gffs")
    paras["figs"] = os.path.join(args.work_dir, "figs")
    paras["motifs"] = os.path.join(args.work_dir, "motifs")
    paras["segs"] = os.path.join(args.work_dir, "segs")
    paras["bins"] = os.path.join(args.work_dir, "bins")
    paras["profiles"] = os.path.join(args.work_dir, "profiles")
    paras["hosts"] = os.path.join(args.work_dir, "hosts")
    paras["RM"] = os.path.join(args.work_dir, "RM_systems")
    paras["RM_prefix"] = os.path.join(paras["RM"], "all_ctgs_RM")
    paras["RM_genes"] = paras["RM_prefix"] + ".rm.genes.tsv"

    paras["kmer_bin"] = os.path.join(sys.path[0], "src", "get_control_IPD")
    paras["ctg_list_file"] = os.path.join(args.work_dir, "contigs_list.txt")
    paras["all_motifs"] = os.path.join(args.work_dir, "all.motifs.csv")
    paras["total_profile"] = os.path.join(args.work_dir, "motif_profile.csv")
    paras["heatmap"] = os.path.join(args.work_dir, "motif_heatmap.pdf")
    paras["depth_file"] = os.path.join(args.work_dir, "mean_depth.csv")
    paras["depth_plot"] = os.path.join(args.work_dir, "depth_distribution.pdf")
    paras["log"] = os.path.join(args.work_dir, "log.txt")

    ## build the directory if not exist
    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(paras["bam_dir"], exist_ok=True)
    os.makedirs(paras["ctg_dir"], exist_ok=True)
    os.makedirs(paras["ipd_dir"], exist_ok=True)
    os.makedirs(paras["control"], exist_ok=True)
    os.makedirs(paras["ipd_ratio"], exist_ok=True)
    os.makedirs(paras["gffs"], exist_ok=True)
    if args.visu_ipd:
        os.makedirs(paras["figs"], exist_ok=True)
    os.makedirs(paras["motifs"], exist_ok=True)
    if args.segment:
        os.makedirs(paras["segs"], exist_ok=True)
    if args.binning:
        os.makedirs(paras["bins"], exist_ok=True)
    os.makedirs(paras["profiles"], exist_ok=True)

    return paras

if __name__ == "__main__":

    args = parse_arguments()
    motif_maker_bin, pbmm2_bin, pbindex_bin = load_binaries()
    ctg_depth_dict = {}


    logger.info("🔬 Running pipeline with the following parameters:")

    paras = get_paras(args)


    # logging.basicConfig(filename = paras["log"],\
    # format='[%(asctime)s-%(filename)s-%(levelname)s:%(message)s]', level = logging.INFO,filemode='w')

    # Create formatter
    formatter = logging.Formatter('[%(asctime)s - %(filename)s - %(levelname)s: %(message)s]')

    # === File handler ===
    file_handler = logging.FileHandler(paras["log"], mode='w')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.INFO)

    # === Console handler ===
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.INFO)

    # === Clear previous handlers to prevent duplicates ===
    if logger.hasHandlers():
        logger.handlers.clear()

    # Add both handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # === Example usage ===
    logger.info("Logging is now working to both file and terminal.")

    ## output all the parameters to the log file
    logger.info("Pipeline parameters:")
    for k, v in vars(args).items():
        logger.info(f"  {k}: {v}")


    # === Insert your pipeline logic below ===
    logger.info("\n🔬 Pipeline execution starts here...\n")

    # === Handle unaligned BAM alignment if needed ===
    if args.unaligned_bam:
        logger.info("Unaligned BAM provided. Running alignment with pbmm2...")
        
        # Create align_bam subfolder
        align_bam_dir = os.path.join(args.work_dir, "align_bam")
        os.makedirs(align_bam_dir, exist_ok=True)
        aligned_bam_path = os.path.join(align_bam_dir, "aligned.bam")
        if not os.path.exists(aligned_bam_path):
            record_resource_usage(
                "Aligning BAM with pbmm2",
                align_bam_with_pbmm2,
                args.unaligned_bam, args.whole_ref, aligned_bam_path, args.read_type, args.threads
            )
        else:
            logger.info(f"Aligned BAM already exists at {aligned_bam_path}, skipping alignment.")
        
        # Set the aligned BAM as the whole_bam for the rest of the pipeline
        args.whole_bam = aligned_bam_path
        logger.info(f"Using aligned BAM file for pipeline: {args.whole_bam}")
    else:
        logger.info("Using provided aligned BAM file for pipeline.")

    if "split" in args.run_steps:
        ctg_depth_dict = record_resource_usage(
            "Splitting BAM files",
            split_bam,
            args.whole_bam, args.work_dir, args.whole_ref, pbindex_bin, args.threads, args.min_len, 100000000, args.min_ctg_cov, args.min_iden
        )
        record_resource_usage(
            "Depth analysis",
            depth_analysis,
            paras, ctg_depth_dict
        )
        
    if "load" in args.run_steps:
        record_resource_usage(
            "Loading IPD data",
            lambda: list(load_ipd_parallel(args, paras))  # Wrap generator in a list to execute fully
        )

    if "control" in args.run_steps:
        record_resource_usage(
            "Generating control files",
            get_control_parallele,
            args, paras
        )

    if "compare" in args.run_steps:
        record_resource_usage(
            "Calculating IPD ratios",
            lambda: list(compare_ipd_parallel(args, paras))  # Wrap generator in a list to execute fully
        )

    if "motif" in args.run_steps:
        ## first depth
        if args.segment:
            logger.info("segment contigs in parallel...")
            record_resource_usage(
                "segment contigs by depth",
                depth_parallel,
                args, paras
            )
        logger.info("Identifying motifs in parallel...")
        ctg_depth_dict = record_resource_usage(
            "Motif identification",
            motif_parallel,
            args, paras
        )
        # record_resource_usage(
        #     "Depth analysis",
        #     depth_analysis,
        #     paras, ctg_depth_dict
        # )
    if "profile" in args.run_steps:
        record_resource_usage(
            "Collecting motifs",
            collect_motifs_worker,
            args, paras
        )
        record_resource_usage(
            "Profiling motifs",
            lambda: list(profile_parallel(args, paras))  
        )



    if "merge" in args.run_steps:
        record_resource_usage(
            "Merging motif profiles",
            run_merge_profile,
            args, paras
        )
    
    if "host" in args.run_steps:
        record_resource_usage(
            "Predicting host",
            predict_host_worker,
            args, paras
        )

    if args.annotate_rm or "anno" in args.run_steps:
        logger.info("Annotating RM systems...")
        ## check if paras["RM"] exists, if not, create it
        os.makedirs(paras["RM"], exist_ok=True)
        cmd = [
            "MicrobeMod", "annotate_rm",
            "-f", args.whole_ref,  # or your desired fasta file
            "-o", paras["RM_prefix"],
            "-t", str(args.threads)
        ]
        if args.rm_gene_file:
            print (f"Using RM gene annotation file: {args.rm_gene_file}")
            paras["RM_genes"] = args.rm_gene_file
        else:
            logger.info(f"Running command: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)

        RM_main(
            anno_file=paras["RM_genes"],
            motif_dir=paras["motifs"],
            RM_dir=paras["RM"],
            bin_file=args.bin_file
        )
    
    if args.clean:
        logger.info("Cleaning up intermediate files...")
        for folder in [paras["bam_dir"], paras["ctg_dir"], paras["ipd_dir"], paras["control"], paras["ipd_ratio"], paras["segs"]]:
            if os.path.exists(folder):
                shutil.rmtree(folder)


    logger.info("\n[Placeholder] Pipeline execution ends here...\n")
    logger.info("🔬 Pipeline execution completed.")

# if __name__ == "__main__":
#     main()
