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

from split_bam import split_bam
from standard_load7 import IPD_load_worker
from generate_contig_list import get_contig_list
from comp_ipd_ratio import get_ipd_ratio
from segment_genome import process_depth_and_gff
from collect_motifs import collect_motifs
from motif_profile import motif_profile_worker
from merge_profile import merge_profile_worker
from cal_invasion_score import batch_MGE_invade


# -------------------------------
# Parameters
# -------------------------------
DEPTH_THRESHOLD = 5
MAX_GAP = 10

motif_maker_bin = "/home/shuaiw/smrtlink/motifMaker"
## pbindex

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Set level for the logger


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
    parser = argparse.ArgumentParser(
        description="Run methylation-based MGE-host linkage discovery pipeline."
    )

    parser.add_argument("--whole_bam", type=str, required=True,
                        help="Input BAM file with kinetic data (HiFi or subreads).")
    parser.add_argument("--whole_ref", type=str, required=True,
                        help="Reference FASTA file for contigs.")
    parser.add_argument("--work_dir", type=str, required=True,
                        help="Working directory for all output files.")
    parser.add_argument("--maxAlignments", type=int, default=100000,
                        help="Maximum number of alignments to process.")
    parser.add_argument("--read_type", choices=["subreads", "hifi"], default="subreads",
                        help="Type of reads in BAM file.")
    parser.add_argument("--max_NM", type=int, default=None,
                        help="Maximum number of mismatches allowed (None to disable).")
    parser.add_argument("--min_len", type=int, default=50000,
                        help="Minimum contig length to process.")
    parser.add_argument("--min_cov", type=int, default=3,
                        help="Minimum read coverage required to retain a base.")
    parser.add_argument("--kmer_mean_db", type=str, default=None,
                        help="Path to optional k-mer mean IPD database.")
    parser.add_argument("--kmer_num_db", type=str, default=None,
                        help="Path to optional k-mer count database.")
    parser.add_argument("--clean", action="store_true", dest="clean", help="Enable cleaning step")
    parser.add_argument("--min_frac", type=float, default=0.5,
                        help="Minimum methylation fraction to retain a motif.")
    parser.add_argument("--min_sites", type=int, default=100,
                        help="Minimum number of methylated sites per motif.")
    parser.add_argument("--min_score", type=int, default=30,
                        help="Minimum score for modification calling.")
    parser.add_argument("--plasmid_file", type=str, default="NA",
                        help="Optional plasmid FASTA file (set to 'NA' if not used).")
    parser.add_argument("--bin_file", type=str, required=False, default=None,
                        help="Path to the binning file containing contig-to-bin mappings.")
    parser.add_argument("--threads", type=int, default=64,
                        help="Number of threads to use for processing.")
    parser.add_argument("--up", type=int, default=7,
                        help="Number of upstream bases to consider for k-mer analysis.")
    parser.add_argument("--down", type=int, default=3,
                        help="Number of downstream bases to consider for k-mer analysis.")
    ## add a new parameter to control the whether to use detect misassembly
    parser.add_argument("--detect_misassembly", action="store_true",
                        help="Enable detection of misassembly in the pipeline.")
    parser.add_argument("--visu_ipd", action="store_true",
                        help="Enable visulization of IPD distribution.")
    parser.add_argument("--binning", action="store_true",
                        help="Enable binning based on methylation.")
    parser.add_argument(
        "--run_steps",
        nargs="+",
        choices=[
            "split", "load", "control", "compare", "motif", "profile", "merge", "host"
        ],
        default=["split", "load", "control", "compare", "motif", "profile", "merge", "host"],
        help="Steps to run in the pipeline (default: all), for easy test."
    )

    

    return parser.parse_args()

def load_ipd_parallel(args, paras):
    logger.info ("Loading IPD data in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for bam in os.listdir(paras["bam_dir"]):
            if not bam.endswith(".bam"):
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
                max_mismatch=args.max_NM,
                read_type=args.read_type,
                ipd_file = ipd_file,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()

            yield finish_code

def get_control_parallele(args, paras):

    get_contig_list(paras["ipd_dir"], paras["ctg_dir"], paras["ctg_list_file"])
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

def motif_worker(ctg_name, bam, fasta, gff, seg_ref, seg_gff, motif, threads, min_score):
    ## run samtools depth
    depth_file = f"{bam}.depth"
    cmd = [
        "samtools", "depth", "-aa", bam
    ]

    with open(depth_file, "w") as depth_output:
        logger.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, stdout=depth_output, stderr=subprocess.DEVNULL, check=True)

    mean_depth = process_depth_and_gff(depth_file, fasta, gff, seg_ref, seg_gff,
                        depth_threshold=DEPTH_THRESHOLD, max_gap=MAX_GAP)

    cmd = [
        motif_maker_bin,
        "find",
        "-f", seg_ref,
        "-g", seg_gff,
        "-o", motif,
        "-j", str(1),
        "-m", str(min_score),
    ]
    logger.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    return ctg_name, mean_depth

def motif_parallel(args, paras):
    logger.info ("Detect motif in parallel...")
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
            seg_ref = os.path.join(paras["segs"], f"{ctg_name}.seg.fa")
            seg_gff = os.path.join(paras["segs"], f"{ctg_name}.seg.gff")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            motifs = os.path.join(paras["motifs"], f"{ctg_name}.motifs.csv")

            future = executor.submit(
                motif_worker,
                ctg_name = ctg_name,
                bam=bam,
                fasta=fasta,
                gff=gff,
                seg_ref=seg_ref,
                seg_gff=seg_gff,
                motif = motifs,
                threads = args.threads,
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
            if os.path.exists(profile):
                # logger.info(f"Profile file {profile} already exists. Skipping...")
                continue

            
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
            plasmid_file = args.plasmid_file,
            bin_file = args.bin_file,
            threads = args.threads,
            bin_flag = args.binning,
        )

def depth_analysis(paras, ctg_depth_dict):
    ## change the dict to df
    depth_df = (
        pd.DataFrame.from_dict(ctg_depth_dict, orient="index", columns=["depth"])
        .reset_index()
        .rename(columns={"index": "contig"})
    )
    depth_df.to_csv(paras["depth_file"], index=False)
    ## plot the depth distribution
    sns.set(style="whitegrid")
    sns.histplot(depth_df, x="depth")
    ## save the plot

    plt.savefig(paras["depth_plot"])

def predict_host_worker(args, paras):
    os.makedirs(paras["hosts"], exist_ok = True)
    if args.plasmid_file != 'NA' and os.path.exists(args.plasmid_file):
        batch_MGE_invade(args.plasmid_file, paras["profiles"], paras["hosts"], bin_file=args.bin_file, min_frac=0.5, threads=args.threads)

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

    paras["kmer_bin"] = os.path.join(sys.path[0], "src", "test")
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
    os.makedirs(paras["segs"], exist_ok=True)
    if args.binning:
        os.makedirs(paras["bins"], exist_ok=True)
    os.makedirs(paras["profiles"], exist_ok=True)

    return paras

if __name__ == "__main__":

    args = parse_arguments()
    ctg_depth_dict = {}

    if args.max_NM is None:
        if args.read_type == "hifi":
            args.max_NM = 200
        else:
            args.max_NM = 1000

    logger.info("🔬 Running MGE-host linkage pipeline with the following parameters:")

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
    logger.info("\n[Placeholder] Pipeline execution starts here...\n")



    if "split" in args.run_steps:
        record_resource_usage(
            "Splitting BAM files",
            split_bam,
            args.whole_bam, args.work_dir, args.whole_ref, args.threads, args.min_len, args.max_NM
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
        ctg_depth_dict = record_resource_usage(
            "Motif identification",
            motif_parallel,
            args, paras
        )
        record_resource_usage(
            "Depth analysis",
            depth_analysis,
            paras, ctg_depth_dict
        )
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
    
    if args.clean:
        logger.info("Cleaning up intermediate files...")
        for folder in [paras["bam_dir"], paras["ctg_dir"], paras["ipd_dir"], paras["control"], paras["ipd_ratio"], paras["segs"]]:
            if os.path.exists(folder):
                shutil.rmtree(folder)


    logger.info("\n[Placeholder] Pipeline execution ends here...\n")
    logger.info("🔬 Pipeline execution completed.")

# if __name__ == "__main__":
#     main()
