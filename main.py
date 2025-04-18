import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import os
import subprocess
import sys

from split_bam import split_bam
from standard_load7 import IPD_load_worker
from generate_contig_list import get_contig_list
from comp_ipd_ratio import get_ipd_ratio
from segment_genome import process_depth_and_gff
from collect_motifs import collect_motifs
from motif_profile import motif_profile_worker
from merge_profile import merge_profile_worker


# -------------------------------
# Parameters
# -------------------------------
DEPTH_THRESHOLD = 5
MAX_GAP = 10

motif_maker_bin = "/home/shuaiw/smrtlink/motifMaker"
## pbindex


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
    parser.add_argument("--no-clean", action="store_false", dest="clean", help="Disable cleaning step")
    parser.add_argument("--min_frac", type=float, default=0.5,
                        help="Minimum methylation fraction to retain a motif.")
    parser.add_argument("--min_sites", type=int, default=100,
                        help="Minimum number of methylated sites per motif.")
    parser.add_argument("--min_score", type=int, default=30,
                        help="Minimum score for modification calling.")
    parser.add_argument("--plasmid_file", type=str, default="NA",
                        help="Optional plasmid FASTA file (set to 'NA' if not used).")
    parser.add_argument("--threads", type=int, default=64,
                        help="Number of threads to use for processing.")

    return parser.parse_args()

def load_ipd_parallel(args, paras):
    print ("Loading IPD data in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for bam in os.listdir(paras["bam_dir"]):
            if not bam.endswith(".bam"):
                continue
            bam_path = os.path.join(paras["bam_dir"], bam)
            ctg_name = bam.replace(".bam", "")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            outputfile = os.path.join(paras["ipd_dir"], f"{ctg_name}.count")

            future = executor.submit(
                IPD_load_worker,
                fasta=fasta,
                subread_bam=bam_path,
                outputfile=outputfile,
                max_mismatch=args.max_NM,
                read_type=args.read_type,
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
    ]

    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def compare_ipd_parallel(args, paras):
    print ("Detect modified bases in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for bam in os.listdir(paras["ipd"]):
            if not bam.endswith(".ipd1.csv"):
                continue
            ctg_name = bam.replace(".ipd1.csv", "")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            control = os.path.join(paras["control"], f"{ctg_name}.ipd2.csv")
            ipd_ratio_file = os.path.join(paras["ipd_ratio"], f"{ctg_name}.ipd3.csv")
            gff = os.path.join(paras["gffs"], f"{ctg_name}.gff")
            figure_file = os.path.join(paras["figs"], f"{ctg_name}.png")

            future = executor.submit(
                get_ipd_ratio,
                csv = control,
                output = ipd_ratio_file,
                gff = gff,
                ref = fasta,
                figure_file = figure_file,
                min_cov = args.min_cov,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()

            yield finish_code

def motif_worker(bam, fasta, gff, seg_ref, seg_gff, motif, threads, min_score):
    ## run samtools depth
    depth_file = f"{bam}.depth"
    cmd = [
        "samtools", "depth", "-aa", bam
    ]

    with open(depth_file, "w") as depth_output:
        print("Running command:", " ".join(cmd))
        subprocess.run(cmd, stdout=depth_output, stderr=subprocess.DEVNULL, check=True)

    process_depth_and_gff(depth_file, fasta, gff, seg_ref, seg_gff,
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
    print("Running command:", " ".join(cmd))
    subprocess.run(cmd, check=True)

    return 0

def motif_parallel(args, paras):
    print ("Detect motif in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for gff in os.listdir(paras["gffs"]):
            if not gff.endswith(".gff"):
                continue
            if gff.endswith(".reprocess.gff"):
                continue
            # print ("####", paras["gffs"], gff)
            ctg_name = gff.replace(".gff", "")
            bam = os.path.join(paras["bam_dir"], f"{ctg_name}.bam")
            gff = os.path.join(paras["gffs"], gff)
            seg_ref = os.path.join(paras["segs"], f"{ctg_name}.seg.fa")
            seg_gff = os.path.join(paras["segs"], f"{ctg_name}.seg.gff")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            motifs = os.path.join(paras["motifs"], f"{ctg_name}.motif.csv")

            future = executor.submit(
                motif_worker,
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
            finish_code = future.result()

            yield finish_code

def collect_motifs_worker(args, paras):
    collect_motifs(
        folder=paras["motifs"],
        all_motif=paras["all_motifs"],
        MIN_FRAC=args.min_frac,
        MIN_detect=args.min_sites,)
    print ("Motif collection done.")

def profile_parallel(args, paras):
    print ("Profile motif in parallel...")
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []

        for gff in os.listdir(paras["gffs"]):
            if not gff.endswith(".gff"):
                continue
            if gff.endswith(".reprocess.gff"):
                continue
            # print ("####", paras["gffs"], gff)
            ctg_name = gff.replace(".gff", "")
            bam = os.path.join(paras["bam_dir"], f"{ctg_name}.bam")
            gff = os.path.join(paras["gffs"], gff)
            seg_ref = os.path.join(paras["segs"], f"{ctg_name}.seg.fa")
            seg_gff = os.path.join(paras["segs"], f"{ctg_name}.seg.gff")
            fasta = os.path.join(paras["ctg_dir"], f"{ctg_name}.fa")
            motifs = os.path.join(paras["motifs"], f"{ctg_name}.motif.csv")
            profile = os.path.join(paras["profiles"], f"{ctg_name}.motifs.profile.csv")
            ipd_ratio = os.path.join(paras["ipd_ratio"], f"{ctg_name}.ipd3.csv")

            future = executor.submit(
                motif_profile_worker,
                my_ref = fasta,
                gff = gff,
                all_motifs = paras["all_motifs"],
                profile = profile,
                ipd_ratio_file = ipd_ratio,
                min_frac = args.min_frac,
                min_sites = args.min_sites,
                score_cutoff = args.min_score,
                min_cov = args.min_cov,
            )
            futures.append(future)


        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()

            yield finish_code

def run_merge_profile(args, paras):
        merge_profile_worker(
            work_dir = args.work_dir, 
            heat_map = paras["heatmap"],
            profile_list = paras["profiles"], 
            total_profile = paras["total_profile"], 
            min_frac = args.min_frac, 
            whole_ref = args.whole_ref, 
            plasmid_file = args.plasmid_file,
        )

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

    paras["kmer_bin"] = os.path.join(sys.path[0], "src", "test")
    paras["ctg_list_file"] = os.path.join(args.work_dir, "contigs_list.txt")
    paras["all_motifs"] = os.path.join(args.work_dir, "all.motifs.csv")
    paras["total_profile"] = os.path.join(args.work_dir, "motif_profile.csv")
    paras["heatmap"] = os.path.join(args.work_dir, "motif_heatmap.pdf")

    ## build the directory if not exist
    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(paras["bam_dir"], exist_ok=True)
    os.makedirs(paras["ctg_dir"], exist_ok=True)
    os.makedirs(paras["ipd_dir"], exist_ok=True)
    os.makedirs(paras["control"], exist_ok=True)
    os.makedirs(paras["ipd_ratio"], exist_ok=True)
    os.makedirs(paras["gffs"], exist_ok=True)
    os.makedirs(paras["figs"], exist_ok=True)
    os.makedirs(paras["motifs"], exist_ok=True)
    os.makedirs(paras["segs"], exist_ok=True)
    os.makedirs(paras["bins"], exist_ok=True)
    os.makedirs(paras["profiles"], exist_ok=True)

    return paras

def main():
    args = parse_arguments()

    if args.max_NM is None:
        if args.read_type == "hifi":
            args.max_NM = 200
        else:
            args.max_NM = 1000

    print("🔬 Running MGE-host linkage pipeline with the following parameters:")
    for k, v in vars(args).items():
        print(f"  {k}: {v}")

    paras = get_paras(args)

    # === Insert your pipeline logic below ===
    # print("\n[Placeholder] Pipeline execution starts here...\n")
    # split_bam(args.whole_bam, args.work_dir, args.whole_ref, args.threads, args.min_len, args.max_NM)
    # print ("Splitting BAM files done.")

    for result in load_ipd_parallel(args, paras):
        print(f"IPD loading finished with code: {result}")

    get_control_parallele(args, paras)
    print ("Control file generation done.")

    for result in compare_ipd_parallel(args, paras):
        print(f"IPD ratio calculation finished with code: {result}")
    print ("IPD ratio calculation done.")

    for result in motif_parallel(args, paras):
        print(f"Motif finished with code: {result}")
    print ("Motif done.")

    collect_motifs_worker(args, paras)

    for result in profile_parallel(args, paras):
        print(f" collection finished with code: {result}")
    print ("Motif profile done.")

    run_merge_profile(args, paras)
    print ("Motif profile merging done.")
    print("\n[Placeholder] Pipeline execution ends here...\n")
    print("🔬 Pipeline execution completed.")

if __name__ == "__main__":
    main()
