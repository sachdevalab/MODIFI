import subprocess
import os
import sys
from Bio import SeqIO
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_unique_motifs, My_sample, Isolation_sample

def cluster_mge_95ani(input_mge, output_dir, threads=8, ani_threshold=95, min_tcov=0.85):
    """
    Cluster MGE sequences based on 95% ANI similarity.
    
    Args:
        input_mge: Path to input MGE FASTA file
        output_dir: Directory for output files
        threads: Number of threads to use (default: 8)
        ani_threshold: ANI threshold for clustering (default: 95)
        min_tcov: Minimum target coverage (default: 0.85)
    
    Returns:
        dict: Paths to output files
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output file paths
    blast_out = os.path.join(output_dir, "megablast.95ani.tsv")
    ani_out = os.path.join(output_dir, "megablast.ani.95ani.tsv")
    cluster_out = os.path.join(output_dir, "megablast.cluster.95ani.tsv")
    blastdb = os.path.join(output_dir, "ani95_blast_db")
    
    try:
        # Create BLAST database
        print("Creating BLAST database...")
        subprocess.run([
            "makeblastdb", "-in", input_mge, "-dbtype", "nucl", "-out", blastdb
        ], check=True)
        
        # Run BLAST
        print("Running BLAST...")
        subprocess.run([
            "blastn", "-query", input_mge, "-db", blastdb, 
            "-outfmt", "6 std qlen slen", "-max_target_seqs", "10000",
            "-out", blast_out, "-num_threads", str(threads)
        ], check=True)
        
        # Clean up BLAST database files
        for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto"]:
            db_file = blastdb + ext
            if os.path.exists(db_file):
                os.remove(db_file)
        
        # Calculate ANI
        print("Calculating ANI...")
        subprocess.run([
            "python", "/groups/sachdeva/projects/sag/SAGLink/workflow/anicalc.py",
            "-i", blast_out, "-o", ani_out
        ], check=True)
        
        # Cluster based on ANI
        print("Clustering based on ANI...")
        subprocess.run([
            "python", "/groups/sachdeva/projects/sag/SAGLink/workflow/aniclust.py",
            "--fna", input_mge, "--ani", ani_out, "--out", cluster_out,
            "--min_ani", str(ani_threshold), "--min_tcov", str(min_tcov)
        ], check=True)
        
        print("Clustering completed successfully!")
        
        return {
            "blast_out": blast_out,
            "ani": ani_out,
            "cluster": cluster_out
        }
        
    except subprocess.CalledProcessError as e:
        print(f"Error during clustering: {e}")
        raise

def main(all_dir, fig_dir, input_mge):
    i = 0
    ## using biopython, output all MGE sequences to input_mge
    with open(input_mge, "w") as f:
        for prefix in os.listdir(all_dir):
            sample_obj = My_sample(prefix, all_dir)
            sample_obj.read_MGE()
            # check if sample_obj.mge_dict is empty
            if not sample_obj.mge_dict:
                print (f">>>No MGE found in sample {prefix}")
                continue
            for mge in sample_obj.mge_dict:
                ## using biopython to read sample_obj.reference_fasta
                ## retrieve the record with id mge
                # Create a dictionary to store all records from the reference fasta
                cmd = "samtools faidx " + sample_obj.reference_fasta + " " + mge + " >> " + input_mge
                os.system(cmd)
                i += 1
    print (f">>>Total {i} MGE sequences written to {input_mge}")





if __name__ == "__main__":
    # Example usage
    input_mge = "/home/shuaiw/borg/paper/network/all_mge.fa"
    output_dir = "/home/shuaiw/borg/paper/MGE/cluster"
    all_dir = "/home/shuaiw/borg/paper/run2/"
    fig_dir = "../../tmp/figures/multi_env_linkage/"
    main(all_dir, fig_dir, input_mge)
    # results = cluster_mge_95ani(input_mge, output_dir, threads=32)
    # print("Output files:", results)