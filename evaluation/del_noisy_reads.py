import pysam
import numpy as np
from tqdm import tqdm

class NoisyReadFilter:
    def __init__(self, input_bam, output_bam, max_nm, percentile):
        self.input_bam = input_bam
        self.output_bam = output_bam
        self.max_nm = max_nm
        self.percentile = percentile

    def compute_avg_ipd_per_read(self):
        """
        Reads a BAM file and computes the average IPD per read, filtering by NM tag.
        """
        reads_with_ipd = []
        with pysam.AlignmentFile(self.input_bam, "rb") as bam:
            for read in tqdm(bam.fetch()):
                try:
                    if read.get_tag("NM") > self.max_nm:
                        continue
                    fi = np.array(read.get_tag("fi")[::-1])  # Reverse fi
                    ri = np.array(read.get_tag("ri"))
                    full_ipd = np.concatenate([fi, ri])
                    avg_ipd = np.mean(full_ipd)
                    reads_with_ipd.append((read.query_name, avg_ipd))
                except KeyError:
                    continue  # Skip reads without fi/ri tags
        return reads_with_ipd

    def determine_ipd_cutoff(self, reads_with_ipd):
        """
        Determines the IPD cutoff based on the given percentile.
        """
        all_ipds = [r[1] for r in reads_with_ipd]
        return np.percentile(all_ipds, self.percentile)

    def filter_reads(self, keep_reads):
        """
        Writes filtered reads to a new BAM file.
        """
        with pysam.AlignmentFile(self.input_bam, "rb") as bam_in, \
             pysam.AlignmentFile(self.output_bam, "wb", template=bam_in) as bam_out:
            for read in bam_in.fetch():
                if read.query_name in keep_reads:
                    bam_out.write(read)
        pysam.index(self.output_bam)

    def run(self):
        """
        Main function to filter noisy reads from a BAM file.
        """
        # Step 1: Compute average IPD per read
        reads_with_ipd = self.compute_avg_ipd_per_read()

        # Step 2: Determine IPD cutoff
        cutoff = self.determine_ipd_cutoff(reads_with_ipd)
        print(f"IPD cutoff ({self.percentile}th percentile): {cutoff:.2f}")

        # Step 3: Make a set of read names to keep
        keep_reads = {qname for qname, avg in reads_with_ipd if avg <= cutoff}

        # Step 4: Write filtered reads to new BAM
        self.filter_reads(keep_reads)
        print(f"Filtered BAM written to {self.output_bam}")


# Parameters
input_bam = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/bams/E_coli_H10407_1.bam"
output_bam = "/home/shuaiw/borg/bench/test/filtered.bam"
MAX_NM = 500
PERCENTILE = 95

if __name__ == "__main__":
    filter_instance = NoisyReadFilter(input_bam, output_bam, MAX_NM, PERCENTILE)
    filter_instance.run()
