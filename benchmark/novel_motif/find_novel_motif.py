import re
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

class MotifCompare:
    IUPAC_CODES = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
        'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
    }

    # def __init__(self):
    #     x = 1

    def is_subset_or_reverse_complement(self, seq1, seq2):
        """
        Check if one sequence is a subset of the other, considering degenerate bases,
        or if they are reverse complements.

        Parameters:
            seq1 (str): The first sequence.
            seq2 (str): The second sequence.

        Returns:
            bool: True if seq1 is a subset of seq2, seq2 is a subset of seq1, or they are reverse complements.
        """
        # Check if seq1 is a subset of seq2
        if self.is_degenerate_subset(seq1, seq2):
            return True

        # Check if seq2 is a subset of seq1
        if self.is_degenerate_subset(seq2, seq1):
            return True

        # Check reverse complement
        reverse_seq1 = str(Seq(seq1).reverse_complement())
        if self.is_degenerate_subset(reverse_seq1, seq2):
            return True

        if self.is_fully_iupac_compatible(seq1, seq2):
            return True

        return False

    def is_degenerate_subset(self, seq1, seq2):
        """
        Check if seq1 is a subset of seq2, considering degenerate bases.

        Parameters:
            seq1 (str): The first sequence.
            seq2 (str): The second sequence.

        Returns:
            bool: True if seq1 is a subset of seq2, False otherwise.
        """
        for base1, base2 in zip(seq1, seq2):
            if not set(self.IUPAC_CODES[base1]).issubset(self.IUPAC_CODES[base2]):
                return False
        return True

    def deduplicate_by_core_and_host_meth(self, motif_list):
        """
        Remove redundant motifs by collapsing ones with the same core motif or reverse complement,
        keeping only the one with the highest host_meth value.

        Parameters:
            motif_list (list of dict): List of motif dictionaries.

        Returns:
            list of dict: Deduplicated list of motifs.
        """
        core_map = {}
        for motif in motif_list:
            core = self.extract_core_simple(motif[self.key])
            reverse_core = self.extract_core_simple(str(Seq(core).reverse_complement()))
            if core not in core_map and reverse_core not in core_map:
                core_map[core] = motif
            elif core in core_map:
                if motif[self.value_key] > core_map[core][self.value_key]:
                    core_map[core] = motif
            elif reverse_core in core_map:
                if motif[self.value_key] > core_map[reverse_core][self.value_key]:
                    core_map[reverse_core] = motif
        return list(core_map.values())

    def extract_core_simple(self, motif):
        """
        Simplify a motif by collapsing IUPAC codes and extracting the longest non-N segment.

        Parameters:
            motif (str): The motif string.

        Returns:
            str: The core motif.
        """
        # Replace degenerate bases with 'N'
        simplified = re.sub(r'[BDHVNRYSWKM]', 'N', motif)
        # Collapse consecutive 'N's
        collapsed = re.sub(r'N+', 'N', simplified)
        # Split into segments and return the longest as core
        segments = re.split(r'N+', collapsed)
        return max(segments, key=len) if segments else motif

    def is_fully_iupac_compatible(self, seq1, seq2):
        """
        Check if two IUPAC motifs are fully compatible base by base.

        Parameters:
            seq1 (str): First motif
            seq2 (str): Second motif

        Returns:
            bool: True if all positions are compatible under IUPAC rules
        """
        IUPAC = {
            'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
            'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
            'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
            'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
        }

        for b1, b2 in zip(seq1, seq2):
            if not (IUPAC.get(b1, set()) & IUPAC.get(b2, set())):
                return False
        return True


def collect_rebase_motif(rebase_protein, all_rebase_motif):
    rebase_motif = set()
    with open(rebase_protein, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            match = re.search(r"RecSeq:([ACGTRYWSMKHBVDN\-x]+)\t", line)
            if match:
                motif = match.group(1)
                # print("Motif:", motif)
                rebase_motif.add(motif)
    print (f"Total unique motifs found: {len(rebase_motif)}")
    ## sort the motifs
    sorted_motifs = sorted(rebase_motif)
    with open(all_rebase_motif, "w") as out:
        for motif in sorted_motifs:
            out.write(motif + "\n")

def load_rebase_motif(all_rebase_motif):
    """
    Load the rebase motif file and return a set of motifs.
    """
    with open(all_rebase_motif, "r") as f:
        motifs = set()
        for line in f:
            motif = line.strip()
            if motif:
                motifs.add(motif)
    return motifs

def check_if_motif_in_rebase(motif, rebase_motifs):
    flag = False
    if motif in rebase_motifs:
        flag =  True
    for rebase_motif in rebase_motifs:
        motif_compare = MotifCompare()
        if len(motif) <= len(rebase_motif) and motif_compare.is_subset_or_reverse_complement(motif, rebase_motif):
            flag = True
            break
    return flag
    

def load_sample_motif(all_rebase_motif, sample_motif_file):
    """
    Load the sample motif file and return a set of motifs.
    """
    rebase_motifs = load_rebase_motif(all_rebase_motif)
    novel = 0
    with open(sample_motif_file, "r") as f:
        motifs = set()
        for line in f:
            if line.startswith("motifString"):
                continue
            motif = line.strip().split(",")[0]
            exit_flag = check_if_motif_in_rebase(motif, rebase_motifs)
            if not exit_flag:
                print (motif, line, end = '')
                novel += 1
    print (f"Total novel motifs found: {novel}")
            



rebase_protein = "/home/shuaiw/methylation/data/rebase/protein_org_seqs.txt"
all_rebase_motif = "/home/shuaiw/borg/paper/novel_motif/rebase.sort.motif.txt"
# collect_rebase_motif(rebase_protein, all_rebase_motif)  
# 
# sample_motif_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/all.motifs.csv"  
sample_motif_file = "/home/shuaiw/borg/bench/soil/run1/all.motifs.csv"
load_sample_motif(all_rebase_motif, sample_motif_file)

## xxx