import re
from Bio.Seq import Seq  # BioPython is used for reverse complement functionality

from Bio.Data import IUPACData

# Precompute IUPAC base mappings
IUPAC_BASES = {k: set(v) for k, v in IUPACData.ambiguous_dna_values.items()}

def iupac_compatible(b1, b2):
    """Return True if IUPAC bases b1 and b2 are compatible (non-empty intersection)."""
    s1 = IUPAC_BASES.get(b1.upper(), set(b1.upper()))
    s2 = IUPAC_BASES.get(b2.upper(), set(b2.upper()))
    return bool(s1 & s2)

def center_align_iupac_similarity(motif1, cent1, motif2, cent2):
    """
    Center-align two motifs and compute IUPAC compatibility score.
    
    Returns:
        aligned1: Padded motif1
        aligned2: Padded motif2
        match_pct: % of positions that are compatible
    """
    ## count number of N in two motifs
    n1 = motif1.count('N')
    n2 = motif2.count('N')
    if len(motif1) != len(motif2) or n1 > 2 or n2 > 2:
        return None, None, 0.0
    # Pad left and right to align centers
    left_pad = max(cent2 - cent1, 0)
    right_pad = max(cent1 - cent2, 0)
    
    padded1 = 'N' * left_pad + motif1 + 'N' * right_pad
    padded2 = 'N' * right_pad + motif2 + 'N' * left_pad

    # Equalize length by padding to the right
    max_len = max(len(padded1), len(padded2))
    padded1 = padded1.ljust(max_len, 'N')
    padded2 = padded2.ljust(max_len, 'N')

    # Compare base by base
    compatible_count = 0
    for b1, b2 in zip(padded1, padded2):
        if iupac_compatible(b1, b2):
            compatible_count += 1

    match_pct = 100.0 * compatible_count / max_len

    return padded1, padded2, match_pct


class MotifFilter:
    IUPAC_CODES = {
        'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
        'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
        'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
    }

    def __init__(self, motif_list, key="motif", score_key="host_total", value_key="host_meth"):
        """
        Initialize the MotifFilter class.

        Parameters:
            motif_list (list of dict): List of motif dictionaries.
            key (str): The dictionary key that contains the motif string.
            score_key (str): The key used to prioritize motifs when redundant.
            value_key (str): The key used to prioritize motifs during deduplication.
        """
        self.motif_list = motif_list
        self.key = key
        self.score_key = score_key
        self.value_key = value_key

    def filter(self):
        """
        Removes redundant motifs, keeping the one with the higher host_total if motifs overlap,
        are reverse complements, or one is a subset of the other (including degenerate cases).

        Returns:
            list of dict: Filtered list of non-redundant motifs.
        """
        filtered = []

        for motif in self.motif_list:
            mseq = motif[self.key]
            keep = True
            to_remove = []

            for existing in filtered:
                eseq = existing[self.key]

                # Check if motifs are related (one contains the other, reverse complement, or degenerate subset)
                if self.is_subset_or_reverse_complement(mseq, eseq):
                    # Compare host_total scores
                    if motif[self.score_key] > existing[self.score_key]:
                        to_remove.append(existing)  # Mark existing for removal
                    else:
                        keep = False  # Don't keep the current motif
                        break

            # Remove lower-scoring overlapping motifs
            for rem in to_remove:
                filtered.remove(rem)

            if keep:
                filtered.append(motif)

        # Deduplicate by core and host meth
        filtered = self.deduplicate_by_core_and_host_meth(filtered)
        return filtered

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

## for motif list, remove redundant motifs, keeping the one with the higher host_total if motifs overlap,
def uniq_similar_motifs(motif_list):
    ## using center_align_iupac_similarity
    filtered = []
    for motif in motif_list:
        mseq = motif['motif']
        keep = True
        to_remove = []

        for existing in filtered:
            eseq = existing['motif']
            a1, a2, sim = center_align_iupac_similarity(mseq, motif['centerPos'], eseq, existing['centerPos'])
            if sim > 0.8:
                # Compare host_total scores
                if motif['host_meth'] > existing['host_meth']:
                    to_remove.append(existing)
                else:
                    keep = False
                    break
        # Remove lower-scoring overlapping motifs
        for rem in to_remove:
            filtered.remove(rem)
        if keep:
            filtered.append(motif)
    # print ("xx", motif_list, filtered, len(motif_list), len(filtered))
    return filtered


if __name__ == "__main__":
    # Example usage
    # motif_list = [
    #     {'motif': 'GSCAGCNNNNNBNNV', 'centerPos': 4, 'host_total': 7114, 'host_meth': 3927, 'plasmid_total': 10, 'plasmid_meth': 6},
    #     {'motif': 'GAAGAG', 'centerPos': 5, 'host_total': 1250, 'host_meth': 932, 'plasmid_total': 27, 'plasmid_meth': 1},
    #     {'motif': 'HNCTCTTCNNNVNNNNB', 'centerPos': 3, 'host_total': 823, 'host_meth': 554, 'plasmid_total': 15, 'plasmid_meth': 0},
    #     {'motif': 'CGBCAG', 'centerPos': 5, 'host_total': 7352, 'host_meth': 4933, 'plasmid_total': 30, 'plasmid_meth': 14},
    #     {'motif': 'GGCAGC', 'centerPos': 4, 'host_total': 3057, 'host_meth': 2617, 'plasmid_total': 16, 'plasmid_meth': 16},
    #     {'motif': 'CGCCAG', 'centerPos': 5, 'host_total': 3803, 'host_meth': 3313, 'plasmid_total': 3, 'plasmid_meth': 0}
    # ]

    # motif_list = [{'motif': 'GAAGAGNNNNNNNNNS', 'centerPos': 5, 'host_total': 969, 'host_meth': 724, 'plasmid_total': 15, 'plasmid_meth': 0}, {'motif': 'CGBCAG', 'centerPos': 5, 'host_total': 7352, 'host_meth': 4933, 'plasmid_total': 30, 'plasmid_meth': 14}, {'motif': 'GGCAGCNNNNNNNNV', 'centerPos': 4, 'host_total': 2745, 'host_meth': 2359, 'plasmid_total': 11, 'plasmid_meth': 11}, {'motif': 'VNNCGCCAG', 'centerPos': 8, 'host_total': 3359, 'host_meth': 2929, 'plasmid_total': 2, 'plasmid_meth': 0}, {'motif': 'GSCAGCNNNNNBNNV', 'centerPos': 4, 'host_total': 7114, 'host_meth': 3927, 'plasmid_total': 10, 'plasmid_meth': 6}]

    # motif_filter = MotifFilter(motif_list)
    # filtered_motifs = motif_filter.filter()
    # print(filtered_motifs)
    # print (len(filtered_motifs), "motifs after filtering.")

    # seq1 = "GSCAGCNNNNNBNNV"
    # seq2 = "GGCAGCNNNNNNNNV"
    # print (motif_filter.is_subset_or_reverse_complement(seq1, seq2))
    # print (motif_filter.is_subset_or_reverse_complement(seq2, seq1))

    motif1 = "RGATCY"
    cent1 = 2  # A in position 2

    motif2 = "GCGATC"
    cent2 = 3  # A in position 3

    a1, a2, sim = center_align_iupac_similarity(motif1, cent1, motif2, cent2)

    print("Motif 1 aligned:", a1)
    print("Motif 2 aligned:", a2)
    print(f"Compatibility: {sim:.1f}%")
