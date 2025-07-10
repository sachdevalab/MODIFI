import itertools
from collections import OrderedDict
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser

class Calc_gc:

    def __init__(self, seq):
        self.seq = seq

    def get_gc_count(self) -> int:
        gc_count = self.seq.count("G") + self.seq.count("C")

        return gc_count

    def get_gc(self) -> float:
        gc = (self.seq.count("G") + self.seq.count("C")) / self.get_seq_len()

        return gc

    def get_seq_len(self) -> int:
        seq_len = len(self.seq)

        return seq_len

class Calc_kmer_freq:
        
    def get_kmer_count(self, seq, kmer_len=4):

        kmer_combos = itertools.product(*(["ATCG"] * kmer_len))
        kmer_combos = sorted(kmer_combos)

        kmer_to_count = OrderedDict()
        kmer_to_count_norm = OrderedDict()

        for kmer in kmer_combos:
            kmer = "".join(kmer)

            kmer_count = seq.count(kmer)
            kmer_count_rc = reverse_complement(seq).count(kmer)

            kmer_count_sum = kmer_count + kmer_count_rc

            kmer_to_count[kmer] = kmer_count_sum

        kmer_count_total = sum(kmer_to_count.values())
        # print(self.ctg.id, kmer_count_total)

        for kmer in kmer_combos:
            kmer = "".join(kmer)

            kmer_count = kmer_to_count[kmer]

            if kmer_count_total > 0:
                kmer_count_norm = kmer_count / kmer_count_total

            elif kmer_count_total == 0:
                kmer_count_norm = 0

            kmer_to_count_norm[kmer] = kmer_count_norm

        return kmer_to_count, kmer_to_count_norm

    def get_kmer_count_bin(self, seq_list, kmer_len=4):

        kmer_combos = itertools.product(*(["ATCG"] * kmer_len))
        kmer_combos = sorted(kmer_combos)

        kmer_to_count = OrderedDict()
        kmer_to_count_norm = OrderedDict()

        for kmer in kmer_combos:
            kmer = "".join(kmer)
            kmer_count = 0
            kmer_count_rc = 0

            for seq in seq_list: ## for each sequence in bin
                kmer_count += seq.count(kmer)
                kmer_count_rc += reverse_complement(seq).count(kmer)

            kmer_count_sum = kmer_count + kmer_count_rc

            kmer_to_count[kmer] = kmer_count_sum

        kmer_count_total = sum(kmer_to_count.values())
        # print(self.ctg.id, kmer_count_total)

        for kmer in kmer_combos:
            kmer = "".join(kmer)

            kmer_count = kmer_to_count[kmer]

            if kmer_count_total > 0:
                kmer_count_norm = kmer_count / kmer_count_total

            elif kmer_count_total == 0:
                kmer_count_norm = 0

            kmer_to_count_norm[kmer] = kmer_count_norm

        return kmer_to_count, kmer_to_count_norm

    def cosine_similarity(self, kmer_to_count_norm, kmer_to_count_norm2):
        kmer_to_count_norm = list(kmer_to_count_norm.values())
        kmer_to_count_norm2 = list(kmer_to_count_norm2.values())
        kmer_to_count_norm = np.array(kmer_to_count_norm).reshape(1, -1)
        kmer_to_count_norm2 = np.array(kmer_to_count_norm2).reshape(1, -1)
        cos_sim = cosine_similarity(np.array(kmer_to_count_norm), np.array(kmer_to_count_norm2))[0][0]
        return cos_sim

    def get_ctg_sim(self, seq1, seq2):
        kmer_to_count, kmer_to_count_norm = self.get_kmer_count(seq1)
        kmer_to_count2, kmer_to_count_norm2 = self.get_kmer_count(seq2)
        cos_sim = self.cosine_similarity(kmer_to_count_norm, kmer_to_count_norm2)
        return cos_sim
    
    def get_bin_sim(self, seq_list1, seq_list2):
        kmer_to_count, kmer_to_count_norm = self.get_kmer_count_bin(seq_list1)
        kmer_to_count2, kmer_to_count_norm2 = self.get_kmer_count_bin(seq_list2)
        cos_sim = self.cosine_similarity(kmer_to_count_norm, kmer_to_count_norm2)
        return cos_sim        

def get_seq(contig_file):
    ## contig file only has one sequence
    for header, seq in SimpleFastaParser(open(contig_file)):
        seq = seq.upper()
        return seq

def gc_content_worker(contig_name, work_dir):
    contig_file = os.path.join(work_dir, "contigs", contig_name+'.fa')

    seq = get_seq(contig_file)
    kmer_freq = Calc_gc(seq)
    gc_count = kmer_freq.get_gc_count()
    gc = kmer_freq.get_gc()
    seq_len = kmer_freq.get_seq_len()

    return gc_count, gc, seq_len


def kmer_freq_sim_worker(contig_name1, contig_name2, work_dir):
    contig_file1 = os.path.join(work_dir, "contigs", contig_name1+'.fa')
    contig_file2 = os.path.join(work_dir, "contigs", contig_name2+'.fa')

    seq1 = get_seq(contig_file1)
    seq2 = get_seq(contig_file2)

    kmer_freq = Calc_kmer_freq()
    cos_sim = kmer_freq.get_ctg_sim(seq1, seq2)

    return cos_sim

def kmer_freq_sim_bin_worker(bin_name1, bin_name2, bin_ctg_dict, work_dir):
    # print (f"Bin {bin_name1} and {bin_name2} kmer frequency similarity")
    # print (bin_ctg_dict[bin_name1], bin_ctg_dict[bin_name2])
    seq_list_1, seq_list_2 = [], []

    bin_1_len, ben_1_gc = 0, 0
    bin_2_len, ben_2_gc = 0, 0
    for contig_name1 in bin_ctg_dict[bin_name1]:
        contig_file1 = os.path.join(work_dir, "contigs", contig_name1+'.fa')
        ## check if contig file exists
        if not os.path.exists(contig_file1):
            continue
        seq1 = get_seq(contig_file1)
        seq_list_1.append(seq1)

        kmer_freq = Calc_gc(seq1)
        gc_count = kmer_freq.get_gc_count()
        seq_len = kmer_freq.get_seq_len()
        bin_1_len += seq_len
        ben_1_gc += gc_count
    bin_1_gc = ben_1_gc / bin_1_len
    
    for contig_name2 in bin_ctg_dict[bin_name2]:
        contig_file2 = os.path.join(work_dir, "contigs", contig_name2+'.fa')
        if not os.path.exists(contig_file2):
            continue
        seq2 = get_seq(contig_file2)
        seq_list_2.append(seq2)
        kmer_freq = Calc_gc(seq2)
        gc_count = kmer_freq.get_gc_count()
        seq_len = kmer_freq.get_seq_len()
        bin_2_len += seq_len
        ben_2_gc += gc_count
    bin_2_gc = ben_2_gc / bin_2_len
    print(f"Bin {bin_name1} GC content: {bin_1_gc}")
    print(f"Bin {bin_name2} GC content: {bin_2_gc}")

    kmer_freq = Calc_kmer_freq()
    cos_sim = kmer_freq.get_bin_sim(seq_list_1, seq_list_2)

    return bin_name1, round(bin_1_gc, 2), round(bin_2_gc, 2), round(cos_sim, 2)

if __name__ == "__main__":
    work_dir = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30_rec4/"

    # contig_name = "E_coli_H10407_1"
    # kmer_len = 4
    # gc_count, gc, seq_len = gc_content_worker(contig_name, work_dir)
    # print(f"GC count: {gc_count}")
    # print(f"GC content: {gc}")
    # print(f"Sequence length: {seq_len}")

    # contig_name1 = "E_coli_H10407_1"
    # contig_name2 = "E_coli_H10407_2"
    # cos_sim = kmer_freq_sim_worker(contig_name1, contig_name2, work_dir)
    # print(f"Cosine similarity: {cos_sim}")

    bin_ctg_dict = {"1":["E_coli_H10407_1", "E_coli_H10407_2"], "2":["B_cepacia_UCB-717_1"]}
    cos_sim = kmer_freq_sim_bin_worker("1", "2", bin_ctg_dict, work_dir)
    print(f"Cosine similarity: {cos_sim}")