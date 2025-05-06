from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import sys
from collections import defaultdict
import os

import matplotlib
matplotlib.use('Agg')  # MUST come before importing pyplot

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator



def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = str(record.seq)
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos, modified_loci, min_cov, ipd_info_dict):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1

    valid_loci_num = 0   ## number of motif sites with >=min_cov 
    motif_loci_num = 0
    motif_modify_num = 0
    for_loci_num = 0
    rev_loci_num = 0
    for_modified_num = 0
    rev_modified_num = 0

    for r, contig in REF.items():
        # contig = str(contig)
        for site in nt_search(contig, motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            # tag = r + ":" + str(site+exact_pos) + "+"
            tag = f"{r}:{site + exact_pos}+"
            if tag in ipd_info_dict:
                if ipd_info_dict[tag]['coverage'] >= min_cov:
                    valid_loci_num += 1
            motif_loci_num += 1
            for_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                for_modified_num += 1
                modified_loci[tag].add(motif_new)

        for site in nt_search(contig, Seq(motif_new).reverse_complement())[
            1:
        ]:
            # tag = r + ":" + str(site+rev_exact_pos) + "-"
            tag = f"{r}:{site + rev_exact_pos}-"
            if tag in ipd_info_dict:
                if ipd_info_dict[tag]['coverage'] >= min_cov:
                    valid_loci_num += 1
            motif_loci_num += 1
            rev_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                rev_modified_num += 1
                modified_loci[tag].add(motif_new)
    # print ("for_loci_num", for_loci_num, for_modified_num, "forward modified ratio", for_modified_num/for_loci_num)
    # print ("rev_loci_num", rev_loci_num, rev_modified_num, "reverse modified ratio", rev_modified_num/rev_loci_num)

    # print ("motif_loci_num", motif_loci_num)
    # print ("motif_modify_num", motif_modify_num)
    if motif_loci_num == 0:
        ratio = 0
    else:
        ratio = motif_modify_num/motif_loci_num
    if valid_loci_num == 0:
        valid_ratio = 0
    else:
        valid_ratio = motif_modify_num/valid_loci_num
    if for_loci_num == 0:
        for_ratio = 0
    else:
        for_ratio = for_modified_num/for_loci_num
    if rev_loci_num == 0:
        rev_ratio = 0
    else:
        rev_ratio = rev_modified_num/rev_loci_num
    if len(modified_loci) == 0:
        proportion_all_modified = 0
    else:
        proportion_all_modified = motif_modify_num/len(modified_loci)

    if ratio > 0.4 and proportion_all_modified > 0.1:
        print (motif_new, "modified_ratio", ratio, motif_modify_num, motif_loci_num, proportion_all_modified)

    return [for_loci_num, for_modified_num,for_ratio,\
            rev_loci_num, rev_modified_num, rev_ratio,\
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified, valid_loci_num, valid_ratio], modified_loci

def get_modified_ratio(gff, score_cutoff):
    ## read the gff file
    f = open(gff, "r")
    modified_loci = {}
    all_modified_loci = {}
    for line in f:
        if line[0] == "#":
            continue
        line = line.strip().split("\t")

        ref = line[0]
        pos = int(line[3]) 
        strand = line[6]
        score = int(line[5])
        all_modified_loci[ref + ":" + str(pos) + strand] = set()
        if score <= score_cutoff:
            continue
        modified_loci[ref + ":" + str(pos) + strand] = score
    print ("no. of modified loci", len(modified_loci))
    return modified_loci, all_modified_loci

def get_reprocess_gff(gff, all_modified_loci, anno):
    reprocess_gff = gff[:-4] + ".reprocess.gff"
    ## read the gff file
    out = open(reprocess_gff, "w")
    f = open(gff, "r")
    i = 0
    for line in f:
        line= line.strip()
        if line[0] == "#":
            print (line, file=out)
            continue
        else:
            if i == 0:
                print (anno, file=out, end = "")
        
        field = line.split("\t")

        ref = field[0]
        pos = int(field[3]) 
        strand = field[6]
        score = int(field[5])
        # if score <= score_cutoff:
        #     continue
        tag = ref + ":" + str(pos) + strand
        
        if tag in all_modified_loci and len(all_modified_loci[tag]) > 0:
            motifs = ",".join(list(all_modified_loci[tag]))
            print (line + ";motif=" + motifs, file=out)
        else:
            print (line, file=out)
        i += 1
    out.close()

def count_motifs(modified_loci, all_modified_loci, score_cutoff):
    non_motif_loci = 0
    motig_modified_loci = defaultdict(int)
    for tag in modified_loci:
        motifs = all_modified_loci[tag]
        if len(motifs) == 0:
            non_motif_loci += 1
        else:
            for motif in motifs:
                motig_modified_loci[motif] += 1
    if len(modified_loci) > 0:
        non_motif_loci_ratio = non_motif_loci/len(modified_loci)
    else:
        non_motif_loci_ratio = 0
    anno = ''
    anno += f"##no. of modified loci: {len(all_modified_loci)}\n"
    anno += f"##no. of filtered modified loci with score >{score_cutoff}: {len(modified_loci)} \n"
    anno += f"##no. of filtered modified loci without motifs: {non_motif_loci} \n"
    anno += f"##ratio of filtered modified loci without motifs: {non_motif_loci_ratio} \n"
    for motif, num in motig_modified_loci.items():
        if num > 50:
            anno += f"##no. of modified loci with motif {motif}: {num} \n"
    return anno

def reload_motif_sites(REF, df):
    motif_sites = {}
    for index, motif in df.iterrows():
        motif_new = motif["motifString"]
        exact_pos = motif["centerPos"]

        motif_len = len(motif_new)
        rev_exact_pos = motif_len - exact_pos + 1
        

        for r, contig in REF.items():
            for site in nt_search(str(contig), motif_new)[1:]:
                # for i in range(site, site + motif_len):
                #     motif_sites[r + ":" + str(i) + "+"] = motif_new
                tag = r + ":" + str(site+exact_pos) + "+"
                motif_sites[tag] = motif_new

            for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
                1:
            ]:
                tag = r + ":" + str(site+rev_exact_pos) + "-"
                motif_sites[tag] = motif_new
    return motif_sites

def read_ipd_ratio(ipd_ratio_file):
    ipd_info_dict = {}
    ## if file is empty, return empty dict
    if os.stat(ipd_ratio_file).st_size == 0:
        print ("ipd ratio file is empty")
        return ipd_info_dict
    else:

        # df_ipd = pd.read_csv(ipd_ratio_file)

        dtype_map = {
            "refName": "category",
            "strand": "int8",
            "tpl": "int32",
            "base": "category",
            "coverage": "int16",
            "tMean": "float32",
            "tErr": "float32",
            "control": "float32",
            "ipd_ratio": "float32",
            "kmer_count": "int32",
            "pvalue": "float32",
            "score": "int8"
        }

        df_ipd = pd.read_csv(ipd_ratio_file, dtype=dtype_map)

        df_ipd['strand'] = df_ipd['strand'].map({1: '-', 0: '+'})
        df_ipd['tag'] = df_ipd['refName'].astype(str) + ":" + (df_ipd['tpl'] + 1).astype(str) + df_ipd['strand'].astype(str)
        ipd_info_dict = df_ipd.set_index('tag').T.to_dict()


        ## convert strand to - and +
        # df_ipd['strand'] = df_ipd['strand'].replace({1: '-', 0: '+'})
        # for index, row in df_ipd.iterrows():
        #     tag = row['refName'] + ":" + str(row['tpl']+1) + row['strand']
        #     ipd_info_dict [tag] = row
    return ipd_info_dict

def count_ipd_ratio(ipd_info_dict, motif_sites, ipd_ratio_file):
    data = []
    # df_ipd = pd.read_csv(ipd_ratio_file)
    # for index, row in df_ipd.iterrows():
    #     if row['strand'] == 1:
    #         strand_string = "-"
    #     else:
    #         strand_string = "+"
    #     tag = row['refName'] + ":" + str(row['tpl']+1) + strand_string
    for tag in ipd_info_dict:
        if tag in motif_sites:
            row = ipd_info_dict[tag]
            data.append([row['refName'], str(row['tpl']+1) , row['strand'], motif_sites[tag], row['ipd_ratio']])
    site_df = pd.DataFrame(data, columns = ["refName", "tpl", "strand", "motif", "ipd_ratio"])
    # print (site_df)
    ### plot the site df using subplot using seaborn
    if len(site_df) > 0:
        ## tpl is numeric, convert it to int
        site_df['tpl'] = site_df['tpl'].astype(int)
        # Apply rolling mean smoothing, smooth it for each strand and each motif separately
        site_df['ipd_ratio'] = site_df.groupby(['refName', 'strand', 'motif'])['ipd_ratio'].transform(lambda x: x.rolling(window=100, min_periods=1).mean())

        try:
            fig, ax = plt.subplots(2, 1, figsize=(15, 7))
            ## use grid
            sns.set(style="whitegrid")

            # Get all unique motifs
            motifs = site_df['motif'].unique()
            palette = sns.color_palette("tab10", n_colors=len(motifs))  # or any other palette
            motif_palette = dict(zip(motifs, palette))


            # Plot for strand "+"
            sns.lineplot(data=site_df[site_df['strand'] == "+"], x="tpl", y="ipd_ratio", hue="motif", palette=motif_palette, ax=ax[0])
            ax[0].set_title("+")
            ## remove x label
            ax[0].set_xlabel("")
            ax[0].set_ylabel("IPD Ratio")

            # Plot for strand "-"
            sns.lineplot(data=site_df[site_df['strand'] == "-"], x="tpl", y="ipd_ratio", hue="motif", palette=motif_palette, ax=ax[1])
            ax[1].set_title("-")
            ax[1].set_xlabel("")
            ax[1].set_ylabel("IPD Ratio")

            if ax[0].legend_ is not None:
                ax[0].legend_.remove()
            if ax[1].legend_ is not None:
                ax[1].legend_.remove()

            handles, labels = ax[0].get_legend_handles_labels()

            # Create shared legend below the plots
            fig.legend(
                handles, labels,
                loc='lower center',
                bbox_to_anchor=(0.5, -0.05),
                ncol=5,
                frameon=False
            )

            # Set fixed number of X-Ticks to 5 while still hiding them
            ax[0].xaxis.set_major_locator(MaxNLocator(nbins=5)) 
            ax[1].xaxis.set_major_locator(MaxNLocator(nbins=5))

            ax[0].grid(True)
            ax[1].grid(True)
            fig = ipd_ratio_file[:-9] + ".pdf"
            plt.savefig(fig, bbox_inches='tight')
            plt.close()
        except Exception as e:
            print(f"Plotting failed: {e}")
        
    else:
        print ("no motif sites in the ipd file")

def motif_profile_worker(my_ref, gff, all_motifs, profile, ipd_ratio_file, min_frac, min_sites, score_cutoff, min_cov, misassembly = True):
    try:
        motifs = pd.read_csv(all_motifs)
        print (f"No. of raw motifs {len(motifs)}", flush=True)
        if len(motifs) == 0:
            print ("no motifs found", flush=True)
            ## get a empty profile file
            df = pd.DataFrame([], columns = ["motifString", "centerPos", "for_loci_num", "for_modified_num", "for_modified_ratio",\
                                                "rev_loci_num", "rev_modified_num", "rev_modified_ratio",\
                                                "motif_loci_num", "motif_modified_num", "all_modified_ratio", "proportion", "valid_loci_num", "motif_modified_ratio"])
            df.to_csv(profile, index=False)

            ## stop the program
            # sys.exit(0)
            return 0
        REF = read_ref(my_ref)
        # print (REF)
        modified_loci, all_modified_loci = get_modified_ratio(gff, score_cutoff)
        ipd_info_dict = read_ipd_ratio(ipd_ratio_file)
        print ("ipd ratio is loaded", flush=True)

        data = []
        for index, motif in motifs.iterrows():
            # if index % 100 == 0:
            #     print (index, motif)
            motif_new = motif["motifString"]
            exact_pos = motif["centerPos"]
            motif_profile, all_modified_loci = get_motif_sites(REF, motif_new, exact_pos, all_modified_loci, min_cov, ipd_info_dict)
            data.append([motif_new, exact_pos] + motif_profile)

        df = pd.DataFrame(data, columns = ["motifString", "centerPos", "for_loci_num", "for_modified_num", "for_modified_ratio",\
                                            "rev_loci_num", "rev_modified_num", "rev_modified_ratio",\
                                            "motif_loci_num", "motif_modified_num", "all_modified_ratio", "proportion", "valid_loci_num", "motif_modified_ratio"])
        df = df.round(4)
        df.to_csv(profile, index=False)
        anno = count_motifs(modified_loci, all_modified_loci, score_cutoff)
        get_reprocess_gff(gff, all_modified_loci, anno)


        ## filter df, keep the motifs with motif_modified_num > 1000
        df = df[df["motif_modified_num"] >= min_sites]
        df = df[df["motif_modified_ratio"] >= min_frac]
        print (len(df), "motifs left after filtering", flush=True)
        if len(df) > 0:
            if misassembly: ## detect misassembly    
                motif_sites = reload_motif_sites(REF, df)
                count_ipd_ratio(ipd_info_dict, motif_sites, ipd_ratio_file)
        else:
            print ("no motif sites left after filtering", flush=True)
        return 0
    except Exception as e:
        print(f"Error in motif_profile_worker: {e} {profile}", flush=True)


         
if __name__ == "__main__":
    # motif_new = "CTGCAG"
    # exact_pos = 5
    score_cutoff = 30
    min_cov = 10

    # my_ref = "/home/shuaiw/borg/all_test//contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_19121_L.fa"
    # gff = "/home/shuaiw/borg/all_test//gffs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_19121_L.gff"
    # all_motifs = "/home/shuaiw/borg/all_test/test_motifs.csv"
    # profile = "/home/shuaiw/borg/all_test/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_19121_L.csv"

    my_ref = sys.argv[1]
    gff = sys.argv[2]
    all_motifs = sys.argv[3]
    profile = sys.argv[4]
    ipd_ratio_file = sys.argv[5]
    min_frac = float(sys.argv[6])
    min_sites = int(sys.argv[7])
    score_cutoff = int(sys.argv[8])
    min_cov = int(sys.argv[9])

    motif_profile_worker(my_ref, gff, all_motifs, profile, ipd_ratio_file, min_frac, min_sites, score_cutoff, min_cov)






# python motif_profile.py /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s10/contigs/E_coli_H10407_1.fa /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s10/gffs/E_coli_H10407_1.gff /home/shuaiw/borg/bench/zymo_new_ref_NM3/all.motifs.csv /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s10/test.csv /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s10/ipd_ratio/E_coli_H10407_1.ipd3.csv 0.1 10 10 1