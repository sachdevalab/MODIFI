import numpy as np
from math import log
import pandas as pd
import os
import argparse
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from derep_motifs import MotifFilter

# IGNORE_MOTIFS = [
#     "GATC",
#     "CCWGG",
# ]

IGNORE_MOTIFS = [
]

def invasion_score_from_counts(motif_data, min_frac=0.5, neutral_score=1.0, max_sites=50000):
    """
    Adds confidence weighting based on motif site counts.
    """
    scores = []
    weights = []
    total_sites = 0
    total_host_sites = 0

    for m in motif_data:
        h_total = m['host_total']
        h_meth = m['host_meth']
        p_total = m['plasmid_total']
        p_meth = m['plasmid_meth']

        if h_total == 0:
            continue

        weight = h_total
        weight = min(weight, 1000)   ### set weight to 1000, once the h_total is larger than 1000
        total_sites += h_total + p_total
        total_host_sites += h_total

        f_host = h_meth / h_total
        if p_total == 0:
            motif_score = neutral_score
        else:
            f_plasmid = p_meth / p_total
            # motif_score = 1 - abs(f_host - f_plasmid)
            ## whether f_host and f_plasmid both are larger than 0.5
            if f_plasmid > min_frac:
                motif_score = 1
            elif f_plasmid > f_host:
                motif_score = 1
            else:
                motif_score = 1 - abs(f_host - f_plasmid)/f_host   ## if the f_host is only 0.5, so divided by f_host to normalize the score
        # print (m['motif'], motif_score, weight)
        scores.append(motif_score * log(weight))
        weights.append(log(weight))
    if not scores:
        return {'invasion_score': 0.0, 'confidence': 0.0, 'final_score': 0.0}

    invasion_score = sum(scores) / sum(weights)

    # Confidence scaling (logarithmic)
    confidence = log(1 + total_sites) / log(1 + max_sites)
    if confidence > 1:
        confidence = 1

    # motif confidence
    motif_confidence = log(1+len(motif_data))/log(1+3)   
    if motif_confidence > 1:
        motif_confidence = 1
    final_score = invasion_score * confidence * motif_confidence

    return {
        'invasion_score': round(invasion_score, 4),
        'confidence': round(confidence, 4),
        'final_score': round(final_score, 4),
        'total_sites': total_sites,
        'motif_confidence': round(motif_confidence, 4),
        'host_motif_num': len(motif_data)
    }


def extract_motif_data(host_df, plasmid_profile, min_frac = 0.5, min_detect = 100):
    # host_df = pd.read_csv(host_profile)
    plasmid_df = pd.read_csv(plasmid_profile)
    ## make a new df, by extract the motifString, centerPos, motif_loci_num, motif_modified_num, motif_modified_ratio from the original two dfs
    motif_data = pd.merge(host_df, plasmid_df, on = ['motifString', 'centerPos'], suffixes=('_host', '_plasmid'))
    motif_data = motif_data[(motif_data['motif_modified_ratio_host'] > min_frac) & (motif_data['motif_modified_num_host'] > min_detect)]
    motif_data = motif_data[['motifString', 'centerPos', 'motif_loci_num_host', 'motif_modified_num_host', 'motif_loci_num_plasmid', 'motif_modified_num_plasmid']]
    motif_data.columns = ['motif', 'centerPos', 'host_total', 'host_meth', 'plasmid_total', 'plasmid_meth']
    # print (motif_data)
    motif_data = motif_data.to_dict('records')
    return motif_data

def merge_bin_profile(ctg_profile_dict, bin_name, bin_ctg_dict):
    """
    Merge the profile of all contigs in a bin into a single dataframe.
    """
    ctg_list = bin_ctg_dict[bin_name]
    ## remove the ctg not in ctg_profile_dict
    ctg_list = [ctg for ctg in ctg_list if ctg in ctg_profile_dict]
    if len(ctg_list) == 0:
        return pd.DataFrame()
    elif len(ctg_list) == 1:
        return ctg_profile_dict[ctg_list[0]]
    else:
        bin_df = ctg_profile_dict[ctg_list[0]]
        ## only retain motifString,centerPos,motif_loci_num, motif_modified_num, and motif_modified_ratio
        bin_df = bin_df[['motifString', 'centerPos', 'motif_loci_num', 'motif_modified_num', 'motif_modified_ratio']]
        for ctg_name in ctg_list[1:]:
            if ctg_name not in ctg_profile_dict:
                print (f"{ctg_name} has no methylation profile.")
                continue
            df = ctg_profile_dict[ctg_name]
            merge_df = pd.merge(bin_df, df, on = ['motifString', 'centerPos'], suffixes=('_x', '_y'))
            merge_df['motif_loci_num'] = merge_df['motif_loci_num_x'] + merge_df['motif_loci_num_y']
            merge_df['motif_modified_num'] = merge_df['motif_modified_num_x'] + merge_df['motif_modified_num_y']
            ## avoid divided by 0, keep original value
            merge_df['motif_modified_ratio'] = merge_df['motif_modified_num'] / merge_df['motif_loci_num'].replace(0, 1)
            bin_df = merge_df[['motifString', 'centerPos', 'motif_loci_num', 'motif_modified_num', 'motif_modified_ratio']]
        return bin_df

def merge_bin_motif_file(ctg_motif_dict, bin_name, bin_ctg_dict):
    ctg_list = bin_ctg_dict[bin_name]
    if len(ctg_list) == 1:
        if ctg_list[0] not in ctg_motif_dict:
            return pd.DataFrame()
        return ctg_motif_dict[ctg_list[0]]
    else:
        motif_df = pd.DataFrame()
        for ctg_name in ctg_list:
            if ctg_name not in ctg_motif_dict:
                # print (f"{ctg_name} has no motif file.")
                continue
            df = ctg_motif_dict[ctg_name]
            if not df.empty:
                motif_df = pd.concat([motif_df, df], axis=0)
        motif_df = motif_df.drop_duplicates()
        return motif_df

def merge_bin_motif(bin_ctg_dict, ctg_motif_dict, ctg_profile_dict):
    bin_df_dict = {}
    bin_motif_dict = {}
    for bin_name in bin_ctg_dict:
        # if bin_name != "RuReacBro_20230708_7_40h_PC_r3_metabinner_4_rmcirc":
        #     continue
        bin_df = merge_bin_profile(ctg_profile_dict, bin_name, bin_ctg_dict)
        bin_motif = merge_bin_motif_file(ctg_motif_dict, bin_name, bin_ctg_dict)
        bin_df_dict[bin_name] = bin_df
        bin_motif_dict[bin_name] = bin_motif
    return bin_df_dict, bin_motif_dict
        

def for_each_plasmid(bin_df_dict, bin_motif_dict, bin_ctg_dict, ctg_profile_dict, ctg_motif_dict,\
                      plasmid_name, profile_dir, host_dir, min_frac = 0.5, MGE_dict={}):
    plasmid_profile = f"{profile_dir}/{plasmid_name}.motifs.profile.csv"
    ## check if the plasmid profile exists
    if not os.path.exists(plasmid_profile):
        print (f"{plasmid_profile} does not exist.")
        return
    data = []
    motif_data_dict = {}
    for bin_name in bin_ctg_dict:
        # print (f"Processing {bin_name}...")
        bin_df = bin_df_dict[bin_name]
        if len(bin_df) == 0:
            continue
        bin_motif = bin_motif_dict[bin_name]

        # bin_df = merge_bin_profile(ctg_profile_dict, bin_name, bin_ctg_dict)
        # if len(bin_df) == 0:
        #     continue
        # bin_motif = merge_bin_motif_file(ctg_motif_dict, bin_name, bin_ctg_dict)
        # print (bin_motif)

        motif_data = extract_motif_data(bin_df, plasmid_profile, min_frac)
        motif_data = filter_motifs(bin_motif, motif_data)

        motif_filter = MotifFilter(motif_data)
        motif_data = motif_filter.filter()

        if bin_name == "RuReacBro_20230708_7_48h_PC_r3_vamb_18_rmcirc":
            print (motif_data)
        motif_data_dict[bin_name] = motif_data
        result = invasion_score_from_counts(motif_data, min_frac, 0.5)
        result['host'] = bin_name
        data.append(result)
    ## convert the data to a df, and sort by final_score
    data = pd.DataFrame(data)
    data = data.sort_values(by = 'final_score', ascending = False, ignore_index = True)
    ## host column the first, final_score the second, invasion_score the third, confidence the forth, total_sites the fifth
    data = data[['host', 'final_score', 'invasion_score', 'confidence', 'motif_confidence', 'total_sites', 'host_motif_num']]
    ## print top ten
    ## remove the rows with final_score = 0
    data = data[data['final_score'] > 0]
    print (data.head(10))
    ## output the data to a csv file
    host_prediction = os.path.join(host_dir, f"{plasmid_name}.host_prediction.csv")
    data.to_csv(host_prediction, index = False)
    ## output the motif data to a csv file
    motif_data_file = os.path.join(host_dir, f"{plasmid_name}.motif_data.csv")
    with open(motif_data_file, 'w') as f:
        for index, rows in data.iterrows():
            print (rows, motif_data_dict[rows['host']], file = f)

def read_genomad(genomad_file):
    MGE_dict = {}
    print (f"Reading {genomad_file}...")
    genomad = pd.read_csv(genomad_file, sep = "\t")
    for i, row in genomad.iterrows():
        
        if re.search('\|provirus', row['seq_name']):
            continue
        # print (row['seq_name'])
        MGE_dict[row['seq_name']] = row
    return MGE_dict

def filter_motifs(host_motif, motif_data):
    # host_motifs = pd.read_csv(host_motif)
    host_motifs = host_motif
    ## filter motifs in motif data, only keep the motifs in host_motifs by examing the motifString and centerPos
    motif_tag_dict = {}
    for index, row in host_motifs.iterrows():
        motif_tga = row['motifString'] + "_tga_" + str(row['centerPos'])
        motif_tag_dict[motif_tga] = 1
    new_motif_data = []

    for m in motif_data:
        ## remove the ignored motifs
        ignore_flag = False
        for ignore_motif in IGNORE_MOTIFS:
            if ignore_motif in m['motif']:
                ignore_flag = True
                break
        if ignore_flag:
            continue

        tag = m['motif'] + "_tga_" + str(m['centerPos'])
        if tag in motif_tag_dict:
            new_motif_data.append(m)
    return new_motif_data

def count_MGE_with_motif(plasmid_name, profile_dir):
    plasmid_motif_file = os.path.join(profile_dir,"../motifs", f"{plasmid_name}.motifs.csv")
    if not os.path.exists(plasmid_motif_file):
        print (f"{plasmid_motif_file} does not exist.")
        return 0
    plasmid_motif = pd.read_csv(plasmid_motif_file)
    return len(plasmid_motif)

def summary_host(host_dir):
    data = []
    for file in os.listdir(host_dir):
        if file.endswith(".host_prediction.csv"):
            plasmid_name = file.split(".")[0]
            host_prediction = os.path.join(host_dir, file)
            df = pd.read_csv(host_prediction)
            ## extract the first row
            ## add new column for plasmid_name at the start
            if len(df) > 0:
                df.insert(0, 'plasmid', plasmid_name)
                data.append(df.iloc[0])
    data = pd.DataFrame(data)
    ## sort by final_score
    data = data.sort_values(by = 'final_score', ascending = False, ignore_index = True)
    ## output the data to a csv file
    host_summary = os.path.join(host_dir, "../", "host_summary.csv")
    data.to_csv(host_summary, index = False)

def batch_MGE_invade(plasmid_file, profile_dir, host_dir, bin_file=None, min_frac = 0.5, threads = 1):

    MGE_dict = read_genomad(plasmid_file)
    ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict = load_ctg_motifs(profile_dir, MGE_dict, threads)
    if bin_file is not None:
        if not os.path.exists(bin_file):
            print (f"{bin_file} does not exist.")
        ctg_bin_dict, bin_ctg_dict = load_bin(bin_file)
    else:
        ctg_bin_dict, bin_ctg_dict = ctg_single_dict, single_ctg_dict
    print (f"Loading is done, {len(bin_ctg_dict)} bins.")
    bin_df_dict, bin_motif_dict = merge_bin_motif(bin_ctg_dict, ctg_motif_dict, ctg_profile_dict)
    print ("threads number for linkage prediction: ", threads)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for plasmid_name in MGE_dict:
            MGE_motif_num = count_MGE_with_motif(plasmid_name, profile_dir)
            if MGE_motif_num == 0:
                print (f"Skip {plasmid_name} with {MGE_motif_num} motifs.")
                continue
            print (f"Processing {plasmid_name} with {MGE_motif_num} original motifs.")
            future = executor.submit(
                for_each_plasmid,
                bin_df_dict = bin_df_dict,
                bin_motif_dict = bin_motif_dict,
                bin_ctg_dict = bin_ctg_dict,
                ctg_profile_dict = ctg_profile_dict,
                ctg_motif_dict = ctg_motif_dict,
                plasmid_name = plasmid_name,
                profile_dir = profile_dir,
                host_dir = host_dir,
                min_frac = min_frac,
                MGE_dict={},
            )
            futures.append(future)
        result = []
        for future in tqdm(as_completed(futures), total=len(futures)):
            finish_code = future.result()
            result.append(finish_code)
        
            # for_each_plasmid(bin_ctg_dict, ctg_profile_dict, ctg_motif_dict, plasmid_name, profile_dir, host_dir, min_frac, {})
    summary_host(host_dir)

def batch_MGE_invade_single(plasmid_file, profile_dir, host_dir, bin_file=None, min_frac = 0.5):

    MGE_dict = read_genomad(plasmid_file)
    ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict = load_ctg_motifs(profile_dir, MGE_dict)
    if bin_file is not None:
        if not os.path.exists(bin_file):
            print (f"{bin_file} does not exist.")
        ctg_bin_dict, bin_ctg_dict = load_bin(bin_file)
    else:
        ctg_bin_dict, bin_ctg_dict = ctg_single_dict, single_ctg_dict
    print (f"Loading is done, {len(bin_ctg_dict)} bins.")
    for plasmid_name in MGE_dict:
        MGE_motif_num = count_MGE_with_motif(plasmid_name, profile_dir)
        if MGE_motif_num == 0:
            print (f"Skip {plasmid_name} with {MGE_motif_num} motifs.")
            continue
        print (f"Processing {plasmid_name} with {MGE_motif_num} original motifs.")
        for_each_plasmid(bin_ctg_dict, ctg_profile_dict, ctg_motif_dict, plasmid_name, profile_dir, host_dir, min_frac, {})
    summary_host(host_dir)

def load_bin(bin_file):
    bin_df = pd.read_csv(bin_file, sep= "\t", header = 0)
    ## add header for the bin_df
    bin_df.columns = ['contig', 'bin_id']
    # print (bin_df)
    ctg_bin_dict = {}
    bin_ctg_dict = defaultdict(list)

    for i, row in bin_df.iterrows():
        contig = row['contig']
        bin_id = row['bin_id']
        ctg_bin_dict[contig] = bin_id
        bin_ctg_dict[bin_id].append(contig)
    return ctg_bin_dict, bin_ctg_dict

def load_ctg_motifs_bk(profile_dir, MGE_dict):
    ctg_single_dict = {}
    single_ctg_dict = defaultdict(list)
    ctg_profile_dict = {}
    ctg_motif_dict = {}
    ## enumerage the files in the profile_dir
    for file in os.listdir(profile_dir):
        if file.endswith(".motifs.profile.csv"):
            ## extract the host name
            ctg_name = file.split(".")[0]
            if ctg_name in MGE_dict: ## skip other MGEs, only focus on chromosomal host
                continue
            host_motif = os.path.join(profile_dir,"../motifs", f"{ctg_name}.motifs.csv")
            if not os.path.exists(host_motif):
                print (f"{host_motif} does not exist.")
                continue
            motif_df = pd.read_csv(host_motif)
            ctg_motif_dict[ctg_name] = motif_df
            host_profile = os.path.join(profile_dir, file)
            ## check if the host profile and host_motif exists
            if not os.path.exists(host_profile):
                print (f"{host_profile} does not exist.")
                continue
            host_df = pd.read_csv(host_profile)
            ctg_single_dict[ctg_name] = ctg_name
            single_ctg_dict[ctg_name].append(ctg_name)
            ctg_profile_dict[ctg_name] = host_df
    return ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict


def process_file_with_args(args):
    file, profile_dir, MGE_dict = args
    if file.endswith(".motifs.profile.csv"):
        ctg_name = file.split(".")[0]
        if ctg_name in MGE_dict:  # Skip other MGEs, only focus on chromosomal host
            return None
        host_motif = os.path.join(profile_dir, "../motifs", f"{ctg_name}.motifs.csv")
        if not os.path.exists(host_motif):
            print(f"{host_motif} does not exist.")
            return None
        motif_df = pd.read_csv(host_motif)
        host_profile = os.path.join(profile_dir, file)
        if not os.path.exists(host_profile):
            print(f"{host_profile} does not exist.")
            return None
        host_df = pd.read_csv(host_profile)
        return {
            "ctg_name": ctg_name,
            "motif_df": motif_df,
            "host_df": host_df
        }
    return None

def load_ctg_motifs(profile_dir, MGE_dict, threads=4):
    # Use ProcessPoolExecutor to process files in parallel
    results = []
    with ProcessPoolExecutor(max_workers=threads) as executor:
        args = [(file, profile_dir, MGE_dict) for file in os.listdir(profile_dir)]
        for result in executor.map(process_file_with_args, args):
            if result is not None:
                results.append(result)

    # Combine results into dictionaries
    ctg_single_dict = {}
    single_ctg_dict = defaultdict(list)
    ctg_profile_dict = {}
    ctg_motif_dict = {}

    for result in results:
        ctg_name = result["ctg_name"]
        ctg_single_dict[ctg_name] = ctg_name
        single_ctg_dict[ctg_name].append(ctg_name)
        ctg_profile_dict[ctg_name] = result["host_df"]
        ctg_motif_dict[ctg_name] = result["motif_df"]

    return ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="get invasion score of MGE", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--work_dir", type=str, help="<str> work_dir", metavar="\b")
    # required.add_argument("--all_profiles", type=str, help="<str> separate by space", nargs="+", metavar="\b")
    optional.add_argument("--plasmid_file", type=str, help="<str> *_plasmid_summary.tsv by genomad.", metavar="\b")
    optional.add_argument("--plasmid", type=str, help="<str> plasmid contig name.", metavar="\b")
    optional.add_argument("--bin_file", type=str, help="<str> binning file: contig bin_id", metavar="\b")
    optional.add_argument("--min_frac", type=float, help="<str> output motif summary.", default=0.5, metavar="\b")
    ## add threads
    optional.add_argument("--threads", type=int, help="<int> number of threads.", default=1, metavar="\b")

    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())



    min_frac = args["min_frac"]
    # plasmid_contig = "E_coli_H10407_4"
    # work_dir = "/home/shuaiw/borg/bench/zymo6_NM3/"
    # genomad_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_genomad/merged2_summary/merged2_plasmid_summary.tsv"
    work_dir = args["work_dir"]
    profile_dir = os.path.join(work_dir, "profiles")
    host_dir = os.path.join(work_dir, "hosts")
    bin_file = args["bin_file"]
    os.makedirs(host_dir, exist_ok = True)


    if args["plasmid_file"]:
        batch_MGE_invade(args["plasmid_file"], profile_dir, host_dir, bin_file, min_frac, args["threads"])
    elif args["plasmid"]:
        plasmid_name = args["plasmid"]
        MGE_motif_num = count_MGE_with_motif(plasmid_name, profile_dir)
        print (f"Processing {plasmid_name} with {MGE_motif_num} original motifs.")
        ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict = load_ctg_motifs(profile_dir, {}, args["threads"])
        print (f"Loading is done, {len(single_ctg_dict)} contigs.", len(ctg_motif_dict))
        if bin_file is not None:
            if not os.path.exists(bin_file):
                print (f"{bin_file} does not exist.")
            ctg_bin_dict, bin_ctg_dict = load_bin(bin_file)
        else:
            ctg_bin_dict, bin_ctg_dict = ctg_single_dict, single_ctg_dict
        bin_df_dict, bin_motif_dict = merge_bin_motif(bin_ctg_dict, ctg_motif_dict, ctg_profile_dict)
        for_each_plasmid(bin_df_dict, bin_motif_dict, bin_ctg_dict, ctg_profile_dict, ctg_motif_dict, args["plasmid"], profile_dir, host_dir, min_frac, {})
    else:
        print ("Please provide either --plasmid_file or --plasmid.")

    ## extract the info saved in the host_dir
    summary_host(host_dir)



    """

    min_frac = 0.5
    # host_profile = "/home/shuaiw/borg/bench/zymo6_NM3/profiles/E_coli_H10407_1.motifs.profile.csv"
    host_profile = "/home/shuaiw/borg/all_test_ccs3/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_4094_L.motifs.profile.csv"
    host_motif = "/home/shuaiw/borg/all_test_ccs3/motifs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_4094_L.motifs.csv"
    plasmid_profile = "/home/shuaiw/borg/all_test_ccs3/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L.motifs.profile.csv"

    motif_data = extract_motif_data(host_profile, plasmid_profile, min_frac)
    print (motif_data)
    motif_data = filter_motifs(host_motif, motif_data)
    print (motif_data)

    result = invasion_score_from_counts(motif_data, min_frac)
    print(result)"
    """

# python cal_invasion_score.py --work_dir /home/shuaiw/borg/bench/zymo_new_ref --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_genomad/merged2_summary/merged2_plasmid_summary.tsv
# python cal_invasion_score.py  --work_dir /home/shuaiw/borg/all_test_ccs3 --plasmid_file /home/shuaiw/borg/contigs/borg_pure.txt
# python cal_invasion_score.py  --work_dir /home/shuaiw/borg/bench/all_ccs_1k --plasmid_file /home/shuaiw/borg/contigs/borg_pure.txt
# python cal_invasion_score.py  --work_dir /home/shuaiw/borg/bench/all_ccs_1k --plasmid_file /home/shuaiw/borg/contigs/genomad/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_summary/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_plasmid_summary.tsv
# python cal_invasion_score.py --work_dir /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30 --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list

