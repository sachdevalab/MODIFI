import profile
from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys, os, re
import numpy as np
from collections import defaultdict
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name




def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos, modified_loci):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1
    motif_sites = {}
    motif_loci_num = 0
    motif_modify_num = 0
    for_loci_num = 0
    rev_loci_num = 0
    for_modified_num = 0
    rev_modified_num = 0

    # motif_ipd_ratio = []
    record_modified_sites = {}

    for r, contig in REF.items():
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            tag = r + ":" + str(site+exact_pos) + "+"
            record_modified_sites[tag] = motif_new
            ##  also mark the nearby sites
            # for i in range(site-100, site+100):
            #     record_modified_sites[r + ":" + str(i) + "+"] = motif_new
            # print (tag)
            motif_loci_num += 1
            for_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                for_modified_num += 1
            # if tag in ipd_ratio_dict:
            #     motif_ipd_ratio.append(ipd_ratio_dict[tag][1])

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            record_modified_sites[tag] = motif_new
            # for i in range(site-100, site+100):
            #     record_modified_sites[r + ":" + str(i) + "+"] = motif_new
            motif_loci_num += 1
            rev_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                rev_modified_num += 1
            # if tag in ipd_ratio_dict:
            #     motif_ipd_ratio.append(ipd_ratio_dict[tag][1])
    # print ("for_loci_num", for_loci_num, for_modified_num, "forward modified ratio", for_modified_num/for_loci_num)
    # print ("rev_loci_num", rev_loci_num, rev_modified_num, "reverse modified ratio", rev_modified_num/rev_loci_num)

    ## plot the distribution of motif_ipd_ratio
    # print (len(motif_ipd_ratio))
    # plt.hist(motif_ipd_ratio, bins=100)
    # plt.xlabel("IPD ratio")
    # plt.ylabel("Frequency")
    # plt.title("IPD ratio distribution of motif " + motif_new)
    # plt.savefig("../tmp/" + motif_new + ".png")
    # plt.close()

    # print ("motif_loci_num", motif_loci_num)
    # print ("motif_modify_num", motif_modify_num)
    if motif_loci_num == 0:
        ratio = 0
    else:
        ratio = motif_modify_num/motif_loci_num
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

    # if ratio > 0.4 and proportion_all_modified > 0.1:
    #     print (motif_new, "modified_ratio", ratio, motif_modify_num, motif_loci_num, proportion_all_modified)

    return [for_loci_num, for_modified_num,for_ratio,\
            rev_loci_num, rev_modified_num, rev_ratio,\
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified], record_modified_sites

def get_modified_ratio(gff, score_cutoff = 30):
    ## read the gff file
    f = open(gff, "r")
    modified_loci = {}
    for line in f:
        if line[0] == "#":
            continue
        line = line.strip().split("\t")

        ref = line[0]
        pos = int(line[3]) 
        strand = line[6]
        score = int(line[5])
        if score <= score_cutoff:
            continue
        modified_loci[ref + ":" + str(pos) + strand] = score
    print ("no. of modified loci", len(modified_loci))
    return modified_loci

def read_ipd_ratio(ipd_ratio_file):
    ipd_ratio_dict = {}
    df = pd.read_csv(ipd_ratio_file, nrows = 1000)
    for index, row in df.iterrows():
        # print (row['refName'], row['tpl'], row['strand'], row['coverage'], row['ipd_ratio'])
        if row['strand'] == 1:
            strand_string = "-"
        else:
            strand_string = "+"
        tag = row['refName'] + ":" + str(row['tpl']+1) + strand_string
        ipd_ratio_dict[tag] = [row['ipd_ratio'], row['tMean']]
    # print (len(ipd_ratio_dict))
    ## print some of the ipd_ratio_dict
    for i, (k, v) in enumerate(ipd_ratio_dict.items()):
        # print (k, v)
        if i > 10:
            break
    return ipd_ratio_dict

def bioreactor():
    prefix_list = [["cow_bioreactor_1", "cow_bioreactor_1_636_C"], \
                   ["cow_bioreactor_2", "cow_bioreactor_2_1062_L"], \
                   ["cow_bioreactor_2", "cow_bioreactor_2_601_C"], \
                   ["cow_bioreactor_4", "cow_bioreactor_4_1750_C"],\
                    ["cow_bioreactor_5", "cow_bioreactor_5_1162_C"]]
    # prefix = "cow_bioreactor_4"
    # contig = "cow_bioreactor_4_1750_C"
    data = []
    motif_list = [["GATC", 2], ["ACNCAG", 5], ["GAAATC", 4], ["ACTNNNNNNRGTC", 1], ["GGCATC", 4]]
    for prefix, contig in prefix_list:
        my_ref = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/contigs/{contig}.fa"
        print (my_ref)
        gff = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/gffs/{contig}.gff"
        ipd_ratio_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/ipd_ratio/{contig}.ipd3.csv"
        REF = read_ref(my_ref)
        # print (REF)
        modified_loci = get_modified_ratio(gff)
        # motifs = pd.read_csv(all_motifs)
        # ipd_ratio_dict = read_ipd_ratio(ipd_ratio_file)
        
        for motif_new, exact_pos in motif_list:
            all_record = {}
            motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
            # print (motif_profile)
            data.append([contig,motif_new, motif_profile[-2]])
    df = pd.DataFrame(data, columns = ["contig", "motifString", "fraction"])
    ## plot heatmap with df
    df_pivot = df.pivot(index='contig', columns='motifString', values='fraction').fillna(0)
    import seaborn as sns
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(df_pivot, cmap="YlGnBu", annot=True, fmt=".2f", cbar_kws={'label': 'Fraction'})
    ax.set_title("Motif Fraction Heatmap")
    ax.set_xlabel("Motif String")
    ax.set_ylabel("Contig")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()
    plt.savefig("../../tmp/results2/motif_fraction_heatmap.pdf")

def find_species(closed_genome, GTDB_file, outdir, prefix):
    df = pd.read_csv(GTDB_file, sep="\t")
    df = df[df["closest_genome_reference"] == closed_genome]
    species_contigs = df["user_genome"].values
    contig_list = [[prefix, contig] for contig in species_contigs]
    #### collect all motifs of these contigs
    motif_set = set()
    for contig in species_contigs:
        motif_file = f"{outdir}/{prefix}_methylation3/motifs/{contig}.motifs.csv"
        ## check if file exists
        if not os.path.exists(motif_file):
            continue
        motif_df = pd.read_csv(motif_file)
        for index, row in motif_df.iterrows():
            motif = str(row["motifString"]) + "_" + str(row["centerPos"])
            motif_set.add(motif)
    return motif_set, contig_list

def profile_heatmap(prefix_list, motif_list, all_dir, tmp_res_file, cluster_obj, data_type="meta", score_cutoff = 30):
    data = []
    for prefix, contig in prefix_list:
        ctg_obj = My_contig(prefix, all_dir, contig, data_type)
        ## skip if files do not exist
        if not os.path.exists(ctg_obj.gff) or not os.path.exists(ctg_obj.ipd_ratio_file) or not os.path.exists(ctg_obj.ctg_ref):
            continue
        REF = read_ref(ctg_obj.ctg_ref)
        # print (REF)
        modified_loci = get_modified_ratio(ctg_obj.gff, score_cutoff)
        # motifs = pd.read_csv(all_motifs)
        # ipd_ratio_dict = read_ipd_ratio(ctg_obj.ipd_ratio_file)

        for motif_new, exact_pos in motif_list:
            motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
            # print (motif_profile)
            data.append([contig,motif_new + "_" + str(exact_pos), motif_profile[-2]])
    df = pd.DataFrame(data, columns = ["contig", "motifString", "fraction"])
    cluster_obj.get_profile(df, tmp_res_file)

    return cluster_obj

def clean_profile(prefix_list, motif_list, all_dir, data_type="meta", score_cutoff = 30):
    data = []
    for prefix, contig in prefix_list:
        ctg_obj = My_contig(prefix, all_dir, contig, data_type)
        if not os.path.exists(ctg_obj.gff) or not os.path.exists(ctg_obj.ipd_ratio_file) or not os.path.exists(ctg_obj.ctg_ref):
            continue
        REF = read_ref(ctg_obj.ctg_ref)
        modified_loci = get_modified_ratio(ctg_obj.gff, score_cutoff)


        for motif_new, exact_pos in motif_list:
            motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
            # print (motif_profile)
            data.append([contig,motif_new + "_" + str(exact_pos), motif_profile[-2], motif_profile[-4], motif_profile[-3]])
    df = pd.DataFrame(data, columns = ["contig", "motifString", "fraction", "motif_loci_num", "motif_modified_num"])
    return df
    
def given_species(all_dir, closed_genome, seq_dir):
    
    all_motif_set = set()
    all_contig_list = []
    plot_name = f"/home/shuaiw/borg/paper/motif_change/plot/motif_fraction_heatmap_{closed_genome}.pdf"
    close_genome_dir = os.path.join(seq_dir, closed_genome)
    if not os.path.exists(close_genome_dir):
        os.makedirs(close_genome_dir)
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        outdir = f"{all_dir}/{prefix}/"
        GTDB_file = f"{outdir}/GTDB/gtdbtk.all.summary.tsv"
        
        ## skip if GTDB_file does not exist
        if not os.path.exists(GTDB_file):
            continue
        motif_set, contig_list = find_species(closed_genome, GTDB_file, outdir, prefix)
        for prefix, contig in contig_list:
            ref = f"{outdir}/{prefix}_methylation3/contigs/{contig}.fa"
            # print (ref)
            ## copy ref to close_genome_dir
            if os.path.exists(ref):
                os.system(f"cp {ref} {close_genome_dir}/")
        all_motif_set.update(motif_set)
        all_contig_list.extend(contig_list)
        # if len(all_contig_list) > 5:
        #     break

    motif_list = []
    for motif in all_motif_set:
        motif_name, motif_pos = motif.split("_")
        motif_list.append([motif_name, int(motif_pos)])
    print (motif_list)
    print (all_contig_list)
    if len(all_contig_list) > 1:
        profile_heatmap(all_contig_list, motif_list, all_dir, plot_name, closed_genome)

def collect_species(all_dir):
    genomes_list = []
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        outdir = f"{all_dir}/{prefix}/"
        GTDB_file = f"{outdir}/GTDB/gtdbtk.all.summary.tsv"
        if not os.path.exists(GTDB_file):
            continue
        df = pd.read_csv(GTDB_file, sep="\t")
        for index, row in df.iterrows():
            # if re.search("s__", row["classification"]):
            #     continue
            if len(row["classification"].split(";")) < 7:
                continue
            ## skip if closest_genome_reference is nan
            if pd.isna(row["closest_genome_reference"]):
                continue
            if row["closest_genome_reference"] in genomes_list:
                continue
            genomes_list.append(row["closest_genome_reference"])
    return genomes_list

def by_GTDB():
    closed_genome =  "GCF_900104055.1" #"GCF_000735435.1"
    all_dir = "/home/shuaiw/borg/paper/run2/"
    seq_dir = "/home/shuaiw/borg/paper/motif_change/seq/"
    genomes_list = collect_species(all_dir)
    print (genomes_list)
    print (len(genomes_list))

    i = 0
    for closed_genome in genomes_list[111:]:

        print (f"Processing {closed_genome}...", i)
        given_species(all_dir, closed_genome, seq_dir)
        i += 1

def read_drep_cluster(drep_clu_file, depth_dict):
    
    drep_clu_dict = defaultdict(list)
    df = pd.read_csv(drep_clu_file)
    for index, row in df.iterrows():
        # drep_clu_dict[row['genome'][:-3]] = row['secondary_cluster']
        contig = row['genome'][:-3]
        # if contig not in depth_dict:
        #     continue
        drep_clu_dict[row['secondary_cluster']].append(row['genome'][:-3])
    return drep_clu_dict

def find_species_drep(contig, all_dir, prefix, data_type="meta", min_frac=0.3, min_sites=100):
    species_contigs = [contig]
    contig_list = [[prefix, contig] for contig in species_contigs]
    #### collect all motifs of these contigs
    motif_set = set()
    for contig in species_contigs:
        ctg_obj = My_contig(prefix, all_dir, contig, data_type)
        print (ctg_obj.motif_file)
        motif_df = ctg_obj.read_motif(min_frac, min_sites)
        if motif_df is None:
            continue
        for index, row in motif_df.iterrows():
            motif = str(row["motifString"]) + "_" + str(row["centerPos"])
            motif_set.add(motif)
    return motif_set, contig_list

def given_species_drep(all_dir, members, seq_dir, cluster, fig_dir, 
                       tmp_res, data_type="meta", min_frac=0.3, min_sites=100, score_cutoff = 30):
    
    all_motif_set = set()
    all_contig_list = []
    plot_name = f"{fig_dir}/{cluster}.pdf"
    tmp_res_file = f"{tmp_res}/{cluster}.csv"
    cluster_obj = My_cluster(cluster, members)
    if isinstance(members[0], str):
        new_members = [(contig,"_".join(contig.split("_")[:-2])) for contig in members]
        # print ("&&&", new_members)
        members = new_members
    print (members)
    for contig, prefix in members:
        # print (contig)
        # prefix = "_".join(contig.split("_")[:-2])
        ctg_obj = My_contig(prefix, all_dir, contig, data_type)
        ref = ctg_obj.ctg_ref

        close_genome_dir = os.path.join(seq_dir, cluster)
        if not os.path.exists(close_genome_dir):
            os.makedirs(close_genome_dir)

        if os.path.exists(ref):
            os.system(f"cp {ref} {close_genome_dir}/")

        motif_set, contig_list = find_species_drep(contig, all_dir, prefix, data_type, min_frac, min_sites)

        all_motif_set.update(motif_set)
        all_contig_list.extend(contig_list)

    motif_list = []
    for motif in all_motif_set:
        motif_name, motif_pos = motif.split("_")
        motif_list.append([motif_name, int(motif_pos)])
    if cluster == "161_4":
        motif_list = [['CRTANNNNNNRTG', 4], ['ATGCAT', 5], ['CAANNNNNNNTAYG', 3], \
            ['TCAYNNNNNNTTG', 3], ['CAANNNNNNRTGA', 3], \
                ['CRTANNNNNNNTTG', 4], ['CAYNNNNNNTAYG', 2]]
        all_contig_list += [["infant_14", "infant_14_230_C"], ["infant_14", "infant_14_88_C"],\
                            ["infant_2", "infant_2_267_C"], ["infant_2", "infant_2_60_C"],\
                            ["infant_3", "infant_3_289_C"], ["infant_3", "infant_3_287_C"],\
                            ["infant_4", "infant_4_271_C"], ["infant_4", "infant_4_61_C"]]
    print ("motif_list", motif_list)
    # print ("all_contig_list", all_contig_list)
    # if len(all_contig_list) > 1:
    cluster_obj = profile_heatmap(all_contig_list, motif_list, all_dir, tmp_res_file, cluster_obj,data_type, score_cutoff)
    
    return cluster_obj

def depth_filter(all_dir, min_depth = 10):
    depth_dict = {}
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        depth_file = f"{all_dir}/{prefix}/{prefix}_methylation3/mean_depth.csv"
        if not os.path.exists(depth_file):
            continue
        depth_df = pd.read_csv(depth_file)
        
        for index, row in depth_df.iterrows():
            if row['depth'] < min_depth:
                continue
            depth_dict[row['contig']] = row['depth']
    return depth_dict

def plot_variation_fraction(variation_data, paper_fig_dir):
    ## count the proportion of clusters with motif variation
    ## x-axis shows cluster member cutoff, for 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    cutoff_dict = defaultdict(lambda: [0,0]) # key: cutoff, value: [total clusters, clusters with variation]
    for cluster, member_num, motif_variation_flag in variation_data:
        for cutoff in range(2, 11):
            if member_num >= cutoff:
                cutoff_dict[cutoff][0] += 1
                if motif_variation_flag == "variation":
                    cutoff_dict[cutoff][1] += 1
    ## plot the proportion barplot
    data = []
    for cutoff in range(2, 11):
        total_clusters, variation_clusters = cutoff_dict[cutoff]
        if total_clusters == 0:
            proportion = 0
        else:
            proportion = variation_clusters / total_clusters
        data.append([cutoff, proportion, variation_clusters, total_clusters - variation_clusters])
    
    ## Save data to CSV
    df = pd.DataFrame(data, columns=['cutoff', 'proportion', 'variation_clusters', 'no_variation_clusters'])
    df.to_csv(f"{paper_fig_dir}/motif_variation_data.csv", index=False)
    
    plt.figure(figsize=(4,4))
    sns.barplot(x=df['cutoff'], y=df['proportion'], color="skyblue")
    plt.xlabel("Cluster member cutoff")
    plt.ylabel("Proportion of clusters with motif variation")
    # plt.ylim(0,1)
    # plt.title("Proportion of clusters with motif variation vs. cluster member cutoff")
    plt.savefig(f"{paper_fig_dir}/motif_variation_proportion.pdf")
    ## also plot the stacked barplot of all clusters
    plt.figure(figsize=(4,4))
    all_clusters_array = np.array([[row[2], row[3]] for row in data])
    bottom = np.zeros(len(all_clusters_array))
    labels = ["Variation", "No Variation"]
    colors = ["skyblue", "lightgray"]
    for i in range(all_clusters_array.shape[1]):
        plt.bar(df['cutoff'], all_clusters_array[:, i], bottom=bottom, color=colors[i], label=labels[i])
        bottom += all_clusters_array[:, i]
    plt.xlabel("Cluster member cutoff")
    plt.ylabel("Number of clusters")
    plt.legend()
    plt.savefig(f"{paper_fig_dir}/motif_variation_stacked_bar.pdf")

def plot_similarity_dist(similarity_data, paper_fig_dir):
    ## similarity_data: list of [cosine_similarity, jaccard_similarity]
    cosine_similarities = [x[0] for x in similarity_data]
    jaccard_similarities = [x[1] for x in similarity_data]
    plt.figure(figsize=(8,4))
    plt.subplot(1,2,1)
    sns.histplot(cosine_similarities, bins=30, kde=True, color="skyblue")
    plt.xlabel("Cosine Similarity")
    plt.ylabel("Frequency")
    plt.title("Distribution of Cosine Similarity between pairs of contigs")
    plt.subplot(1,2,2)
    sns.histplot(jaccard_similarities, bins=30, kde=True, color="lightgreen")
    plt.xlabel("Jaccard Similarity")
    plt.ylabel("Frequency")
    plt.title("Distribution of Jaccard Similarity between pairs of contigs")
    plt.tight_layout()
    plt.savefig(f"{paper_fig_dir}/similarity_distributions.pdf")

def plot_clade_size(clade_data, paper_fig_dir):
    clade_df = pd.DataFrame(clade_data, columns=["clade_num", "cluster", "member_num", "phylum", "species"])
    ## plot distribution of clade sizes
    plt.figure(figsize=(4,4))
    sns.histplot(clade_df['clade_num'], bins=30, kde=True, color="salmon")
    plt.xlabel("Number of Motif Clades")
    plt.ylabel("Frequency")
    plt.title("Distribution of Motif Clade Sizes across Clusters")
    plt.savefig(f"{paper_fig_dir}/clade_size_distribution.pdf")
    ## select the top 20 species with highest clade_num, plot barplot, show species name, 
    # cluster name,  and member_num. rotate whole plot for better visualization
    top20_df = clade_df.sort_values(by="clade_num", ascending=False).head(10)
    
    # Create labels that include cluster name and member count
    top20_df['combined_label'] = top20_df.apply(lambda row: f"{row['species']}\n({row['cluster']}, n={row['member_num']})", axis=1)

    plt.figure(figsize=(8,4))
    ax = sns.barplot(x="clade_num", y="combined_label", data=top20_df, palette="viridis")
    
    # Adjust layout to accommodate longer labels
    plt.subplots_adjust(left=0.4)  # Make more space for y-axis labels
    
    # Set font size for better readability
    plt.tick_params(axis='y', labelsize=10)
    plt.tick_params(axis='x', labelsize=10)
    
    plt.xlabel("Number of Motif Clades", fontsize=12)
    plt.ylabel("Species (Cluster, Members)", fontsize=12)
    plt.title("Top 10 Species with Highest Motif Clade Numbers", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{paper_fig_dir}/top20_species_clade_numbers.pdf", bbox_inches='tight')

def collect_represent_ctgs(dereplicated_genomes_dir):
    represent_ctg_set = set()
    for file in os.listdir(dereplicated_genomes_dir):
        if file.endswith(".fa"):
            contig = file[:-3]
            represent_ctg_set.add(contig)
    print (len(represent_ctg_set), "representative contigs collected")
    return represent_ctg_set

def select_represent(represent_ctg_set, members):
    for contig in members:
        if contig in represent_ctg_set:
            return contig
    return None

def main_meta():
    all_dir = "/home/shuaiw/borg/paper/run2/"
    ctg_taxa_dict = get_ctg_taxa(all_dir)
    # print (ctg_taxa_dict)
    # ANI = 95
    # for ANI in [95, 97, 98, 99]:
    for ANI in [99]:
        drep_clu_file = f"/home/shuaiw/borg/paper/specificity/dRep_{ANI}_out/data_tables/Cdb.csv"
        dereplicated_genomes_dir = f"/home/shuaiw/borg/paper/specificity/dRep_{ANI}_out/dereplicated_genomes/"
        seq_dir = "/home/shuaiw/borg/paper/motif_change/seq_drep/"
        fig_dir = f"/home/shuaiw/borg/paper/motif_change/plot_drep2_{ANI}/"
        tmp_res = f"/home/shuaiw/borg/paper/motif_change/result_drep2_{ANI}/"
        paper_fig_dir = f"../../tmp/figures/strain_diff/drep_{ANI}/"
        ## create fig_dir and tmp_res if not exist
        os.makedirs(fig_dir, exist_ok=True)
        os.makedirs(tmp_res, exist_ok=True)
        os.makedirs(paper_fig_dir, exist_ok=True)


        # depth_dict = depth_filter(all_dir, min_depth = 10)
        drep_clu_dict = read_drep_cluster(drep_clu_file, {})
        represent_ctg_set = collect_represent_ctgs(dereplicated_genomes_dir)
        # drep_clu_dict = read_drep_cluster(drep_clu_file, depth_dict)
        
        print (len(drep_clu_dict), "drep clusters")
        ## count how many clusters have more than 10 members
        count = 0
        cutoff = 1
        for cluster, members in drep_clu_dict.items():
            if len(members) > cutoff:
                count += 1
        print (count, f"clusters have more than {cutoff} members")

        variation_data = []
        similarity_data = []
        clade_data = []
        for cluster, members in drep_clu_dict.items():
            if cluster != "161_4":
                continue
            if len(members) > cutoff:
                print ("cluster", cluster, len(members), len(variation_data))

                cluster_obj = given_species_drep(all_dir, members, seq_dir, cluster,
                                                fig_dir, tmp_res, min_frac=0.5, min_sites=500)
                
                # cluster_obj = My_cluster(cluster, members) 
                # cluster_obj.load_df(tmp_res)
                # cluster_obj.manual_filter_motifs()

                if len(cluster_obj.profile_df) < 2:
                    print ("skip cluster with less than 2 contigs with motif profiles")
                    continue

                represent_ctg = select_represent(represent_ctg_set, members)
                represent_ctg_lineage = ctg_taxa_dict[represent_ctg] if represent_ctg in ctg_taxa_dict else "NA"
                cluster_phylum = classify_taxa(represent_ctg_lineage, "phylum")
                cluster_species = classify_taxa(represent_ctg_lineage, "species")
                cluster_name = get_detail_taxa_name(represent_ctg_lineage)

                motif_variation_flag = cluster_obj.check_diff_motifs()

                
                similarity_data_cluster, motif_clade_num = cluster_obj.pairwise_compare(bin_freq=0.6)
                variation_data.append([cluster, len(members), motif_variation_flag])
                similarity_data += similarity_data_cluster
                clade_data += [[motif_clade_num, cluster, len(members), cluster_phylum, cluster_name]]
                print ("###############################", represent_ctg_lineage)
                if motif_variation_flag == "variation" or motif_clade_num > 1:
                    plot_name = f"{fig_dir}/{cluster}.pdf"
                    cluster_obj.plot_profile(cluster, plot_name, cluster_name)
                # if len(variation_data) > 10:
                #     break
        # plot_variation_fraction(variation_data, paper_fig_dir)
        # plot_similarity_dist(similarity_data, paper_fig_dir)
        # plot_clade_size(clade_data, paper_fig_dir)
        # print ("all done")

def main_isolation():
    all_dir = "/home/shuaiw/borg/paper/isolation/batch2_results/"
    ctg_taxa_dict = get_ctg_taxa(all_dir, 'isolation')
    # print (ctg_taxa_dict)
    # ANI = 95
    # for ANI in [95, 97, 98, 99]:
    for ANI in [99]:
        drep_clu_file = f"/home/shuaiw/borg/paper/specificity/iso_{ANI}_out/data_tables/Cdb.csv"
        dereplicated_genomes_dir = f"/home/shuaiw/borg/paper/specificity/iso_{ANI}_out/dereplicated_genomes/"
        seq_dir = "/home/shuaiw/borg/paper/motif_change/iso_seq_drep/"
        edit_dir = "/home/shuaiw/borg/paper/motif_change/iso_edit_dist/"
        fig_dir = f"/home/shuaiw/borg/paper/motif_change/iso_plot_drep2_{ANI}/"
        tmp_res = f"/home/shuaiw/borg/paper/motif_change/iso_result_drep2_{ANI}/"
        paper_fig_dir = f"../../tmp/figures/strain_diff/iso_drep_{ANI}/"
        ## create fig_dir and tmp_res if not exist
        os.makedirs(fig_dir, exist_ok=True)
        os.makedirs(tmp_res, exist_ok=True)
        os.makedirs(paper_fig_dir, exist_ok=True)

        """
        # depth_dict = depth_filter(all_dir, min_depth = 10)
        drep_clu_dict = read_drep_cluster(drep_clu_file, {})
        represent_ctg_set = collect_represent_ctgs(dereplicated_genomes_dir)
        # drep_clu_dict = read_drep_cluster(drep_clu_file, depth_dict)
        
        print (len(drep_clu_dict), "drep clusters")
        ## count how many clusters have more than 10 members
        count = 0
        cutoff = 1
        for cluster, members in drep_clu_dict.items():
            if len(members) > cutoff:
                count += 1
        print (count, f"clusters have more than {cutoff} members")

        variation_data = []
        similarity_data = []
        clade_data = []
        cluster_id = 0
        dnadiff_list_all = []
        dna_motif_corr = []
        for cluster, members in drep_clu_dict.items():
            if len(members) > cutoff:
                # print ("cluster", cluster, len(members), len(variation_data))

                # cluster_obj = given_species_drep(all_dir, members, seq_dir, cluster,
                #                                 fig_dir, tmp_res, "isolation", min_frac=0.5, min_sites=500)
                
                cluster_obj = My_cluster(cluster, members) 
                dnadiff_list, dnadiff_mat = cluster_obj.pairwise_edit_distance(seq_dir, edit_dir)
                print ("#####cluster_id", cluster, cluster_id, len(drep_clu_dict))
                cluster_id += 1
                # break
                dnadiff_list_all += dnadiff_list
                cluster_obj.load_df(tmp_res)
                print (">>>>", cluster_obj.profile_df)
                cluster_obj.manual_filter_motifs()

                if len(cluster_obj.profile_df) < 2:
                    print ("skip cluster with less than 2 contigs with motif profiles")
                    continue

                ### get taxa of cluster
                represent_ctg = select_represent(represent_ctg_set, members)
                represent_ctg_lineage = ctg_taxa_dict[represent_ctg] if represent_ctg in ctg_taxa_dict else "NA"
                cluster_phylum = classify_taxa(represent_ctg_lineage, "phylum")
                cluster_species = classify_taxa(represent_ctg_lineage, "species")

                motif_variation_flag = cluster_obj.check_diff_motifs()
                similarity_data_cluster, motif_clade_num, jaccard_matrix = cluster_obj.pairwise_compare(bin_freq=0.6)
                dna_motif_corr = update_corr(dnadiff_mat, jaccard_matrix, dna_motif_corr, cluster_obj)
                # variation_data.append([cluster, len(members), motif_variation_flag])
                # similarity_data += similarity_data_cluster
                # clade_data += [[motif_clade_num, cluster, len(members), cluster_phylum, cluster_species]]

                # if motif_variation_flag == "variation" or motif_clade_num > 1:
                #     plot_name = f"{fig_dir}/{cluster}.pdf"
                #     cluster_obj.plot_profile(cluster, plot_name, cluster_species)
                # if len(variation_data) > 10:
                #     break
        # plot_variation_fraction(variation_data, paper_fig_dir)
        # plot_similarity_dist(similarity_data, paper_fig_dir)
        # plot_clade_size(clade_data, paper_fig_dir)
        ## dna_motif_corr to df 
        dna_motif_corr_df = pd.DataFrame(dna_motif_corr, columns=["cluster", "contig1", "contig2", "dnadiff", "jaccard"])
        dna_motif_corr_df.to_csv(os.path.join(paper_fig_dir, "dna_motif_corr.csv"), index=False)
        """
        dna_motif_corr_df = pd.read_csv(os.path.join(paper_fig_dir, "dna_motif_corr.csv"))
        plot_dna_motif_corr(dna_motif_corr_df, paper_fig_dir)
        
        # print ("all done")
    # batch_dnadiff(dnadiff_list_all)
    
# 

def update_corr(dnadiff_mat, jaccard_matrix, dna_motif_corr, cluster_obj):
    for i in range(len(cluster_obj.members)):
        for j in range(i+1, len(cluster_obj.members)):
            dna_motif_corr.append((cluster_obj.cluster, cluster_obj.members[i], cluster_obj.members[j], dnadiff_mat[i][j], jaccard_matrix[i][j]))
            if dnadiff_mat[i][j] == 0 and jaccard_matrix[i][j] < 1:
                print ("!!!>>>>>", cluster_obj.cluster, cluster_obj.members[i], cluster_obj.members[j],\
                        dnadiff_mat[i][j], jaccard_matrix[i][j])
    return dna_motif_corr

def plot_dna_motif_corr(dna_motif_corr_df, paper_fig_dir):
    
    output_file = os.path.join(paper_fig_dir, "dna_motif_corr.pdf")
    if len(dna_motif_corr_df) == 0:
        print(f"[⚠️] No data to plot for {output_file}")
        return
    dnadiff_values = dna_motif_corr_df["dnadiff"]
    jaccard_values = dna_motif_corr_df["jaccard"]
    plt.figure(figsize=(6, 6))
    plt.scatter(dnadiff_values, jaccard_values, alpha=0.5)
    plt.xlabel("DNA Edit Distance")
    plt.ylabel("Motif Jaccard Index")
    plt.title("DNA vs Motif Correlation")
    plt.savefig(output_file)
    plt.close()

def batch_dnadiff(dnadiff_list_all):
    ## construct a batch script to run all dnadiff commands, use 10 scripts
    n_scripts = 30
    h = open("batch_all.sh", "w")
    scripts = [[] for _ in range(n_scripts)]
    for i, cmd in enumerate(dnadiff_list_all):
        scripts[i % n_scripts].append(cmd)
    for i, script_cmds in enumerate(scripts):
        script_file = f"scripts/dnadiff_batch_{i}.sh"
        with open(script_file, "w") as f:
            f.write("#!/bin/bash\n")
            for cmd in script_cmds:
                f.write(cmd + "\n")
        # os.system(f"chmod +x {script_file}")
        print (f"nohup bash scripts/dnadiff_batch_{i}.sh &", file = h)
    h.close()

def main_asthma():
    all_dir = "/home/shuaiw/borg/paper/run2/"
    seq_dir = "/home/shuaiw/borg/paper/asthma/B_longum"
    B_longum_list = "/home/shuaiw/borg/paper/asthma/test/B_longum.list"
    members = []
    for line in open(B_longum_list):
        contig = line.strip().split("\t")[0]
        members.append(contig)

    # members = ["asthma_11_4_C"]
    cluster = "profile"
    plot_name = os.path.join(seq_dir, f"borg_motif_profile.pdf")
    cluster_species = "B_longum"

    cluster_obj = given_species_drep(all_dir, members, seq_dir, cluster,
                                    seq_dir, seq_dir, min_frac=0.2, 
                                    min_sites=10, score_cutoff = 20)
    cluster_obj.plot_profile(cluster, plot_name, cluster_species)

if __name__ == "__main__":
    # main_meta()
    main_isolation()
    # main_asthma()