import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity
import re
import sys
from Bio.Seq import Seq
from scipy.stats import ttest_ind

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import My_gene, get_unique_motifs, get_detail_taxa_name, My_sample, get_ctg_taxa, classify_taxa, My_contig, Isolation_sample

PHYLUM_COLORS = {
    "Pseudomonadota":    "#d8b365",  # tan / brown
    "Bacillota":         "#f46d43",  # orange
    "Bacillota_A":       "#8da0cb",  # blue-lavender
    "Desulfobacterota":       "#fc8d62",  # light orange / salmon (same family as Bacillota)
    "Actinomycetota":    "#66c2a5",  # teal-green
    "Bacteroidota":      "#e78ac3",  # pink
    "Campylobacterota":  "#a6d854",  # light green
    "Acidobacteriota":   "#1b9e77",  # dark green
    "Verrucomicrobiota": "#7570b3",  # muted purple
    "Chloroflexota":   "#e7298a",  # magenta
    # Archaea phyla
    "Thermoproteota":    "#b15928",  # brown/rust
    "Thermoplasmatota":  "#fdbf6f",  # light orange/yellow
    "Halobacteriota":    "#ff7f00",  # bright orange
    "Aenigmatarchaeota": "#cab2d6",  # light purple
    "Nanoarchaeota":     "#6a3d9a",  # deep purple
    "Hadarchaeota":      "#fb9a99",  # light red/pink
    "Asgardarchaeota":   "#e31a1c",  # red
    "Micrarchaeota":     "#a6cee3",  # light blue

    "Others":            "#e6e6e6"   # very light neutral (only gray)
}



BAC_PHYLUM_COLORS = {
    "Pseudomonadota":    "#d8b365",  # tan / brown
    "Bacillota":         "#f46d43",  # orange
    "Bacillota_A":       "#8da0cb",  # blue-lavender
    "Desulfobacterota":       "#fc8d62",  # light orange / salmon (same family as Bacillota)
    "Actinomycetota":    "#66c2a5",  # teal-green
    "Bacteroidota":      "#e78ac3",  # pink
    "Campylobacterota":  "#a6d854",  # light green
    "Acidobacteriota":   "#1b9e77",  # dark green
    "Verrucomicrobiota": "#7570b3",  # muted purple
    "Chloroflexota":   "#e7298a",  # magenta
    "Bacillota_C":    "#b15928",  # brown/rust
    "Nitrospirota":  "#fdbf6f",  # light orange/yellow
    "Bacillota_I":    "#ff7f00",  # bright orange
    "Desulfobacterota_B": "#cab2d6",  # light purple
    "Patescibacteria":     "#6a3d9a",  # deep purple
    "Gemmatimonadota":      "#fb9a99",  # light red/pink
    "Others":            "#e6e6e6"   # very light neutral (only gray)
}


ARC_PHYLUM_COLORS = {
    # Archaea phyla
    "Thermoproteota":    "#b15928",  # brown/rust
    "Thermoplasmatota":  "#fdbf6f",  # light orange/yellow
    "Halobacteriota":    "#ff7f00",  # bright orange
    "Aenigmatarchaeota": "#cab2d6",  # light purple
    "Nanoarchaeota":     "#6a3d9a",  # deep purple
    "Hadarchaeota":      "#fb9a99",  # light red/pink
    "Asgardarchaeota":   "#e31a1c",  # red
    "Micrarchaeota":     "#a6cee3",  # light blue
    "Nanobdellota":    "#66c2a5",  # teal-green
    "Methanobacteriota":      "#e78ac3",  # pink
    "B1Sed10-29":  "#a6d854",  # light green
    "Others":            "#e6e6e6"   # very light neutral (only gray)

}


header = """TREE_COLORS
#use this template to define branch colors and styles, colored ranges and label colors/font styles/backgrounds
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
SEPARATOR SPACE
#SEPARATOR COMMA

#First 3 fields define the node, type and color
#Possible types are:
#'range': defines a colored range (colored background for labels/clade)
#'clade': defines color/style for all branches in a clade
#'branch': defines color/style for a single branch
#'label': defines font color/style for the leaf label
#'label_background': defines the leaf label background color

#The following additional fields are required:
#for 'range', field 4 defines the colored range label (used in the legend)

#The following additional fields are optional:
#for 'label', field 4 defines the font style ('normal',''bold', 'italic' or 'bold-italic') and field 5 defines the numeric scale factor for the font size (eg. with value 2, font size for that label will be 2x the standard size)
#for 'clade' and 'branch', field 4 defines the branch style ('normal' or 'dashed') and field 5 defines the branch width scale factor (eg. with value 0.5, branch width for that clade will be 0.5 the standard width)

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR
"""

Label_header = """LABELS
#use this template to change the leaf labels, or define/change the internal node names
#additionally, you can specify a custom class for internal nodes, which can be used to automatically collapse them
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#NODE_ID,LABEL,CLASS

#Examples

#note that the class field is optional

#define a name and class for an internal node. Class 'kingdom' will be available when using the #automatic clade collapsing function
#9031|9606,Metazoa,kingdom

#change the label for leaf node 9606\n"""

gradient_header="""
DATASET_GRADIENT
#In gradient datasets, each ID is associated to a single numeric value which is converted to a colored box based on the gradient defined.

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throught this file.
#SEPARATOR TAB
SEPARATOR SPACE
#SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL label 1

#dataset color (can be changed later)
COLOR #ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#


#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark


#LEGEND_TITLE Dataset legend
#LEGEND_SCALE 1
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_HORIZONTAL 0
#LEGEND_SHAPES 1 2 3
#LEGEND_COLORS #ff0000 #00ff00 #0000ff
#LEGEND_LABELS value1 value2 value3
#LEGEND_SHAPE_SCALES 1 1 0.5

#width of the gradient strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the gradient strip
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#automatically create and display a legend based on the color gradient defined below
#AUTO_LEGEND 1

#define the gradient colors. Values in the dataset will be mapped onto the corresponding color gradient.
#COLOR_MIN #ff0000
#COLOR_MAX #0000ff

#you can specify a gradient with three colors (e.g red to yellow to green) by setting 'USE_MID_COLOR' to 1, and specifying the midpoint color
#USE_MID_COLOR 1
#COLOR_MID #ffff00

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 1

#display or hide the dataset label above the gradient strip
#SHOW_LABELS 1

#text label size factor
#SIZE_FACTOR 1

#text label rotation
#LABEL_ROTATION 0

#text label shift in pixels (positive or negative)
#LABEL_SHIFT 0

#align the dataset label to the tree circle; only applies in circular display mode
#LABEL_ALIGN_TO_TREE,0

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#ID2 value2
#9606 10000
#LEAF1|LEAF2 11000"""

strip_heaer="""DATASET_COLORSTRIP
#In colored strip datasets, each ID is associated to a color box/strip and can have an optional label. Color can be specified in hexadecimal, RGB or RGBA notation. When using RGB or RGBA notation, you cannot use COMMA as the dataset separator

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
#SEPARATOR COMMA
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL label1

#dataset color (can be changed later)
COLOR #ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#If COLOR_BRANCHES is set to 1, branches of the tree will be colored according to the colors of the strips above the leaves.
#When all children of a node have the same color, it will be colored the same, ie. the color will propagate inwards towards the root.
COLOR_BRANCHES 0


#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

#LEGEND_TITLE Dataset_legend
#LEGEND_SCALE 1
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_HORIZONTAL 0
#LEGEND_SHAPES 1 1 2 2
#LEGEND_COLORS #ff0000 #00ff00 rgba(0,255,0,0.5) #0000ff
#LEGEND_LABELS value1 value2 value3 value4
#LEGEND_SHAPE_SCALES 1 1 0.5 1

#width of the colored strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the color strip 
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#if set to 1, border will be drawn completely around each colored strip box
#COMPLETE_BORDER 0

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 0


#display or hide the individual label inside each colored strip (when defined in the data below)
#SHOW_STRIP_LABELS 1

#position of the strip label within the box; 'top', 'center' or 'bottom'
#STRIP_LABEL_POSITION center

#strip label size factor (relative to the tree leaf labels)
#STRIP_LABEL_SIZE_FACTOR 1


#rotation of the strip labels; used only in rectangular tree display mode
#STRIP_LABEL_ROTATION 0

#strip label shift in pixels (positive or negative)
#STRIP_LABEL_SHIFT 0

#STRIP_LABEL_COLOR #000000

#draw a black outline around the text (width in pixels)
#STRIP_LABEL_OUTLINE 0.5

#calculate the label color automatically (black or white), based on the darkness of the color strip
#STRIP_LABEL_AUTO_COLOR 0

#display or hide the dataset label above the colored strip
#SHOW_LABELS 1

#dataset label size factor
#SIZE_FACTOR 1

#dataset label rotation
#LABEL_ROTATION 0

#dataset label shift in pixels (positive or negative)
#LABEL_SHIFT 0

#align the dataset label to the tree circle; only applies in circular display mode
#LABEL_ALIGN_TO_TREE,0

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages

#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA

#Examples:
#assign a red colored strip to leaf 9606, with label 'Human'
#9606 #ff0000 Human

#assign a green, semi-transparent (alpha 0.5) strip to an internal node, without any label. If 'Show internal values' is set to 'No', this will only be displayed if the node is collapsed. 
#9606|5664 rgba(0,255,0,0.5)

#soil_1_886_C #377eb8 Soil
"""



def get_best_ctg(depth_file, fai, min_len = 1000000):
    """
    Get the best contig based on length from a fasta file.
    """
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    length_dict = {}
    for index, row in depth_df.iterrows():

        if row['depth'] >= 10:
            good_depth[row['contig']] = row['depth']
    best_ctgs = {}
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            length_dict[ctg] = length
            if ctg[-1] == "C" and length >= min_len and ctg in good_depth:
                best_ctgs[ctg] = length
    # print (f"Total {len(best_ctgs)} contigs with length >= {min_len} found.")
    return good_depth, length_dict

def get_drep_motifs():
    drep_motifs = "/home/shuaiw/borg/bench/soil/run1/test.motifs_drep.csv"
    drep_motif_df = pd.read_csv(drep_motifs)
    drep_motifs_dict = {}
    for index, row in drep_motif_df.iterrows():
        motif = row['motifString']
        if motif not in drep_motifs_dict:
            drep_motifs_dict[motif] = 1
    return drep_motifs_dict

def plot_strain_heatmap(profile_file, min_frac=0.4):
    profiles = pd.read_csv(profile_file, index_col=0)
    ## only keep the columns that are in the best contigs
    best_ctgs = get_best_ctg()
    multiple_strain_ctg, ctg_phylum = read_drep_cluster()
    ## retain only the contigs that are in multiple strain clusters
    best_ctgs = [ctg for ctg in best_ctgs if ctg in multiple_strain_ctg]
    print (f"Total {len(best_ctgs)} contigs")
    
    profiles = profiles.loc[:, profiles.columns.isin(best_ctgs)]
    print ("original shape", profiles.shape)                                    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    print ("filtered shape 1", profiles.shape)
    print (profiles)
    drep_motifs_dict = get_drep_motifs()


    ## filter profiles to keep only the rows with motifString  that are in drep_motifs_dict
    profiles = profiles[profiles.index.isin(drep_motifs_dict.keys())]
    print ("filtered shape 2", profiles.shape)


    df = profiles.T
    # ctg_phylum = get_phylum(best_ctgs)

    phylum_colors = sns.color_palette("husl", len(set(ctg_phylum.values())))
    phylum_color_dict = {phylum: color for phylum, color in zip(set(ctg_phylum.values()), phylum_colors)}
    phylum_colors_list = [phylum_color_dict[ctg_phylum[ctg]] for ctg in df.index]

    # Use clustermap instead of heatmap
    g = sns.clustermap(
        df,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(30, 15),
        row_colors=phylum_colors_list,
        xticklabels=True,
        yticklabels=True
    )

    # Add legend for phylum colors
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=phylum_color_dict[phylum], label=phylum) for phylum in phylum_color_dict]
    plt.legend(handles=handles, title="Strain", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


    plt.savefig("../../tmp/results/heatmap_strain_cluster.png", dpi=300, bbox_inches='tight')
    plt.clf()

def plot_heatmap(profile_file, min_frac=0.4):
    profiles = pd.read_csv(profile_file, index_col=0)
    ## only keep the columns that are in the best contigs
    best_ctgs = get_best_ctg()
    
    profiles = profiles.loc[:, profiles.columns.isin(best_ctgs)]
    print ("original shape", profiles.shape)                                    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    print ("filtered shape 1", profiles.shape)
    print (profiles)
    drep_motifs_dict = get_drep_motifs()


    ## filter profiles to keep only the rows with motifString  that are in drep_motifs_dict
    profiles = profiles[profiles.index.isin(drep_motifs_dict.keys())]
    print ("filtered shape 2", profiles.shape)

    # Filter columns where any value is greater than 0.5
    # profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]
    # print ("filtered shape", profiles.shape)

    df = profiles.T
    ctg_phylum = get_phylum(best_ctgs)

    phylum_colors = sns.color_palette("husl", len(set(ctg_phylum.values())))
    phylum_color_dict = {phylum: color for phylum, color in zip(set(ctg_phylum.values()), phylum_colors)}
    phylum_colors_list = [phylum_color_dict[ctg_phylum[ctg]] for ctg in df.index]

    # Use clustermap instead of heatmap
    g = sns.clustermap(
        df,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(30, 15),
        row_colors=phylum_colors_list,
        xticklabels=True,
        yticklabels=True
    )

    # Add legend for phylum colors
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=phylum_color_dict[phylum], label=phylum) for phylum in phylum_color_dict]
    plt.legend(handles=handles, title="Phylum", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


    plt.savefig("../../tmp/results/heatmap.png", dpi=300, bbox_inches='tight')
    plt.clf()

    # plot_drep_similarity_vs_cos_sim(df)

def plot_drep_similarity_vs_cos_sim(df):
    ## compute pairwise cosine similarity matrix
    cosine_sim_matrix = pd.DataFrame(
        cosine_similarity(df),
        index=df.index,
        columns=df.index
    )
    drep_sim_dict = get_ctg_sim()
    ## for each pair of contigs, get the Euclidean distance and dRep similarity
    drep_sim_list = []
    for ctg1 in df.index:
        for ctg2 in df.index:
            if ctg1 != ctg2:
                euclidean_dist = cosine_sim_matrix.loc[ctg1, ctg2]
                if ctg1 + "&" + ctg2 not in drep_sim_dict:
                    print (f"Not Found dRep similarity for {ctg1} and {ctg2}")
                    continue
                drep_sim = drep_sim_dict[ctg1 + "&" + ctg2]
                drep_sim_list.append((ctg1, ctg2, euclidean_dist, drep_sim))
    drep_sim_df = pd.DataFrame(drep_sim_list, columns=['ctg1', 'ctg2', 'cos_sim', 'drep_sim'])
    ### plot scatter plot of dRep similarity vs Euclidean distance
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=drep_sim_df, x='cos_sim', y='drep_sim', alpha=0.3, s=25)
    plt.xlabel('Cosine Similarity')
    plt.ylabel('dRep Similarity')
    plt.savefig("../../tmp/results/drep_similarity_vs_cos_sim.png", dpi=300, bbox_inches='tight')
    # # print (drep_sim_dict)

def get_ctg_sim():
    drep_sim_dict = {}
    drep_sim = "/home/shuaiw/borg/contigs/dRep/dRep_out4/data_tables/Mdb.csv"
    df = pd.read_csv(drep_sim)
    for index, row in df.iterrows():
        ctg1 = row['genome1'][:-6]
        ctg2 = row['genome2'][:-6]
        drep_sim_dict[ctg1 + "&" + ctg2] = row['similarity']
        drep_sim_dict[ctg2 + "&" + ctg1] = row['similarity']
    return drep_sim_dict

def get_ctg_sim2():
    drep_sim_dict = {}
    drep_sim = "/home/shuaiw/borg/contigs/fastANI.output"
    df = pd.read_csv(drep_sim, header =None, sep="\t")
    for index, row in df.iterrows():
        ## get first colunm as ctg1 and second column as ctg2

        ctg1 = row[0][:-6].split("/")[-1]
        ctg2 = row[1][:-6].split("/")[-1]
        drep_sim_dict[ctg1 + "&" + ctg2] = row[2]
        drep_sim_dict[ctg2 + "&" + ctg1] = row[2]
    # print (drep_sim_dict)
    return drep_sim_dict

def read_drep_cluster(drep_clu_file, best_ctgs):
    
    drep_clu_dict = {}
    df = pd.read_csv(drep_clu_file)
    for index, row in df.iterrows():
        drep_clu_dict[row['genome'][:-6]] = row['secondary_cluster']
    ## count the number of contigs in each cluster
    clu_count = {}
    for ctg, clu in drep_clu_dict.items():
        if ctg not in best_ctgs:
            continue
        if clu not in clu_count:
            clu_count[clu] = 0
        clu_count[clu] += 1
    multiple_strain_ctg = {}
    ctg_phylum = {}
    clu_set = set()
    for ctg, clu in drep_clu_dict.items():
        if clu not in clu_count:
            continue
        if clu_count[clu] > 1:
            multiple_strain_ctg[ctg] = clu
            ctg_phylum[ctg] = clu
            clu_set.add(clu)
    print ("multiple_strain_ctg", multiple_strain_ctg)
    return multiple_strain_ctg, ctg_phylum, len(clu_set)

def jaccard():
    # drep_clu_file = "/home/shuaiw/methylation/data/borg/contigs/dRep/dRep_out/data_tables/Cdb.csv"
    # profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"
    # depth_file = "/home/shuaiw/borg/bench/soil/run1/mean_depth.csv"

    prefix = "ERR12723529"
    work_dir = "/home/shuaiw/borg/paper/run"
    drep_clu_file = os.path.join(work_dir, prefix, "dRep_out", "data_tables", "Cdb.csv")
    profile_file = os.path.join(work_dir, prefix, f"{prefix}_methylation", "motif_profile.csv")
    depth_file = os.path.join(work_dir, prefix, f"{prefix}_methylation", "mean_depth.csv")

    multiple_strain_ctg, ctg_phylum = read_drep_cluster(drep_clu_file)
    ## read motif profile file
    profiles = pd.read_csv(profile_file, index_col=0)
    ## binary it by 0.5
    profiles_bin = (profiles >= 0.5).astype(int)
    profiles_bin = profiles_bin.T
    ## calculate Jaccard similarity in the same cluster and different clusters
    jaccard_similarities = []
    ctg_list = list(multiple_strain_ctg.keys())
    # print (profiles_bin)
    for i in range(len(ctg_list)):
        for j in range(i + 1, len(ctg_list)):
            ctg1 = ctg_list[i]
            ctg2 = ctg_list[j]
            ## skip if ctg1 not in profiles_bin or ctg2 not in profiles_bin
            if ctg1 not in profiles_bin.index or ctg2 not in profiles_bin.index:
                continue
            if ctg1 == ctg2:
                continue
            # print (profiles_bin.loc[ctg1])
            intersection = np.sum((profiles_bin.loc[ctg1] & profiles_bin.loc[ctg2]).astype(int))
            union = np.sum((profiles_bin.loc[ctg1] | profiles_bin.loc[ctg2]).astype(int))
            if union == 0:
                continue
            jaccard_similarity = intersection / union if union > 0 else 0
            if multiple_strain_ctg[ctg1] != multiple_strain_ctg[ctg2]:
                # Calculate Jaccard similarity
                jaccard_similarities.append([jaccard_similarity, "diff"])
            else:
                jaccard_similarities.append([jaccard_similarity, "same"])
    print (f"Total {len(jaccard_similarities)} Jaccard similarities calculated.")
    ## plot the Jaccard similarities
    df_jaccard = pd.DataFrame(jaccard_similarities, columns=['jaccard_similarity', 'cluster_type'])
    plt.figure(figsize=(7, 6))
    sns.boxplot(data=df_jaccard, x='cluster_type', y='jaccard_similarity', palette='Set2')
    plt.title('Jaccard Similarity of Motifs in Multiple Strain Clusters')
    plt.xlabel('Cluster Type')
    plt.ylabel('Jaccard Similarity')
    plt.savefig("../../tmp/results/jaccard_similarity.png", dpi=300, bbox_inches='tight')
    plt.clf()   

def jaccard_batch():
    all_dir = "/home/shuaiw/borg/paper/run2/"
    jaccard_similarities = []
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        print (f"Processing {prefix}...")
        # work_dir = f"{all_dir}/{prefix}/{prefix}_methylation3"

        # prefix = "ERR12723529"
        work_dir = "/home/shuaiw/borg/paper/run2"
        # drep_clu_file = os.path.join(work_dir, prefix, "dRep_out_99", "data_tables", "Cdb.csv")
        # drep_clu_file = os.path.join(work_dir, prefix, "dRep_out", "data_tables", "Cdb.csv")
        drep_clu_file = os.path.join(work_dir, prefix, "dRep_out_97", "data_tables", "Cdb.csv")
        profile_file = os.path.join(work_dir, prefix, f"{prefix}_methylation3", "motif_profile.csv")
        depth_file = os.path.join(work_dir, prefix, f"{prefix}_methylation3", "mean_depth.csv")
        fai = os.path.join(work_dir, prefix, f"{prefix}.hifiasm.p_ctg.rename.fa.fai")
        all_host_file = f"{all_dir}/{prefix}/all_host_ctgs.tsv" 
        ## skip if all_host_file does not exist
        if not os.path.exists(all_host_file):
            print(f"Skipping {prefix} as all_host_file does not exist.")
            continue
        if not os.path.exists(fai):
            continue
        best_ctgs, good_depth = get_best_ctg(depth_file, fai)
        best_ctgs = read_all_host(all_host_file, good_depth)
        ## read length of each ctg in a dict

        ## skip if not exist
        if not os.path.exists(drep_clu_file) or not os.path.exists(profile_file):
            print(f"Skipping {prefix} due to missing files.")
            continue

        multiple_strain_ctg, ctg_phylum, clu_num = read_drep_cluster(drep_clu_file, best_ctgs)
        if clu_num < 2:
            continue
        ## read motif profile file
        profiles = pd.read_csv(profile_file, index_col=0)
        ## binary it by 0.5
        profiles_bin = (profiles >= 0.4).astype(int)
        profiles_bin = profiles_bin.T
        
        # ctg_list = list(multiple_strain_ctg.keys())
        ctg_list = list(ctg_phylum.keys())
        # print (profiles_bin)
        for i in range(len(ctg_list)):
            for j in range(i + 1, len(ctg_list)):
                ctg1 = ctg_list[i]
                ctg2 = ctg_list[j]
                ## skip if ctg1 not in profiles_bin or ctg2 not in profiles_bin
                if ctg1 not in profiles_bin.index or ctg2 not in profiles_bin.index:
                    continue
                if ctg1 == ctg2:
                    continue
                # print (profiles_bin.loc[ctg1])
                intersection = np.sum((profiles_bin.loc[ctg1] & profiles_bin.loc[ctg2]).astype(int))
                union = np.sum((profiles_bin.loc[ctg1] | profiles_bin.loc[ctg2]).astype(int))
                if union == 0:
                    continue
                jaccard_similarity = intersection / union if union > 0 else 0
                if ctg1 in multiple_strain_ctg:
                    cluster1 = multiple_strain_ctg[ctg1]
                else:
                    cluster1 = ctg1
                if ctg2 in multiple_strain_ctg:
                    cluster2 = multiple_strain_ctg[ctg2]
                else:
                    cluster2 = ctg2

                if cluster1 != cluster2:
                    # Calculate Jaccard similarity
                    jaccard_similarities.append([jaccard_similarity, "diff", prefix])
                else:
                    jaccard_similarities.append([jaccard_similarity, "same", prefix])
        # if len(jaccard_similarities) > 5:
        #     break
    print (f"Total {len(jaccard_similarities)} Jaccard similarities calculated.")
    ## plot the Jaccard similarities
    df_jaccard = pd.DataFrame(jaccard_similarities, columns=['jaccard_similarity', 'cluster_type', 'sample'])
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df_jaccard, hue='cluster_type', x='sample', y='jaccard_similarity', palette='Set2')
    plt.xticks(rotation=90)
    plt.title('Jaccard Similarity of Motifs in Multiple Strain Clusters')
    plt.ylabel('Jaccard Similarity')
    plt.savefig("../../tmp/results/jaccard_similarity_batch.pdf", dpi=300, bbox_inches='tight')
    plt.clf()   
    
def get_phylum(best_ctgs):
    import re
    bin_phylum = {}
    df = pd.read_csv(gtdb_file, sep="\t")
    for index, row in df.iterrows():
        if re.search("Unclassified", row["classification"]):
            taxon = "Unclassified"
        else:
            # print (row["classification"])
            taxon = row["classification"].split(";")[1].strip()
        bin_phylum[row["user_genome"]] = taxon
    ctg_bin = {}
    with open(bin_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            ctg, bin_id = line.strip().split("\t")
            ctg_bin[ctg] = bin_id
    
    ctg_phylum = {}
    for ctg in best_ctgs:
        bin_id = ctg_bin[ctg]
        if bin_id in bin_phylum:
            ctg_phylum[ctg] = bin_phylum[bin_id]
        else:
            ctg_phylum[ctg] = "Unclassified"
    return (ctg_phylum)

def out_best_ctgs(ref, best_ref, best_ctgs):
    ## use biopython to extract the best contigs from the reference fasta file
    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    with open(ref, "r") as f_in, open(best_ref, "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            if record.id in best_ctgs:
                SeqIO.write(record, f_out, "fasta")
    print(f"Extracted {len(best_ctgs)} best contigs to {best_ref}")

def out_best_ctgs2(ref, best_ref, best_ctgs, best_ctg_dir):
    ## use biopython to extract the best contigs from the reference fasta file
    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser    
    with open(ref, "r") as f_in:
        for record in SeqIO.parse(f_in, "fasta"):
            if record.id in best_ctgs:
                with open(os.path.join(best_ctg_dir, str(record.id) + ".fasta"), "w") as f_out:
                    SeqIO.write(record, f_out, "fasta")
    print(f"Extracted {len(best_ctgs)} best contigs to {best_ref}")

def count_motifs(best_ctgs, all_dir, prefix, environment, ctg_taxa_dict):
    data = []
    has_motif_ctg_num = 0
    motif_num_list = []
    for ctg in best_ctgs:
        ctg_obj = My_contig(prefix, all_dir, ctg)
        df_motif = ctg_obj.read_motif(min_frac=0.3, min_sites=100)
        # print (df_motif)
        ## only keep themotifs with fraction >= 0.4, and nDetected >=100
        unique_motifs_num, unique_motifs, unique_motifs_identifier = get_unique_motifs(df_motif)
        if unique_motifs_num > 0:
            has_motif_ctg_num += 1
        motif_num_list.append(unique_motifs_num)
        ctg_lineage = ctg_taxa_dict[ctg] if ctg in ctg_taxa_dict else "Unknown"
        ctg_phylum = classify_taxa(ctg_lineage, level='phylum')
        ctg_domain = classify_taxa(ctg_lineage, level='domain')
        ctg_len = ctg_obj.get_ctg_len()
        ctg_RM_dict = ctg_obj.load_RM()
        # print (ctg_RM_dict)
        ctg_RM_num = len(ctg_RM_dict)
        # if ctg_phylum == "Unknown":
        #     os.system(f"cp {ctg_obj.ctg_ref} /home/shuaiw/borg/paper/specificity/ming_second_chr")
        data.append([prefix, unique_motifs_num, environment,ctg, ctg_phylum, ctg_domain, ctg_lineage, ctg_len,ctg_RM_num])

    return data

def count_MT_num(anno_file):
    
    df = pd.read_csv(anno_file, sep="\t")
    # Extract contig name from 'Gene' column and add as a new column
    df['contig'] = df['Gene'].apply(lambda x: '_'.join(str(x).split('_')[:-1]))  
    ctg_MT_num = {}
    for index, row in df.iterrows():
        if row['contig'] not in ctg_MT_num:
            ctg_MT_num[row['contig']] = set()
        if row['Gene type'] == 'MT':
            ctg_MT_num[row['contig']].add(row['HMM'])
    return ctg_MT_num

def count_modified_base(all_dir, prefix, best_ctgs, length_dict, env, ctg_taxa_dict,score_cutoff = 30):
    base_data = []
    for ctg in best_ctgs:
        length = length_dict[ctg]
        ctg_obj = My_contig(prefix, all_dir, ctg)
        ctg_obj.get_mod_ratio(score_cutoff)
        ctg_lineage = ctg_taxa_dict[ctg_obj.contig] if ctg_obj.contig in ctg_taxa_dict else "Unknown"
        ctg_phylum = classify_taxa(ctg_lineage, level='phylum')
        base_data.append([prefix, ctg, length, ctg_obj.modified_num, 
                          ctg_obj.modified_motif_num, ctg_obj.modified_ratio, 
                          ctg_obj.modified_motif_ratio, ctg_obj.motif_ratio, env, ctg_phylum])
    return base_data

def get_stastics():
    df_genome_data = pd.read_csv("../../tmp/results2/genome_data_all_samples.csv")
    df_all_base_data = pd.read_csv("../../tmp/results2/base_count_all_samples.csv")
    df_all_data = pd.read_csv("../../tmp/results2/motif_num_all_samples.csv")
    ## exclude 96plex environment for df_genome_data, df_all_base_data, df_all_data
    df_genome_data = df_genome_data[df_genome_data['environment'] != '96plex']
    df_all_base_data = df_all_base_data[df_all_base_data['environment'] != '96plex']
    df_all_data = df_all_data[df_all_data['environment'] != '96plex']

    ## count sample number for each environment in df_genome_data
    env_groups = df_genome_data.groupby('environment')
    print ("Sample number for each environment:")
    for env, group in env_groups:
        total_samples_env = group.shape[0]
        print (f"Environment: {env}, Total samples: {total_samples_env}")
    print ("*************************")

    print ("Genome data summary:")
    ## count how many samples proportion have map_ratio > 0.5 and samples proportion have map_ratio > 0.7
    total_samples = df_genome_data.shape[0]
    map_ratio_50 = df_genome_data[df_genome_data['map_ratio'] > 0.5].shape[0]
    map_ratio_70 = df_genome_data[df_genome_data['map_ratio'] > 0.7].shape[0]
    print (f"Total samples: {total_samples}, Samples with map_ratio > 0.5: {map_ratio_50} ({map_ratio_50/total_samples*100:.2f}%), Samples with map_ratio > 0.7: {map_ratio_70} ({map_ratio_70/total_samples*100:.2f}%)")

    ## get proportion with motif_num > 0 in df_all_data
    total_ctgs = df_all_data.shape[0]
    ctgs_with_motif = df_all_data[df_all_data['motif_num'] > 0].shape[0]
    print (f"Total contigs: {total_ctgs}, Contigs with motif_num > 0: {ctgs_with_motif} ({ctgs_with_motif/total_ctgs*100:.2f}%)")
    ## count the average motif_num for all contigs
    avg_motif_num = df_all_data['motif_num'].mean()
    print (f"Average motif_num for all contigs: {avg_motif_num:.2f}") 
    ## count the proportion with motif_num > 0 for each environment
    env_groups = df_all_data.groupby('environment')
    for env, group in env_groups:
        total_ctgs_env = group.shape[0]
        ctgs_with_motif_env = group[group['motif_num'] > 0].shape[0]
        ## also print the average motif_num for each environment
        avg_motif_num_env = group['motif_num'].mean()
        print (f"Environment: {env}, Total contigs: {total_ctgs_env}, Contigs with motif_num > 0: {ctgs_with_motif_env}({ctgs_with_motif_env/total_ctgs_env*100:.2f}%), Average motif_num: {avg_motif_num_env:.2f}")
    print ("*************************")
    ## count average regulate_motif_num for each environment
    env_groups = df_genome_data.groupby('environment')
    for env, group in env_groups:
        total_samples_env = group.shape[0]
        avg_regulate_motif_num_env = group['regulate_motif_num'].mean()
        print (f"Environment: {env}, Total samples: {total_samples_env}, Average regulate_motif_num: {avg_regulate_motif_num_env:.2f}")
    ## count the proportion of infant sample has regulate_motif_num > 0
    infant_group = df_genome_data[df_genome_data['environment'] == 'infant']
    total_infant_samples = infant_group.shape[0]
    infant_samples_with_regulate_motif = infant_group[infant_group['regulate_motif_num'] > 0].shape[0]
    print (f"Infant samples: {total_infant_samples}, Infant samples with regulate_motif_num > 0: {infant_samples_with_regulate_motif} ({infant_samples_with_regulate_motif/total_infant_samples*100:.2f}%)")
    print ("*************************")
    ### count the average  modified_ratio for each environment
    env_groups = df_all_base_data.groupby('environment')
    for env, group in env_groups:
        total_ctgs_env = group.shape[0]
        avg_modified_ratio_env = group['modified_ratio'].mean()
        avg_modified_motif_ratio_env = group['modified_motif_ratio'].mean()
        avg_motif_ratio_env = group['motif_ratio'].mean()
        print (f"Environment: {env}, Total contigs: {total_ctgs_env}, Average modified_ratio: {avg_modified_ratio_env:.4f}, Average modified_motif_ratio: {avg_modified_motif_ratio_env:.4f}, Average motif_ratio: {avg_motif_ratio_env:.4f}")
    print ("*************************")

def plot_motif(df_all_data, fig_dir):
    ### sort by sample
    df_all_data = df_all_data.sort_values(by='sample')
    ## plot boxplot where sample is on x-axis and motif_num is on y-axis
    plt.figure(figsize=(14, 6))
    sns.boxplot(data=df_all_data, x='sample', y='motif_num', hue='environment', palette='Set2')
    plt.xticks(rotation=90)
    plt.title('Distribution of Motif Numbers Across Samples')
    plt.xlabel('Sample')
    plt.ylabel('Number of Motifs')
    plt.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f"{fig_dir}/motif_num_distribution_all_samples.png", dpi=300, bbox_inches='tight')

def plot_motif_env(df_all_data, fig_dir):
    ctg_num_cutoff = 10
    ### sort by sample
    df_all_data = df_all_data.sort_values(by='environment')
    ## plot boxplot where sample is on x-axis and motif_num is on y-axis

    plt.figure(figsize=(6, 6))
    env_counts = df_all_data['environment'].value_counts()
    env_with_enough_data = env_counts[env_counts > ctg_num_cutoff].index
    df_filtered_env = df_all_data[df_all_data['environment'].isin(env_with_enough_data)]

    order = df_filtered_env.groupby('environment')['motif_num'].mean().sort_values().index
    sns.boxplot(data=df_filtered_env, x='environment', y='motif_num', order=order)

    plt.xticks(rotation=90)
    plt.xlabel(f'Environment (>{ctg_num_cutoff} contigs)')
    plt.ylabel('Number of Motifs')
    plt.savefig(f"{fig_dir}/motif_num_env_filtered.png", dpi=300, bbox_inches='tight')
    
    ## Filter phyla with > 50 values for better statistical validity
    phylum_counts = df_all_data['phylum'].value_counts()
    phylum_with_enough_data = phylum_counts[phylum_counts > ctg_num_cutoff].index
    df_filtered_phylum = df_all_data[df_all_data['phylum'].isin(phylum_with_enough_data)]
    
    ## also plot one with phylum is on x-axis and motif_num is on y-axis
    plt.figure(figsize=(10, 5))
    order = df_filtered_phylum.groupby('phylum')['motif_num'].mean().sort_values().index
    sns.boxplot(data=df_filtered_phylum, x='phylum', y='motif_num', order=order)
    plt.xticks(rotation=90)
    plt.xlabel(f'Phylum (>{ctg_num_cutoff} contigs)')
    plt.ylabel('Number of Motifs')
    plt.savefig(f"{fig_dir}/motif_num_phylum_filtered.png", dpi=300, bbox_inches='tight')

    ## Filter environments with > 50 values for better statistical validity

    
    ## also plot one with environment as x-axis and motif_num is on y-axis, and domain as hue
    plt.figure(figsize=(8, 3))
    order = df_filtered_env.groupby('environment')['motif_num'].mean().sort_values().index
    sns.boxplot(data=df_filtered_env, x='environment', y='motif_num', hue='domain', order=order)
    plt.xticks(rotation=90)
    plt.xlabel(f'Environment (>{ctg_num_cutoff} samples)')
    plt.ylabel('Number of Motifs')
    plt.legend(title='Domain', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(f"{fig_dir}/motif_num_env_domain_filtered.png", dpi=300, bbox_inches='tight')

def sort_genome_motif(df_all_data, fig_dir):
    ## add column species to df_all_data based on the lineage column, get the last word after semicolon
    df_all_data['species'] = df_all_data['lineage'].apply(get_detail_taxa_name)


    top20_df = df_all_data.sort_values(by="motif_num", ascending=False).head(20)
    
    # Create labels that include cluster name and member count
    top20_df['combined_label'] = top20_df.apply(lambda row: f"{row['species']}\n({row['contig']}, n={row['motif_num']})", axis=1)

    plt.figure(figsize=(8,8))
    ax = sns.barplot(x="motif_num", y="combined_label", data=top20_df, palette="viridis")
    
    # Adjust layout to accommodate longer labels
    plt.subplots_adjust(left=0.4)  # Make more space for y-axis labels
    
    # Set font size for better readability
    plt.tick_params(axis='y', labelsize=10)
    plt.tick_params(axis='x', labelsize=10)
    
    plt.xlabel("Number of Motifs", fontsize=12)
    plt.ylabel("Taxa (contig, motif number)", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/top20_genome_motif.pdf", bbox_inches='tight')


def plot_base(df_all_base_data, fig_dir):
    ### sort by sample
    ctg_num_cutoff = 20
    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    df_all_base_data = df_all_base_data.sort_values(by='sample')
    
    # Filter data for environment and phylum
    env_counts = df_all_base_data['environment'].value_counts()
    env_with_enough_data = env_counts[env_counts > ctg_num_cutoff].index
    df_filtered_env = df_all_base_data[df_all_base_data['environment'].isin(env_with_enough_data)]
    
    phylum_counts = df_all_base_data['phylum'].value_counts()
    phylum_with_enough_data = phylum_counts[phylum_counts > ctg_num_cutoff].index
    df_filtered_phylum = df_all_base_data[df_all_base_data['phylum'].isin(phylum_with_enough_data)]
    
    ## Row 1: modified_ratio
    # Environment
    order1_env = df_filtered_env.groupby('environment')['modified_ratio'].mean().sort_values().index
    sns.boxplot(data=df_filtered_env, x='environment', y='modified_ratio', palette='Set2', ax=axes[0, 0], order=order1_env)
    axes[0, 0].set_title('Modified Ratios by Environment')
    axes[0, 0].set_xlabel(f'Environment (>{ctg_num_cutoff} contigs)')
    axes[0, 0].set_ylabel('Modified Ratio')
    axes[0, 0].tick_params(axis='x', rotation=90)
    
    # Phylum
    order1_phylum = df_filtered_phylum.groupby('phylum')['modified_ratio'].mean().sort_values().index
    sns.boxplot(data=df_filtered_phylum, x='phylum', y='modified_ratio', palette='Set2', ax=axes[0, 1], order=order1_phylum)
    axes[0, 1].set_title('Modified Ratios by Phylum')
    axes[0, 1].set_xlabel(f'Phylum (>{ctg_num_cutoff} contigs)')
    axes[0, 1].set_ylabel('Modified Ratio')
    axes[0, 1].tick_params(axis='x', rotation=90)

    ## Row 2: modified_motif_ratio
    # Environment
    order2_env = df_filtered_env.groupby('environment')['modified_motif_ratio'].mean().sort_values().index
    sns.boxplot(data=df_filtered_env, x='environment', y='modified_motif_ratio', palette='Set2', ax=axes[1, 0], order=order2_env)
    axes[1, 0].set_title('Modified Motif Ratios by Environment')
    axes[1, 0].set_xlabel(f'Environment (>{ctg_num_cutoff} contigs)')
    axes[1, 0].set_ylabel('Modified Motif Ratio')
    axes[1, 0].tick_params(axis='x', rotation=90)
    
    # Phylum
    order2_phylum = df_filtered_phylum.groupby('phylum')['modified_motif_ratio'].mean().sort_values().index
    sns.boxplot(data=df_filtered_phylum, x='phylum', y='modified_motif_ratio', palette='Set2', ax=axes[1, 1], order=order2_phylum)
    axes[1, 1].set_title('Modified Motif Ratios by Phylum')
    axes[1, 1].set_xlabel(f'Phylum (>{ctg_num_cutoff} contigs)')
    axes[1, 1].set_ylabel('Modified Motif Ratio')
    axes[1, 1].tick_params(axis='x', rotation=90)

    ## Row 3: motif_ratio
    # Environment
    order3_env = df_filtered_env.groupby('environment')['motif_ratio'].mean().sort_values().index
    sns.boxplot(data=df_filtered_env, x='environment', y='motif_ratio', palette='Set2', ax=axes[2, 0], order=order3_env)
    axes[2, 0].set_title('Motif Ratios by Environment')
    axes[2, 0].set_xlabel(f'Environment (>{ctg_num_cutoff} contigs)')
    axes[2, 0].set_ylabel('Motif Ratio')
    axes[2, 0].tick_params(axis='x', rotation=90)
    
    # Phylum
    order3_phylum = df_filtered_phylum.groupby('phylum')['motif_ratio'].mean().sort_values().index
    sns.boxplot(data=df_filtered_phylum, x='phylum', y='motif_ratio', palette='Set2', ax=axes[2, 1], order=order3_phylum)
    axes[2, 1].set_title('Motif Ratios by Phylum')
    axes[2, 1].set_xlabel(f'Phylum (>{ctg_num_cutoff} contigs)')
    axes[2, 1].set_ylabel('Motif Ratio')
    axes[2, 1].tick_params(axis='x', rotation=90)

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/base_count_all_samples.png", dpi=300, bbox_inches='tight')
    ### save the dataframe

def plot_meta(df_genome_data, fig_dir):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 12))
    df_genome_data = df_genome_data.sort_values(by='sample')
    df_genome_data.set_index('sample', inplace=True)

    sns.barplot(x=df_genome_data.index, y='map_ratio', hue='environment', data=df_genome_data.reset_index(), ax=ax1, palette='Set2')
    ax1.set_title('Mapping Ratio Across Samples')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Mapping Ratio ')
    ax1.tick_params(axis='x', rotation=90)

    sns.barplot(x=df_genome_data.index, y='linkage_num', hue='environment', data=df_genome_data.reset_index(), ax=ax2, palette='Set2')
    ax2.set_title('Linkage Number Across Samples')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Linkage Number ')
    ax2.tick_params(axis='x', rotation=90)

    sns.barplot(x=df_genome_data.index, y='regulate_motif_num', hue='environment', data=df_genome_data.reset_index(), ax=ax3, palette='Set2')
    ax3.set_title('Regulatory Motif Number Across Samples')
    ax3.set_xlabel('Sample')
    ax3.set_ylabel('Regulatory Motif Number ')
    ax3.tick_params(axis='x', rotation=90)

    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    # Remove legend for ax2 and ax3
    ax2.legend_.remove()
    ax3.legend_.remove()
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/meta_all_samples.png", dpi=300, bbox_inches='tight')

def plot_genome(df_genome_data, fig_dir):
    ### plot bar plot for genome size and N50 in separate subplots, using log scale for y-axis
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 12))
    df_genome_data = df_genome_data.sort_values(by='sample')
    df_genome_data.set_index('sample', inplace=True)
    
    # N50 subplot
    sns.barplot(x=df_genome_data.index, y='N50', hue='environment', data=df_genome_data.reset_index(), ax=ax1, palette='Set2')
    ax1.set_yscale('log')
    ax1.set_title('N50 Across Samples')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('N50 Size (log scale)')
    ax1.tick_params(axis='x', rotation=90)
    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Genome size subplot
    sns.barplot(x=df_genome_data.index, y='genome_size', hue='environment', data=df_genome_data.reset_index(), ax=ax2, palette='Set2')
    ax2.set_yscale('log')
    ax2.set_title('Genome Size Across Samples')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Genome Size (log scale)')
    ax2.tick_params(axis='x', rotation=90)
    ax2.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    sns.barplot(x=df_genome_data.index, y='best_ctg_num', hue='environment', data=df_genome_data.reset_index(), ax=ax3, palette='Set2')
    ax3.set_title('Best Contig Number Across Samples')
    ax3.set_xlabel('Sample')
    ax3.set_ylabel('Best Contig Number')
    ax3.tick_params(axis='x', rotation=90)
    ax3.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/genome_size_N50_all_samples.png", dpi=300, bbox_inches='tight')

def plot_MTase(df_all_data, fig_dir):
    ##scatter plot and box plot for MTase num vs motif num
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Scatter plot with regression line
    sns.scatterplot(data=df_all_data, x='RM_num', y='motif_num', hue='environment', 
                   palette='Set2', alpha=0.6, s=40, ax=ax1)
    # Add regression line
    sns.regplot(data=df_all_data, x='RM_num', y='motif_num', scatter=False, 
                color='red', ax=ax1)
    ax1.set_title('MTase Number vs Motif Number (Scatter)')
    ax1.set_xlabel('MTase Number')
    ax1.set_ylabel('Motif Number')
    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Box plot - each x shows a number of MTase, y shows number of motifs
    # Sort the RM_num values and create ordered categories
    rm_order = sorted(df_all_data['RM_num'].unique())
    df_all_data['RM_num_str'] = df_all_data['RM_num'].astype(str)
    sns.boxplot(data=df_all_data, x='RM_num_str', y='motif_num',
                order=[str(x) for x in rm_order], ax=ax2)
    ax2.set_title('Motif Number Distribution by MTase Count')
    ax2.set_xlabel('MTase Number')
    ax2.set_ylabel('Motif Number')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/MTase_vs_motif_num.png", dpi=300, bbox_inches='tight')

    ## print all the rows with motif_num > RM_num
    # df_motif_more = df_all_data[df_all_data['motif_num'] > df_all_data['RM_num']]
    ## output all the rows with motif_num > RM_num
    ## sort column order, first sample then motif_num and RM_num
    # df_motif_more = df_motif_more[['sample', 'motif_num', 'RM_num', 'environment', 'contig', 'phylum', 'domain', 'lineage','ctg_len']]
    # df_motif_more.to_csv(f"{fig_dir}/motif_num_greater_than_MTase_num.csv", index=False)

def read_metadata(meta_file):
    sample_env_dict = {}
    with open(meta_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            sample = parts[1]
            env = parts[3]
            sample_env_dict[sample] = env
    return sample_env_dict

def main(all_dir, fig_dir, sample_env_dict):
    all_data = []
    all_base_data = []
    genome_data = []
    genome_list = []
    ctg_taxa_dict = get_ctg_taxa(all_dir)
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if re.search("sludge", prefix) or re.search("human", prefix):
            continue
        sample_obj = My_sample(prefix, all_dir)
        if not os.path.exists(sample_obj.depth_file):
            print (f"Skipping {prefix} as depth file not found.")
            continue
        # if prefix != "ocean_1":
        #     continue
        
        print (f"Processing {prefix}...")
        sample_obj.read_depth()
        map_ratio = sample_obj.read_mapping()
        linkage_num =  sample_obj.read_host()
        regulate_motif_num = sample_obj.read_orphan()
        N50, genome_size = sample_obj.get_N50_size()
        cpu_time_hours, wallclock_time_hours, peak_memory_gb = sample_obj.get_time()
        print(f"{prefix}: N50 size: {N50}, Genome size: {genome_size}")
        my_genome_list, best_ctgs = sample_obj.get_final_best_ctg2()
        # analyze_gene(all_dir, meta_dir, prefix, best_ctgs)
        # print ("best ctgs num:", len(best_ctgs))
        genome_list += my_genome_list
        print (f">>>Total {len(genome_list)} high depth contigs collected so far.")

        genome_data.append([prefix, N50, genome_size, sample_env_dict[prefix], map_ratio, linkage_num, \
                            regulate_motif_num, len(best_ctgs), cpu_time_hours, wallclock_time_hours, peak_memory_gb])

        # print (f"Total {len(best_ctgs)} best contigs with depth >= 10 found.")
        all_data += count_motifs(best_ctgs, all_dir, prefix, sample_env_dict[prefix], ctg_taxa_dict)
        # all_base_data += count_modified_base(all_dir, prefix, best_ctgs, sample_obj.length_dict, sample_env_dict[prefix],ctg_taxa_dict)
        # break
    print ("start plot...")
    df_all_data = pd.DataFrame(all_data, columns=['sample', 'motif_num', 'environment', 'contig', 'phylum', 'domain', 'lineage','ctg_len', 'RM_num'])
    df_genome_data = pd.DataFrame(genome_data, columns=['sample', 'N50', 'genome_size', 'environment', 'map_ratio', 'linkage_num', 'regulate_motif_num','best_ctg_num','cpu_time_hours','wallclock_time_hours','peak_memory_gb'])
    # df_all_base_data = pd.DataFrame(all_base_data, columns=['sample', 'ctg', 'length', 'modified_num', 'modified_motif_num', 'modified_ratio', 'modified_motif_ratio', 'motif_ratio', 'environment','phylum'])

    df_genome_data.to_csv(f"{fig_dir}/genome_data_all_samples.csv", index = False)
    # df_all_base_data.to_csv(f"{fig_dir}/base_count_all_samples.csv", index=False)
    df_all_data.to_csv(f"{fig_dir}/motif_num_all_samples.csv", index=False)    
    
    # df_all_data = pd.read_csv(f"{fig_dir}/motif_num_all_samples.csv")
    # df_genome_data = pd.read_csv(f"{fig_dir}/genome_data_all_samples.csv")
    # df_all_base_data = pd.read_csv(f"{fig_dir}/base_count_all_samples.csv")

    # plot_motif_env(df_all_data, fig_dir)
    # sort_genome_motif(df_all_data, fig_dir)
    # plot_motif(df_all_data, fig_dir)
    # plot_genome(df_genome_data, fig_dir)
    # plot_meta(df_genome_data, fig_dir)
    # # plot_base(df_all_base_data, fig_dir)
    # plot_MTase(df_all_data, fig_dir)
    # ## save genome data

    # with open(genome_list_file, "w") as f:
    #     for genome in genome_list:
    #         f.write(genome + "\n")

def rerun(fig_dir):
    df_all_data = pd.read_csv(f"{fig_dir}/motif_num_all_samples.csv")
    plot_motif_env(df_all_data, fig_dir)

def analyze_gene(all_dir, meta_dir, prefix, contig_list):
    
    for contig in contig_list:
        

        contig_obj = My_contig(prefix, all_dir, contig)
        cds_gff = meta_dir + "/" + contig + "/prokka/" + f"{contig}.gff"
        modified_gff = contig_obj.reprocess_gff
        genome_file = contig_obj.ctg_ref
        count_dir = f"{meta_dir}/{contig}/count/"
        count_file = os.path.join(count_dir, f"{contig}_region_count.csv")

        
        # Create count directory if needed
        os.makedirs(count_dir, exist_ok=True)
        print (contig)
        # Initialize gene object
        my_gene = My_gene(cds_gff, contig, genome_file)
        my_gene.collect_regulation_region()
        my_gene.read_modified_gff(modified_gff)
        
        # Perform intersection     
        df = my_gene.intersect_gene()
        print (df)

        # df.to_csv(count_file, index=False)
        # processed += 1
            


def plot_coding( meta_dir, fig_dir, min_len = 1000000):
    ctg_num_cutoff = 10
    whole_df = pd.DataFrame()
    for count_csv in glob.glob(f"{meta_dir}/*/count/*_region_count.csv"):
        df = pd.read_csv(count_csv)
        whole_df = pd.concat([whole_df, df], ignore_index=True)
    ## filter genome_length >= min_len
    whole_df = whole_df[whole_df['genome_length'] >= min_len]
    print ("Total {0} genomes with length >= {1} bp.".format(len(whole_df), min_len))
    ## the data is like:
    ##  genome  genome_length  regulatory_count  regulatory_length  regulatory_frequency  cds_count  ...  non_coding_count  non_coding_length  non_coding_frequency            sample  environment           phylum
    ## plot boxplot, y is value, x is environment, hue is region_type
    melted_df_raw = pd.melt(whole_df, id_vars=['genome', 'sample', 'environment', 'phylum','genome_length'], value_vars=['regulatory_frequency', 'cds_frequency', 'non_coding_frequency','general_frequency'], var_name='region_type', value_name='frequency')
    print (melted_df_raw)
    ## remove elements with regulatory_frequency
    melted_df = melted_df_raw[(melted_df_raw['region_type'] != 'regulatory_frequency') & (melted_df_raw['region_type'] != 'general_frequency')]

    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(16, 16))

    # Filter environments with > 50 values
    env_counts = melted_df['environment'].value_counts()
    env_with_enough_data = env_counts[env_counts > ctg_num_cutoff].index
    melted_df_env = melted_df[melted_df['environment'].isin(env_with_enough_data)]
    
    # Filter phyla with > 50 values
    phylum_counts = melted_df['phylum'].value_counts()
    phylum_with_enough_data = phylum_counts[phylum_counts > ctg_num_cutoff].index
    
    
    sns.boxplot(data=melted_df_env, x='environment', y='frequency', hue='region_type', palette='Set2', ax=ax1)
    ax1.tick_params(axis='x', rotation=90)
    ax1.set_title('Frequency of Modified Bases in Different Genomic Regions Across Environments')
    ax1.set_xlabel(f'Environment (n>={ctg_num_cutoff})')
    ax1.set_ylabel('Frequency of Modified Bases')

    ## plot a boxplot with envs as x-axis and general frequency as y-axis, sorted by mean frequency
    general_freq_df = melted_df_raw[melted_df_raw['region_type'] == 'general_frequency']
    # Sort environments by mean frequency
    order = general_freq_df.groupby('environment')['frequency'].mean().sort_values().index
    sns.boxplot(data=general_freq_df, x='environment', y='frequency', order=order, palette='Set2', ax=ax2)
    ax2.tick_params(axis='x', rotation=90)
    ax2.set_title('General Frequency of Modified Bases Across Environments')
    ax2.set_xlabel(f'Environment (n>={ctg_num_cutoff})')
    ax2.set_ylabel('Mod Frequency')

    ## plot a boxplot with envs as x-axis and non-coding frequency as y-axis, sorted by mean frequency
    non_coding_freq_df = melted_df_raw[melted_df_raw['region_type'] == 'non_coding_frequency']
    # Sort environments by mean frequency
    order = non_coding_freq_df.groupby('environment')['frequency'].mean().sort_values().index
    sns.boxplot(data=non_coding_freq_df, x='environment', y='frequency', order=order, palette='Set2', ax=ax3)
    ax3.tick_params(axis='x', rotation=90)
    ax3.set_title('Non-coding Frequency of Modified Bases Across Environments')
    ax3.set_xlabel(f'Environment (n>={ctg_num_cutoff})')
    ax3.set_ylabel('Mod Frequency')

    melted_df_phylum = general_freq_df[general_freq_df['phylum'].isin(phylum_with_enough_data)]
    order = melted_df_phylum.groupby('phylum')['frequency'].mean().sort_values().index
    sns.boxplot(data=melted_df_phylum, x='phylum', y='frequency', order=order, palette='Set2', ax=ax4)
    ax4.tick_params(axis='x', rotation=90)
    ax4.set_title('General Frequency of Modified Bases Across Phyla')
    ax4.set_xlabel(f'Phylum (n>={ctg_num_cutoff})')
    ax4.set_ylabel('Frequency of Modified Bases')

    ## also plot two boxplots to regulatory_frequency across environment and phylum
    regulatory_freq_df = melted_df_raw[melted_df_raw['region_type'] == 'regulatory_frequency']
    # Sort environments by mean frequency
    order = regulatory_freq_df.groupby('environment')['frequency'].mean().sort_values().index
    sns.boxplot(data=regulatory_freq_df, x='environment', y='frequency', order=order, palette='Set2', ax=ax5)
    ax5.tick_params(axis='x', rotation=90)
    ax5.set_title('Regulatory Frequency of Modified Bases Across Environments')
    ax5.set_xlabel(f'Environment (n>={ctg_num_cutoff})')
    ax5.set_ylabel('Mod Frequency')

    melted_df_phylum_reg = regulatory_freq_df[regulatory_freq_df['phylum'].isin(phylum_with_enough_data)]
    order = melted_df_phylum_reg.groupby('phylum')['frequency'].mean().sort_values().index
    sns.boxplot(data=melted_df_phylum_reg, x='phylum', y='frequency', order=order, palette='Set2', ax=ax6)
    ax6.tick_params(axis='x', rotation=90)
    ax6.set_title('Regulatory Frequency of Modified Bases Across Phyla')
    ax6.set_xlabel(f'Phylum (n>={ctg_num_cutoff})')
    ax6.set_ylabel('Frequency of Modified Bases')

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/frequency_coding.png", dpi=300, bbox_inches='tight')
    melted_df.to_csv(f"{fig_dir}/frequency_coding.csv", index=False)
    ## also plot box plot coding vs non-coding and add the pvalue of t-test to  plot

    coding_df = melted_df[melted_df['region_type'] == 'cds_frequency']
    non_coding_df = melted_df[melted_df['region_type'] == 'non_coding_frequency']

    # Filter out samples with insufficient data for robust statistics
    print(f"Before filtering: Coding samples: {len(coding_df)}, Non-coding samples: {len(non_coding_df)}")
    
    # Remove NaN values and ensure minimum sample sizes
    coding_df_filtered = coding_df.dropna(subset=['frequency'])
    non_coding_df_filtered = non_coding_df.dropna(subset=['frequency'])
    
    print(f"After NaN filtering: Coding samples: {len(coding_df_filtered)}, Non-coding samples: {len(non_coding_df_filtered)}")
    
    # Only perform t-test if we have sufficient samples (>= 30 for each group)
    if len(coding_df_filtered) >= 30 and len(non_coding_df_filtered) >= 30:
        # Perform t-test
        t_stat, p_value = ttest_ind(coding_df_filtered['frequency'], non_coding_df_filtered['frequency'])
        print(f"T-test between coding and non-coding: t-statistic = {t_stat:.4f}, p-value = {p_value:.2e}")
        
        # Create comparison plot for coding vs non-coding
        plt.figure(figsize=(4, 4))
        sns.boxplot(data=melted_df, x='region_type', y='frequency', palette='Set2')
        plt.title(f'p={p_value:.2e}\n(n_coding={len(coding_df_filtered)}, n_non-coding={len(non_coding_df_filtered)})')
        plt.ylabel('Frequency of Modified Bases')
        plt.xticks([0, 1], ['CDS', 'Non-coding'])
        
        
        # Add significance bar
        y_max = max(coding_df_filtered['frequency'].max(), non_coding_df_filtered['frequency'].max())
        plt.plot([0, 1], [y_max * 1.05, y_max * 1.05], 'k-', linewidth=1)
        
        plt.tight_layout()
        plt.savefig(f"{fig_dir}/coding_vs_noncoding_comparison.pdf", dpi=300, bbox_inches='tight')
        plt.show()
    else:
        print("Warning: Insufficient sample sizes for t-test. Need at least 30 samples per group.")
        print(f"Current sample sizes - Coding: {len(coding_df_filtered)}, Non-coding: {len(non_coding_df_filtered)}")

    ## print top 10 genomes with highest coding frequency and also top 10 genomes with highest non-coding frequency
    top10_coding = coding_df.sort_values(by='frequency', ascending=False).head(20)
    top10_non_coding = non_coding_df.sort_values(by='frequency', ascending=False).head(20)
    print ("Top 10 genomes with highest coding frequency:")
    print (top10_coding[['genome', 'frequency', 'environment', 'phylum']])
    print ("Top 10 genomes with highest non-coding frequency:")
    print (top10_non_coding[['genome', 'frequency', 'environment', 'phylum']])
    ## also report these with highest general frequency, whole_df has general_frequency column
    top10_all = whole_df.sort_values(by='general_frequency', ascending=False).head(20)
    print ("Top 10 genomes with highest general frequency:")
    print (top10_all[['genome', 'general_frequency', 'environment', 'phylum']])

def count_good_ctgs(df_all_data):
    ## count how many unique species, phyla in the dataframe
    print (len(df_all_data), "contigs in total.")
    df_all_data['species'] = df_all_data['lineage'].apply(lambda x: x.split(";")[-1])
    unique_species_set = set()
    unique_phyla_set = set()
    no_species_count = 0
    for index, row in df_all_data.iterrows():
        species = row['species']
        phylum = row['phylum']
        ## remove 's__' from species name
        if species.startswith("s__"):
            species = species[3:]
        if species != "":
            unique_species_set.add(species)
        else:
            no_species_count += 1
        unique_phyla_set.add(phylum)
    print (f"Total unique species: {len(unique_species_set)}")
    print (f"Total unique phyla: {len(unique_phyla_set)}")
    print (f"Total contigs with no species annotation: {no_species_count}")
    ## count how many contigs have motif_num > 0
    motif_more_than_0 = df_all_data[df_all_data['motif_num'] > 0]
    print (f"Total contigs with motif number > 0: {len(motif_more_than_0)}")
    ## proportion
    proportion = len(motif_more_than_0) / len(df_all_data) * 100
    print (f"Proportion of contigs with motif number > 0: {proportion:.2f}%")
    ## print the contig name with 0 motifs
    motif_0_ctgs = df_all_data[df_all_data['motif_num'] == 0]['contig'].tolist()
    print (f"Contigs with 0 motifs: {motif_0_ctgs}")
    ## count the proportion of contigs with motif_num > 0 in each environment
    env_groups = df_all_data.groupby('environment')
    for env, group in env_groups:
        total_ctgs = len(group)
        motif_more_than_0_env = group[group['motif_num'] > 0]
        proportion_env = len(motif_more_than_0_env) / total_ctgs * 100
        print (f"Environment: {env}, Total contigs: {total_ctgs}, Contigs with motif number > 0: {len(motif_more_than_0_env)}, Proportion: {proportion_env:.2f}%")

    ## also count the proportion of contigs with motif_num > 0 in each phylum
    phylum_groups = df_all_data.groupby('phylum')
    for phylum, group in phylum_groups:
        total_ctgs = len(group)
        if total_ctgs < 40:
            continue
        motif_more_than_0_phylum = group[group['motif_num'] > 0]
        proportion_phylum = len(motif_more_than_0_phylum) / total_ctgs * 100
        print (f"Phylum: {phylum}, Total contigs: {total_ctgs}, Contigs with motif number > 0: {len(motif_more_than_0_phylum)}, Proportion: {proportion_phylum:.2f}%")

    ## check the average motif number in each environment
    for env, group in env_groups:
        avg_motif_num = group['motif_num'].mean()
        print (f"Environment: {env}, Average motif number: {avg_motif_num:.2f}")
    ## check the average motif number in each phylum, only for phyla with more than 40 contigs
    for phylum, group in phylum_groups:
        total_ctgs = len(group)
        if total_ctgs < 40:
            continue
        avg_motif_num = group['motif_num'].mean()
        print (f"Phylum: {phylum}, Average motif number: {avg_motif_num:.2f}")
    
    ## print the contig and lineage of top 5 contigs with highest motif number
    top5_motif = df_all_data.sort_values(by='motif_num', ascending=False).head(5)
    print ("Top 5 contigs with highest motif number:")
    for index, row in top5_motif.iterrows():
        print (f"Contig: {row['contig']}, Motif number: {row['motif_num']}, Lineage: {row['lineage']}")

    ## count how many are bacteria and how many are archaea
    bacteria_count = len(df_all_data[df_all_data['domain'] == 'd__Bacteria'])
    archaea_count = len(df_all_data[df_all_data['domain'] == 'd__Archaea'])
    print (f"Total Bacteria contigs: {bacteria_count}")
    print (f"Total Archaea contigs: {archaea_count}")


def color_tree(tree_file, tree_results, micro_type = "bacteria"):
    from Bio import Phylo
    if micro_type == "bacteria":
        PHYLUM_COLORS = BAC_PHYLUM_COLORS
    elif micro_type == "archaea":
        PHYLUM_COLORS = ARC_PHYLUM_COLORS
    df_all_data = pd.read_csv(f"{fig_dir}/motif_num_all_samples.csv")
    count_good_ctgs(df_all_data)
    motif_num_anno = gradient_header + "\n"
    env_anno = strip_heaer + "\n"
    phylum_anno = header + "\n"

    ## Use predefined PHYLUM_COLORS for top 8 most abundant phyla
    ## remove p__ from the phylum names in df_all_data
    df_all_data['phylum'] = df_all_data['phylum'].apply(lambda x: x[3:] if x.startswith("p__") else x)
    
    # Read the tree and count leaf nodes per phylum
    print(f"\nReading tree from: {tree_file}")
    tree = Phylo.read(tree_file, "newick")
    tree_leaf_names = [leaf.name for leaf in tree.get_terminals()]
    print(f"Total leaves in tree: {len(tree_leaf_names)}")
    
    # Count phyla based on tree leaf nodes
    tree_phylum_counts = {}
    for leaf_name in tree_leaf_names:
        # Find this leaf in the dataframe
        matching_rows = df_all_data[df_all_data['contig'] == leaf_name]
        if not matching_rows.empty:
            phylum = matching_rows.iloc[0]['phylum']
            tree_phylum_counts[phylum] = tree_phylum_counts.get(phylum, 0) + 1
    
    # Sort and get top 8
    sorted_phyla = sorted(tree_phylum_counts.items(), key=lambda x: x[1], reverse=True)
    top8_phyla = [p[0] for p in sorted_phyla[:14]]
    
    print(f"\nTop 8 phyla by tree node count:")
    for phylum, count in sorted_phyla[:20]:
        print(f"  {phylum}: {count} nodes")
    ## save all the phyla with node counts in to a CSV file
    phylum_counts_df = pd.DataFrame(sorted_phyla, columns=['Phylum', 'Node_Count'])
    phylum_counts_df.to_csv(f"{tree_results}/phylum_node_counts.csv", index=False)
    
    # Assign colors from PHYLUM_COLORS to top 8 phyla
    phylum_color = {}
    for phylum in top8_phyla:
        if phylum in PHYLUM_COLORS:
            phylum_color[phylum] = PHYLUM_COLORS[phylum]
        else:
            print (f"Phylum {phylum} not in predefined colors, assigning from palette.")
            # If phylum not in predefined colors, assign from palette
            color_idx = top8_phyla.index(phylum) % len(PHYLUM_COLORS)
            available_colors = [c for c in PHYLUM_COLORS.values() if c != PHYLUM_COLORS["Others"]]
            phylum_color[phylum] = list(available_colors)[color_idx]
    
    # All other phyla get the 'Others' color
    all_phyla = df_all_data['phylum'].unique()
    for phylum in all_phyla:
        if phylum not in phylum_color:
            phylum_color[phylum] = PHYLUM_COLORS["Others"]

    print(f"Top 8 phyla by node count and color:")
    for phylum in top8_phyla:
        print(f"  {phylum}: {tree_phylum_counts[phylum]} nodes" f" - Color: {phylum_color[phylum]}")

    ## prepare environment color map
    env_list = df_all_data['environment'].unique().tolist()
    env_color = {}
    color_palette_env = sns.color_palette("Set1", n_colors=len(env_list))
    for i, env in enumerate(env_list):
        rgb = color_palette_env[i]
        hex_color = '#%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
        env_color[env] = hex_color
    print (env_color)
    
    # Use top8_phyla as the known phyla list
    known_phyla = top8_phyla
    
    for index, row in df_all_data.iterrows():
        ctg_name = row['contig']
        if row['phylum'] not in known_phyla:
            phylum_anno += f"{ctg_name} range {phylum_color[row['phylum']]} Others\n"
        else:
            phylum_anno += f"{ctg_name} range {phylum_color[row['phylum']]} {row['phylum']}\n"

        if row['motif_num'] > 1:
            motif_num_anno += f"{ctg_name} 2\n"
        else:
            motif_num_anno += f"{ctg_name} {row['motif_num']}\n"
        env_anno += f"{ctg_name} {env_color[row['environment']]} {row['environment']}\n"
        
    with open(f"{tree_results}/motif_num_annotations.txt", "w") as f:
        f.write(motif_num_anno)
    with open(f"{tree_results}/env_annotations.txt", "w") as f:
        f.write(env_anno)

    with open(f"{tree_results}/phylum_annotations.txt", "w") as f:
        f.write(phylum_anno)

def plot_motif_len(fig_dir):
    df_all_data = pd.read_csv(f"{fig_dir}/motif_num_all_samples.csv")
    
    # Create motif number bins
    max_motif = df_all_data['motif_num'].max()
    df_all_data['motif_bin'] = pd.cut(df_all_data['motif_num'], 
                                      bins=range(0, max_motif + 2), 
                                      right=False, 
                                      labels=[f"{i}" for i in range(max_motif + 1)])
    
    # Filter out bins with very few samples (less than 5)
    bin_counts = df_all_data['motif_bin'].value_counts()
    valid_bins = bin_counts[bin_counts >= 5].index
    df_filtered = df_all_data[df_all_data['motif_bin'].isin(valid_bins)]
    
    # Create subplot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Scatter plot with reduced overlap
    sns.scatterplot(data=df_all_data, x='ctg_len', y='motif_num', hue='environment', 
                   palette='Set2', alpha=0.5, s=30, ax=ax1)
    ax1.set_xlabel('Contig Length')
    ax1.set_ylabel('Motif Number')
    ax1.set_title('Motif Number vs Contig Length')
    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Box plot of length distribution by motif number bins
    sns.boxplot(data=df_filtered, x='motif_bin', y='ctg_len', ax=ax2)
    ax2.set_xlabel('Motif Number Bin')
    ax2.set_ylabel('Contig Length')
    ax2.set_title('Contig Length Distribution by Motif Number Bins')
    ax2.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/motif_num_vs_contig_length.png", dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    meta_file = "/home/shuaiw/mGlu/assembly_pipe/prefix_table.tab"
    fig_dir = "../../tmp/figures/multi_env_linkage/"
    genome_list_file =  "/home/shuaiw/borg/paper/specificity/genome.list"
    all_dir = "/home/shuaiw/borg/paper/run2/"
    meta_dir = "/home/shuaiw/borg/paper/gene_anno/meta/"   ## get by gene_anno.smk from assembly_pipe/
    tree_results = "/home/shuaiw/borg/paper/specificity/tree/"
    bacteria_tree = "/home/shuaiw/borg/paper/isolation/GTDB_tree/meta_gtdb_tree/gtdbtk.unrooted.tree"
    archaea_tree = "/home/shuaiw/borg/paper/isolation/GTDB_tree/meta_gtdb_tree/archaea.unrooted.tree"
    archea_tree_results = "/home/shuaiw/borg/paper/specificity/archea_tree/"
    # color_tree(bacteria_tree, tree_results, "bacteria")
    # color_tree(archaea_tree, archea_tree_results, "archaea")
    # plot_motif_len(fig_dir)
    sample_env_dict = read_metadata(meta_file)
    main(all_dir, fig_dir, sample_env_dict)

    # plot_coding(meta_dir, fig_dir)
    # rerun(fig_dir)
    # get_stastics()
    # jaccard()
    # jaccard_batch()















    # fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
    # ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    # best_ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.circular.contigs.fa"
    # best_ctg_dir = "/home/shuaiw/methylation/data/borg/contigs/circular/"
    # profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"
    # anno_file = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_RM.rm.genes.tsv"
    # bin_file = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.tab"
    # gtdb_file = "/home/shuaiw/borg/contigs/GTDB/gtdbtk.all.summary.tsv"

    # work_dir = "/home/shuaiw/borg/bench/soil/run1/"
    # fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
    # prefix = "soil"

    # work_dir = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2"
    # fai = "/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa.fai"
    # prefix = "cow"

    # work_dir = "/home/shuaiw/borg/paper/infant/NANO_2_INF1340011_4PB"
    # fai = "/home/shuaiw/borg/allison/NANO_2_INF1340011_4PB_HR_HIFIASM_META_scaffold_min1000.fa.fai"
    # prefix = "infant"
    # prefix = "SRR23446540"
    # work_dir = f"/home/shuaiw/borg/paper/run/{prefix}/{prefix}_methylation"
    # fai = f"/home/shuaiw/borg/paper/run/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
    # all_host_file = f"/home/shuaiw/borg/paper/run/{prefix}/all_host_ctgs.tsv"


    # jaccard()
    # plot_MT_motif()
    # out_best_ctgs2(ref, best_ref, best_ctgs, best_ctg_dir)  ## split ctgs


    # out_best_ctgs(ref, best_ref, best_ctgs)
    # plot_heatmap(profile_file)
    # plot_strain_heatmap(profile_file)