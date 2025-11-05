# -*- coding: utf-8 -*-
import subprocess
import os
from pathlib import Path
import sys
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

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


legend_header = """DATASET_COLORSTRIP
#In colored strip datasets, each ID is associated to a color box/strip and can have an optional label. Color can be specified in hexadecimal, RGB or RGBA notation. When using RGB or RGBA notation, you cannot use COMMA as the dataset separator

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

SEPARATOR SPACE
#SEPARATOR TAB
#SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL Phylum_Colors

#dataset color (can be changed later)
DATASET_COLOR #ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#If COLOR_BRANCHES is set to 1, branches of the tree will be colored according to the colors of the strips above the leaves.
#When all children of a node have the same color, it will be colored the same, otherwise it will be black.
COLOR_BRANCHES 0

#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#

#Each dataset can have a legend, which is defined using LEGEND_XXX fields
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark

LEGEND_TITLE Phylum
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 0
LEGEND_SHAPES"""

def get_taxa(isolation_taxa, name):
    # sra_id = contig.split(".")[0]
    if name in isolation_taxa:
        return isolation_taxa[name]
    else:
        return "unknown"

def read_gtdb(gatk):
    """
    Read the GTDB summary file and return a dictionary of contig to bin mapping.
    """
    gtdb_df = pd.read_csv(gatk, sep='\t')
    isolation_taxa = {}
    for index, row in gtdb_df.iterrows():
        anno = row['classification']
        sra_id = row['user_genome']
        isolation_taxa[sra_id] = anno
        return isolation_taxa
        
def single_run():
    resultdir = f"/groups/banfield/projects/multienv/methylation_temp/batch2_results/"
    run_taxa_dict = {}
    for folder in os.listdir(resultdir):
        prefix = folder
        outdir = f"{resultdir}/{prefix}"
        reference_fasta = f"{outdir}/{prefix}.hifiasm.p_ctg.rename.fa"
        gtdb = f"{outdir}/GTDB/gtdbtk.bac120.summary.tsv"
        ## skip if gtdb file not found
        if not os.path.exists(gtdb):
            print (f"GTDB file not found for {prefix}, skipping...")
            continue
        isolation_taxa = read_gtdb(gtdb)
        bin_name = f"{prefix}.hifiasm.p_ctg.rename"
        taxa = get_taxa(isolation_taxa, bin_name)
        print (f"Processing {prefix} with taxa {taxa}", bin_name)
        run_taxa_dict[prefix] = taxa

    return run_taxa_dict

def color_phylum():
    ## in iTol, color by phylum using general colors
    ## collect phyllums in run_taxa_dict, and assign colors
    taxa_list = set(list(run_taxa_dict.values()))
    
    # Extract unique phylums
    phylums = set()
    for taxa in run_taxa_dict.values():
        if taxa != "unknown":
            phylum = taxa.split(";")[1] if len(taxa.split(";")) > 1 else "unknown"
            phylum = phylum.replace("p__", "")
            phylums.add(phylum)
    
    # Use seaborn Set1 color palette
    phylum_list = sorted(list(phylums))  # Sort for consistent coloring
    n_colors = len(phylum_list)
    
    # Get Set1 palette colors (Set1 has 9 distinct colors)
    set1_colors = sns.color_palette("Set1", n_colors=max(n_colors, 9))
    
    # Convert RGB to hex format for iTOL
    phylum_color = {}
    for i, phylum in enumerate(phylum_list):
        if i < len(set1_colors):
            rgb = set1_colors[i]
            hex_color = "#{:02x}{:02x}{:02x}".format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            phylum_color[phylum] = hex_color
        else:
            # If more phylums than Set1 colors, cycle through
            rgb = set1_colors[i % len(set1_colors)]
            hex_color = "#{:02x}{:02x}{:02x}".format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            phylum_color[phylum] = hex_color
    
    # Add default color for unknown
    phylum_color["unknown"] = "#cccccc"
    
    print(f"Found {len(phylum_list)} unique phylums using Seaborn Set1 palette:")
    for phylum, color in phylum_color.items():
        print(f"  {phylum}: {color}")
    
    # Generate legend for iTOL

    
    # Add legend shapes and colors for each phylum
    legend_shapes = []
    legend_colors = []
    legend_labels = []
    
    for phylum in sorted(phylum_color.keys()):
        if phylum != "unknown":
            legend_shapes.append("1")  # square shape
            legend_colors.append(phylum_color[phylum])
            legend_labels.append(phylum)
    
    legend_content = legend_header + " " + " ".join(legend_shapes) + "\n"
    legend_content += "LEGEND_COLORS " + " ".join(legend_colors) + "\n"
    legend_content += "LEGEND_LABELS " + " ".join(legend_labels) + "\n\n"
    legend_content += "#=================================================================#\n"
    legend_content += "#     DATA SECTION                                                #\n"  
    legend_content += "#=================================================================#\n"
    legend_content += "DATA\n"
    legend_content += "#ID COLOR LABEL\n"
    
    ## generate annotations file for iTol
    Label_anno = Label_header
    
    # Generate color strip data for legend
    legend_data = ""
    
    with open("/groups/banfield/projects/multienv/methylation_temp/GTDB_tree/phylum_annotations.txt", "w") as f:
        f.write(f"{header}\n")
        for bin_name, taxa in run_taxa_dict.items():
            phylum = taxa.split(";")[1] if len(taxa.split(";")) > 1 else "unknown"
            phylum = phylum.replace("p__", "")
            color = phylum_color.get(phylum, "#cccccc")
            
            # Extract species information (typically at index 6 in GTDB taxonomy)
            species = taxa.split(";")[6] if len(taxa.split(";")) > 6 else "unknown"
            species = species.replace("s__", "")

            f.write(f"{bin_name} label {color}\n")
            node_rename = f"{species}_{bin_name}"
            Label_anno += f"{bin_name},{node_rename}\n"
            
            # Add to legend data
            legend_data += f"{bin_name} {color} {phylum}\n"
    
    # Save legend file with color strip format
    with open("/groups/banfield/projects/multienv/methylation_temp/GTDB_tree/phylum_legend.txt", "w") as f:
        f.write(legend_content + legend_data)

    print(f"\nGenerated iTOL annotation file with {len(run_taxa_dict)} entries using Set1 colors")
    print(f"Generated legend file: phylum_legend.txt")
    print(f"Legend includes {len([p for p in phylum_color.keys() if p != 'unknown'])} phylums")
    print(run_taxa_dict)
    
    ## save label annotation file
    with open("/groups/banfield/projects/multienv/methylation_temp/GTDB_tree/label_annotations.csv", "w") as f:
        f.write(Label_anno)



run_taxa_dict = single_run()
color_phylum()