"""
Given the checkM's result, plot dot plot where x is contamination and y is completeness
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple
from Bio import SeqIO




def read_report(report):
    df = pd.read_csv(report, sep="\t")
    ## plot contamination and completeness
    df = df[["Contamination", "Completeness"]]
    df.columns = ["contamination", "completeness"]
    df = df.dropna()
    ## plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df, x="contamination", y="completeness")
    plt.xlabel("Contamination")
    plt.ylabel("Completeness")
    plt.title("Contamination vs Completeness")
    ## save figure
    plt.savefig(report.replace(".tsv", ".png"))
    plt.close()

    print (len(df))
    high_quality_bins = len(df[(df["completeness"] >= 90) & (df["contamination"] <= 5)])
    print ("Number of high quality bins: ", high_quality_bins,high_quality_bins/len(df)) 
    ## count the number of bins with medium quality bins
    medium_quality_bins = len(df[(df["completeness"] >= 50) & (df["contamination"] <= 10)])
    print ("Number of medium quality bins: ", medium_quality_bins,medium_quality_bins/len(df)) 

    ## count the number of bins with low quality bins
    low_quality_bins = len(df[(df["completeness"] < 50)])
    print ("Number of low quality bins: ", low_quality_bins,low_quality_bins/len(df))

    ## count the number of bins with completeness < 30
    df["completeness"] = df["completeness"].astype(float)
    df["contamination"] = df["contamination"].astype(float)
    df1 = df[df["completeness"] < 30]
    print("Number of bins with completeness < 30: ", len(df1))
    ### count the number of bins with contamination > 20
    df2 = df[df["contamination"] > 20]
    print("Number of bins with contamination > 20: ", len(df2))

    return df

def categorize_quality(df):
    """Categorize bins based on completeness and contamination"""
    # Initialize with default category
    df["quality_category"] = "Low Quality"
    
    # Apply conditions in order (most restrictive first)
    # High quality: completeness >= 90 AND contamination <= 5
    high_quality_mask = (df["completeness"] >= 90) & (df["contamination"] <= 5)
    df.loc[high_quality_mask, "quality_category"] = "High Quality"
    
    # Medium quality: completeness >= 50 AND contamination <= 10 (but not high quality)
    medium_quality_mask = (df["completeness"] >= 50) & (df["contamination"] <= 10) & (~high_quality_mask)
    df.loc[medium_quality_mask, "quality_category"] = "Medium Quality"
    
    # Low quality: everything else (completeness < 50 OR other conditions not met)
    # This is already set as default above
    
    return df

def plot_quality_categories(df1, df2):
    """Compare two methods using bar plot for quality categories"""
    

    
    # Add quality categories to both dataframes
    df1_cat = categorize_quality(df1.copy())
    df2_cat = categorize_quality(df2.copy())
    
    # Add method labels
    df1_cat["method"] = "Methylation-based"
    df2_cat["method"] = "MetaBAT2"
    
    # Combine dataframes
    combined_df = pd.concat([df1_cat, df2_cat], ignore_index=True)
    
    # Count bins in each category for each method
    category_counts = combined_df.groupby(["method", "quality_category"]).size().reset_index(name="count")
    print("\nBin counts by quality category:")
    print(category_counts)
    
    # Create pivot table for easier plotting
    pivot_counts = category_counts.pivot(index="quality_category", columns="method", values="count").fillna(0)
    
    # Create bar plot
    plt.figure(figsize=(12, 8))
    ax = pivot_counts.plot(kind="bar", width=0.8, color=["#2E86C1", "#E74C3C"])
    plt.title("Quality Categories Comparison Between Methods", fontsize=16, fontweight="bold")
    plt.xlabel("Quality Category", fontsize=14)
    plt.ylabel("Number of Bins", fontsize=14)
    plt.legend(title="Method", fontsize=12, title_fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", alpha=0.3)
    
    # Add value labels on bars
    for container in ax.containers:
        ax.bar_label(container, fmt='%d', fontsize=10)
    
    plt.tight_layout()
    plt.savefig("../../tmp/results/quality_categories_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    # Create percentage comparison
    plt.figure(figsize=(12, 8))
    
    # Calculate percentages
    total_counts = combined_df.groupby("method").size()
    percentage_data = []
    
    for method in ["Methylation-based", "MetaBAT2"]:
        method_data = combined_df[combined_df["method"] == method]
        total = len(method_data)
        
        for category in ["High Quality", "Medium Quality", "Low Quality"]:
            count = len(method_data[method_data["quality_category"] == category])
            percentage = (count / total) * 100 if total > 0 else 0
            percentage_data.append({"method": method, "quality_category": category, "percentage": percentage})
    
    percentage_df = pd.DataFrame(percentage_data)
    pivot_percentage = percentage_df.pivot(index="quality_category", columns="method", values="percentage").fillna(0)
    
    ax2 = pivot_percentage.plot(kind="bar", width=0.8, color=["#2E86C1", "#E74C3C"])
    plt.title("Quality Categories Comparison (Percentage)", fontsize=16, fontweight="bold")
    plt.xlabel("Quality Category", fontsize=14)
    plt.ylabel("Percentage of Bins (%)", fontsize=14)
    plt.legend(title="Method", fontsize=12, title_fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", alpha=0.3)
    
    # Add percentage labels on bars
    for container in ax2.containers:
        ax2.bar_label(container, fmt='%.1f%%', fontsize=10)
    
    plt.tight_layout()
    plt.savefig("../../tmp/results/quality_categories_percentage_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    print("\nPercentage breakdown:")
    print(pivot_percentage)

# report = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break/checkM2/quality_report.tsv"
# report = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break/checkM_metabat/quality_report.tsv"
our_report = "/home/shuaiw/borg/bench/soil/run2/checkm_methy/quality_report.tsv"
metabat_report = "/home/shuaiw/borg/bench/soil/run2/checkm_metabat/quality_report.tsv"
our_df = read_report(our_report)
metabat_df = read_report(metabat_report)
plot_quality_categories(our_df, metabat_df)