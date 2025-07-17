import pandas as pd
import re
from Bio import SeqIO
import sys
import os


def get_circular_ctgs(fai):
    """
    Get the best contig based on length from a fasta file.
    """
    circular_ctgs = {}
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            if ctg[-1] == "C":
                circular_ctgs[ctg] = length
    print (f"Total {len(circular_ctgs)} circular contigs found.")
    return circular_ctgs

def get_high_quality_ctgs(checkM_report):
    """
    Get the high quality contigs from a CheckM report.
    """
    high_quality_ctgs = []
    df = pd.read_csv(checkM_report, sep="\t", header=0)
    df = df[df['Completeness'] >= 50]
    df = df[df['Contamination'] <= 5]
    for ctg in df['Name']:
        high_quality_ctgs.append(ctg)

    print (f"Total {len(high_quality_ctgs)} high quality contigs found.")
    return high_quality_ctgs

def get_best_ctg(fai, checkM_report, high_quality_file, classification_dict):
    high_quality_ctgs = get_high_quality_ctgs(checkM_report)
    circular_ctgs = get_circular_ctgs(fai)
    # selected_ctgs = list(set(high_quality_ctgs) | set(circular_ctgs))
    selected_host_ctgs = set()
    novel_elements = {}

    for circular in circular_ctgs:
        if circular not in classification_dict:
            print(f"Warning: {circular} is not found in classification dictionary.")
        else:
            classification = classification_dict[circular]
            if classification == 'Bacteria' or classification == 'Archaea':
                selected_host_ctgs.add(circular)
            else:
                novel_elements[circular] = circular_ctgs[circular]
    selected_host_ctgs = selected_host_ctgs | set(high_quality_ctgs)
    print (f"Total {len(selected_host_ctgs)} host contigs selected.")
    print (f"Total {len(novel_elements)} novel elements found.")
    print (f"Novel elements: {novel_elements}")
    # return selected_ctgs
    f = open(high_quality_file, "w")
    for ctgs in selected_host_ctgs:
        f.write(f"{ctgs}\t{ctgs}\t{classification_dict[ctgs]}\n")
    f.close()
    return selected_host_ctgs, novel_elements

def classify(classification):
    if re.search(r'd__Bacteria', classification) or re.search(r'Unclassified Bacteria', classification):
        return 'Bacteria'
    elif re.search(r'd__Archaea', classification) or re.search(r'Unclassified Archaea', classification):
        return 'Archaea'
    else:
        return classification

def get_bacteria_archea(gtdb_ar53, gtdb_bac120):
    classification_dict = {}
    ## check if the files exist
    if os.path.exists(gtdb_ar53):
        df_ar53 = pd.read_csv(gtdb_ar53, sep="\t")
        for index, row in df_ar53.iterrows():
            classification_dict[row['user_genome']] = classify(row['classification'])
    else:
        print(f"Warning: {gtdb_ar53} does not exist. Skipping GTDB AR53 classification.")
    if os.path.exists(gtdb_bac120):
        df_bac120 = pd.read_csv(gtdb_bac120, sep="\t")
        for index, row in df_bac120.iterrows():
            classification_dict[row['user_genome']] = classify(row['classification'])
    else:
        print(f"Warning: {gtdb_bac120} does not exist. Skipping GTDB BAC120 classification.")
    # print(classification_dict)
    return classification_dict




class mge_obj:

    def __init__(self, name, type, length, method):
        self.name = name
        self.type = type
        self.length = length
        self.methods = [method]
        self.classification = None
    
    def add_method(self, method):
        if method not in self.methods:
            self.methods.append(method)


def read_virsorter2(virsorter2):
    """
    Read the virsorter2 output file and return a DataFrame.
    """
    virus_dict = {}
    df = pd.read_csv(virsorter2, sep="\t", header=0)
    for index, row in df.iterrows():
        seqname=row['seqname']
        type="virus"
        length=row['length']
        field = seqname.split("||")
        if re.search("_partial", field[1]):
            continue
        mge_name = field[0]
        mge = mge_obj(
            name=mge_name,
            type=type,
            length=length,
            method="virsorter2"
        )
        # Do something with the mge object
        virus_dict[mge_name] = mge
    return virus_dict

def read_vibrant(virbrant_fna):
    vibrant_dict = {}
    for record in SeqIO.parse(virbrant_fna, "fasta"):
        seq_name = record.id
        length = len(record.seq)
        field = seq_name.split("_")
        if field[-2] == "fragment":
            continue
        mge = mge_obj(
            name=seq_name,
            type="virus",
            length=length,
            method="vibrant"
        )
        vibrant_dict[seq_name] = mge
    return vibrant_dict

def read_genomad(genomad_plasmid, mge_type="plasmid"):
    print (f"Reading {genomad_plasmid}...")
    genomad_dict = {}
    genomad = pd.read_csv(genomad_plasmid, sep = "\t")
    for i, row in genomad.iterrows():
        
        if re.search('\|provirus', row['seq_name']):
            continue
        if row['seq_name'] == 'seq_name':
            continue
        mge = mge_obj(
            name=row['seq_name'],
            type=mge_type,
            length=row['length'],
            method="genomad"
        )
        genomad_dict[row['seq_name']] = mge
    return genomad_dict



def get_novel_elements(novel_elements):
    novel_dict = {}
    for novel in novel_elements:
        mge = mge_obj(
            name=novel,
            type="novel",
            length=novel_elements[novel],  # Length is not provided in the classification dictionary
            method="in-house"
        )
        novel_dict[novel] = mge
    return novel_dict

def get_mge_union(all_mge_file, classification_dict, vibrant_virus, virsorter2_virus, genomad_virus, genomad_plasmid, novel_dict):
    all_mge = list((vibrant_virus.keys() | virsorter2_virus.keys() | genomad_virus.keys() | genomad_plasmid.keys() | novel_dict.keys()))
    all_mge = set(all_mge)
    data = []
    for mge in all_mge:
        methods = []
        if mge in genomad_virus:
            methods.append("genomad")
            mge = genomad_virus[mge] 
        if mge in genomad_plasmid:
            methods.append("genomad")
            mge = genomad_plasmid[mge]
        if mge in vibrant_virus:
            methods.append("vibrant")
            mge = vibrant_virus[mge]
        if mge in virsorter2_virus:
            methods.append("virsorter2")
            mge = virsorter2_virus[mge]
        if mge in novel_dict:
            methods.append("in-house")
            mge = novel_dict[mge]
        mge.methods = methods
        classification = classification_dict.get(mge.name, "no_info")
        data.append({
            "seq_name": mge.name,
            "type": mge.type,
            "length": mge.length,
            "methods": ",".join(mge.methods),
            "GTDB": classification
        })
    df = pd.DataFrame(data)
    df = df.sort_values(by=['type', 'length'], ascending=[True, False])
    print (df)
    df.to_csv(all_mge_file, sep="\t", index=False)



if __name__ == "__main__":
    workdir = sys.argv[1]
    prefix = sys.argv[2]

    fai = f"{workdir}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
    checkM_report = f"{workdir}/checkM2/quality_report.tsv"
    high_quality_file = f"{workdir}/all_host_ctgs.tsv"
    gtdb_ar53 = f"{workdir}/GTDB/gtdbtk.ar53.summary.tsv"
    gtdb_bac120 = f"{workdir}/GTDB/gtdbtk.bac120.summary.tsv"
    classification_dict = get_bacteria_archea(gtdb_ar53, gtdb_bac120)
    selected_host_ctgs, novel_elements = get_best_ctg(fai, checkM_report, high_quality_file, classification_dict)

    genomad_plasmid = f"{workdir}/Genomad/{prefix}.hifiasm.p_ctg.rename_summary/{prefix}.hifiasm.p_ctg.rename_plasmid_summary.tsv"
    genomad_virus = f"{workdir}/Genomad/{prefix}.hifiasm.p_ctg.rename_summary/{prefix}.hifiasm.p_ctg.rename_virus_summary.tsv"
    virsorter2 = f"{workdir}/virsorter2/final-viral-score.tsv"
    virbrant = f"{workdir}/vibrant/VIBRANT_{prefix}.hifiasm.p_ctg.rename/VIBRANT_phages_{prefix}.hifiasm.p_ctg.rename/{prefix}.hifiasm.p_ctg.rename.phages_combined.fna"
    all_mge_file = f"{workdir}/all_mge.tsv"



    vibrant_virus = read_vibrant(virbrant)
    virsorter2_virus = read_virsorter2(virsorter2)
    genomad_virus = read_genomad(genomad_virus, mge_type="virus")
    genomad_plasmid = read_genomad(genomad_plasmid, mge_type="plasmid")
    novel_dict = get_novel_elements(novel_elements)

    get_mge_union(all_mge_file, classification_dict, vibrant_virus, virsorter2_virus, genomad_virus, genomad_plasmid, novel_dict)
