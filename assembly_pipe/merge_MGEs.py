import pandas as pd
import re
from Bio import SeqIO

class mge_obj:

    def __init__(self, name, type, length):
        self.name = name
        self.type = type
        self.length = length


def read_virsorter2(virsorter2):
    """
    Read the virsorter2 output file and return a DataFrame.
    """
    virus_list= []
    df = pd.read_csv(virsorter2, sep="\t", header=0)
    for index, row in df.iterrows():
        seqname=row['seqname']
        type="virus"
        length=row['length']
        field = seqname.split("||")
        if field[1] == "full":
            continue
        mge_name = field[0]
        mge = mge_obj(
            name=mge_name,
            type=type,
            length=length
        )
        # Do something with the mge object
        virus_list.append(mge)
    return virus_list

def read_vibrant(virbrant_fna):
    vibrant_virus = []
    for record in SeqIO.parse(virbrant_fna, "fasta"):
        seq_name = record.id
        length = len(record.seq)
        field = seq_name.split("_")
        if field[-2] == "fragment":
            continue
        mge = mge_obj(
            name=seq_name,
            type="virus",
            length=length
        )
        vibrant_virus.append(mge)
    return vibrant_virus

def read_genomad(genomad_plasmid, mge_type="plasmid"):
    print (f"Reading {genomad_plasmid}...")
    genomad_mge = []
    genomad = pd.read_csv(genomad_plasmid, sep = "\t")
    for i, row in genomad.iterrows():
        
        if re.search('\|provirus', row['seq_name']):
            continue
        if row['seq_name'] == 'seq_name':
            continue
        mge = mge_obj(
            name=row['seq_name'],
            type=mge_type,
            length=row['length']
        )
        genomad_mge.append(mge)
    return genomad_mge

def get_mge_union(all_mge, all_mge_file):
    data = []
    for mge in all_mge:
        data.append({
            "name": mge.name,
            "type": mge.type,
            "length": mge.length
        })
    df = pd.DataFrame(data)
    df = df.drop_duplicates(subset=['name', 'type'])
    df = df.sort_values(by=['type', 'length'], ascending=[True, False])
    print (df)
    df.to_csv(all_mge_file, sep="\t", index=False)


if __name__ == "__main__":
    workdir = "/home/shuaiw/borg/assembly/96plex/96plex_p5_2"
    prefix = "96plex_p5"

    genomad_plasmid = f"{workdir}/Genomad/{prefix}.final_summary/{prefix}.final_plasmid_summary.tsv"
    genomad_virus = f"{workdir}/Genomad/{prefix}.final_summary/{prefix}.final_virus_summary.tsv"
    virsorter2 = f"{workdir}/virsorter2/final-viral-score.tsv"
    virbrant = f"{workdir}/vibrant/VIBRANT_{prefix}.final/VIBRANT_phages_{prefix}.final/{prefix}.final.phages_combined.fna"
    all_mge_file = f"{workdir}/all_mge.tsv"

    vibrant_virus = read_vibrant(virbrant)
    virsorter2_virus = read_virsorter2(virsorter2)
    genomad_virus = read_genomad(genomad_virus, mge_type="virus")
    genomad_plasmid = read_genomad(genomad_plasmid, mge_type="plasmid")
    all_mge = vibrant_virus + virsorter2_virus + genomad_virus + genomad_plasmid
    get_mge_union(all_mge, all_mge_file)
