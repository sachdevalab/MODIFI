#!/usr/bin/env python3
"""ECE criterion: exclude any ECE whose contig is placed to a PHYLUM by GTDB-Tk (enough bac120/ar53
markers to classify beyond domain) -- such contigs are chromosomal fragments, not extrachromosomal
elements. Domain-only 'Unclassified Bacteria/Archaea' is NOT excluded. Applied to the canonical final
set expanded/filterpass_FINAL.csv (backed up first). GTDB source = the metagenome taxa pickle."""
import pickle, shutil, re, os
import pandas as pd

A = "/home/shuaiw/borg/revision/ece_anno"
FIN = f"{A}/expanded/filterpass_FINAL.csv"
BACKUP = f"{A}/expanded/filterpass_FINAL.prephylumgate.csv"
TAXA = "/home/shuaiw/borg/paper/gene_anno/meta_ctg_taxa_dict.pkl"

def phylum(contig, d):
    t = d.get(contig)
    if not t or t in ("Unknown", "Unclassified") or re.search("Unclassified", str(t)):
        return ""
    for tok in str(t).split(";"):
        tok = tok.strip()
        if tok.startswith("p__") and len(tok) > 3:
            return tok
    return ""

def main():
    d = pickle.load(open(TAXA, "rb"))
    fp = pd.read_csv(FIN)
    fp["gtdb_phylum"] = fp["MGE"].map(lambda c: phylum(c, d))
    dropped = fp[fp.gtdb_phylum != ""]
    kept = fp[fp.gtdb_phylum == ""].copy()

    if not os.path.exists(BACKUP):
        shutil.copy2(FIN, BACKUP)
    kept.drop(columns=["gtdb_phylum"]).to_csv(FIN, index=False)

    print(f"backup -> {BACKUP}")
    print(f"filterpass_FINAL: {len(fp)} -> {len(kept)}  (dropped {len(dropped)} phylum-assigned ECEs)")
    print("by type kept:", kept.MGE_type.value_counts().to_dict())
    print("\ndropped contigs:")
    for _, r in dropped.iterrows():
        print(f"  {r.MGE:20s} {r.MGE_type:7s} len={r.mge_len:>9} {r.gtdb_phylum}")

if __name__ == "__main__":
    main()
