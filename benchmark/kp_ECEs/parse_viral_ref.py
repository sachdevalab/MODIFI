#!/usr/bin/env python3
"""Parse RefSeq viral gbff -> host-annotated reference set for phage-host validation.
Writes viral_ref_hosted.fna (sequences with a documented host) + viral_ref_host.tsv
(accession, organism, host, host_genus)."""
import gzip, re, sys
from Bio import SeqIO

GBFF="/home/shuaiw/borg/revision/linked_eces/viral_ref/viral.1.genomic.gbff.gz"
FNA="/home/shuaiw/borg/revision/linked_eces/viral_ref/viral_ref_hosted.fna"
TSV="/home/shuaiw/borg/revision/linked_eces/viral_ref/viral_ref_host.tsv"

def host_genus(h):
    h=re.sub(r"^(uncultured|Candidatus)\s+","",str(h)).strip()
    return h.split(" ")[0] if h else ""

n=0; kept=0
with gzip.open(GBFF,"rt") as fh, open(FNA,"w") as fo, open(TSV,"w") as ft:
    ft.write("accession\torganism\thost\thost_genus\n")
    for rec in SeqIO.parse(fh,"genbank"):
        n+=1
        host=""
        for feat in rec.features:
            if feat.type=="source":
                for key in ("host","lab_host"):
                    if key in feat.qualifiers:
                        host=feat.qualifiers[key][0]; break
            if host: break
        if not host:      # only keep references with a documented host
            continue
        org=rec.annotations.get("organism","")
        acc=rec.id
        try:
            seq=str(rec.seq)
        except Exception:
            continue
        fo.write(f">{acc}\n{seq}\n")
        ft.write(f"{acc}\t{org}\t{host}\t{host_genus(host)}\n")
        kept+=1
        if n%5000==0: print(f"  parsed {n}, kept {kept}",flush=True)
print(f"DONE parsed {n} viral records, kept {kept} with documented host")
