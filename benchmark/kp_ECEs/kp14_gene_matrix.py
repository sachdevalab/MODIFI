#!/usr/bin/env python3
"""Build gene-profile matrices for the 14 Kp-cluster representatives from eggNOG-mapper output:
  (a) representative x 26 COG functional categories (gene counts)
  (b) representative x curated ECE functional families (gene counts)
Plus per-rep has_toxin (eggNOG toxin-antitoxin terms) and has_virulence (abricate VFDB).
Writes long+wide source-data CSVs consumed by plot_kp14_gene_profile.R.
"""
import os, re, csv
import pandas as pd

GP = "/home/shuaiw/borg/revision/kp_eces/gene_profile"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
EMAP = f"{GP}/kp14.emapper.annotations"
os.makedirs(FIG, exist_ok=True)

COG_NAMES = {  # 26 COG functional categories
 "J":"Translation","A":"RNA processing","K":"Transcription","L":"Replication/repair",
 "B":"Chromatin","D":"Cell cycle/division","Y":"Nuclear structure","V":"Defense",
 "T":"Signal transduction","M":"Cell wall/membrane","N":"Cell motility","Z":"Cytoskeleton",
 "W":"Extracellular","U":"Secretion/trafficking","O":"PTM/chaperones","X":"Mobilome (phage/TE)",
 "C":"Energy","G":"Carbohydrate","E":"Amino acid","F":"Nucleotide","H":"Coenzyme","I":"Lipid",
 "P":"Inorganic ion","Q":"Secondary metab","R":"General function","S":"Unknown"}

# curated ECE functional families (priority order: first match wins)
FAMILIES = [
 ("Replication", r"replicat|\bRepA\b|\bRep\b protein|replication initiat|DnaA|primase|\bRepB\b|RepE|initiator"),
 ("Partitioning/stability", r"partition|ParA|ParB|\bSoj\b|StbA|ParM|segregation|plasmid stab|centromere"),
 ("Conjugation/mobilization", r"conjug|\bTra[A-Z]\b|\bTrb[A-Z]|\bVir[BD]|type IV secret|T4SS|relaxase|\bMob[A-Z]|mobiliz|conjugal|coupling protein|pilus|mating"),
 ("Transposon/integrase", r"transposase|integrase|resolvase|recombinase|\bIS[0-9]|\bTn[0-9]|insertion sequence|mobile element|invertase"),
 ("Toxin-antitoxin", r"toxin|antitoxin|CcdB|CcdA|RelE|RelB|ParE|MazF|MazE|VapC|VapB|HigB|HigA|HipA|\bDoc\b|\bPhd\b|Zeta|addiction"),
 ("Defense (RM/CRISPR/Abi)", r"restriction|modification methylase|DNA methyltransferase|methyltransferase|CRISPR|\bCas[0-9]|abortive infection|\bAbi[A-Z]|anti-restriction|BREX|retron"),
 ("Metal/biocide resistance", r"arsen|mercur|\bmer[A-Z]|copper|\bpco[A-Z]|silver|\bsil[A-Z]|tellur|cobalt-zinc-cadmium|\bczc|quaternary ammonium|\bqac"),
 ("AMR (antibiotic)", r"beta-lactamase|lactamase|aminoglycoside|chloramphenicol|tetracycline resist|macrolide|sulfonamid|dihydrofolate reductase|dihydropteroate|quinolone|antibiotic resist|efflux.*(drug|antibiot)"),
 ("Virulence", r"adhesin|invasin|hemolysin|haemolysin|enterotoxin|cytotoxin|siderophore|aerobactin|yersiniabactin|salmochelin|colibactin|capsul|fimbria|\bpili\b|hemagglutinin|type III secret|T3SS|effector"),
 ("Phage structural", r"capsid|terminase|portal|\btail\b|baseplate|major head|tail fiber|tail spike|phage.*(structural|coat)|prohead|scaffold"),
 ("Metabolism", r"dehydrogenase|kinase|transferase|reductase|synthase|synthetase|hydrolase|oxidase|permease|transporter|metaboli|biosynth|ABC transporter"),
]
FAM_ORDER = [f for f,_ in FAMILIES] + ["Other/hypothetical"]


def contig_of(orf):
    return orf.rsplit("_", 1)[0]


def read_emapper(path):
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t"); continue
            if line.startswith("#") or not line.strip():
                continue
            if header is None:
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return pd.DataFrame(rows)


def classify_family(text, cog):
    for fam, pat in FAMILIES:
        if re.search(pat, text, flags=re.IGNORECASE):
            return fam
    if "X" in str(cog):
        return "Transposon/integrase"     # mobilome COG catch
    return "Other/hypothetical"


def main():
    em = read_emapper(EMAP)
    qcol = "query" if "query" in em.columns else em.columns[0]
    em["rep"] = em[qcol].map(contig_of)
    reps_meta = pd.read_csv(f"{GP}/kp14_representatives.tsv", sep="\t")
    reps = list(reps_meta["representative"])
    em = em[em["rep"].isin(reps)]

    # text field for family classification
    def field(r, c): return str(r.get(c, "") or "")
    em["text"] = (em.get("Description","").astype(str) + " | " + em.get("Preferred_name","").astype(str)
                  + " | " + em.get("PFAMs","").astype(str) + " | " + em.get("KEGG_ko","").astype(str))
    em["family"] = [classify_family(t, c) for t, c in zip(em["text"], em.get("COG_category",""))]

    # (a) COG category matrix
    cog_rows = []
    for _, r in em.iterrows():
        cats = re.findall(r"[A-Z]", str(r.get("COG_category","")))
        for c in cats:
            if c in COG_NAMES:
                cog_rows.append((r["rep"], c))
    cog = pd.DataFrame(cog_rows, columns=["rep","cog"])
    cog_mat = cog.groupby(["rep","cog"]).size().unstack(fill_value=0).reindex(index=reps, fill_value=0)
    cog_mat = cog_mat.reindex(columns=[c for c in COG_NAMES if c in cog_mat.columns], fill_value=0)
    cog_mat.to_csv(f"{FIG}/kp14_cog_matrix.csv")

    # (b) curated family matrix
    fam = em.groupby(["rep","family"]).size().unstack(fill_value=0).reindex(index=reps, fill_value=0)
    fam = fam.reindex(columns=[f for f in FAM_ORDER if f in fam.columns], fill_value=0)
    fam.to_csv(f"{FIG}/kp14_family_matrix.csv")

    # toxin flag (eggNOG TA terms)
    ta = em[em["family"]=="Toxin-antitoxin"].groupby("rep").size()
    # virulence from abricate VFDB
    vir = {}
    ab = f"{GP}/kp14_abricate_vfdb.tsv"
    if os.path.exists(ab):
        adf = pd.read_csv(ab, sep="\t", comment=None)
        adf = adf[adf["#FILE"].notna()] if "#FILE" in adf.columns else adf
        seqcol = "SEQUENCE" if "SEQUENCE" in adf.columns else adf.columns[1]
        genecol = "GENE" if "GENE" in adf.columns else None
        for _, r in adf.iterrows():
            s = str(r[seqcol])
            if s in reps:
                vir.setdefault(s, []).append(str(r[genecol]) if genecol else "VF")

    ann = reps_meta.copy()
    ann["n_genes_annotated"] = ann["representative"].map(em.groupby("rep").size()).fillna(0).astype(int)
    ann["n_toxin_genes"] = ann["representative"].map(ta).fillna(0).astype(int)
    ann["has_toxin"] = ann["n_toxin_genes"] > 0
    ann["vfdb_genes"] = ann["representative"].map(lambda r: ";".join(sorted(set(vir.get(r,[])))))
    ann["has_virulence"] = ann["vfdb_genes"] != ""

    # --- mobilization type (mob_typer) ---
    mob_f = f"{GP}/kp14_mobtyper.tsv"
    if os.path.exists(mob_f):
        mob = pd.read_csv(mob_f, sep="\t")
        mob["rep"] = mob["sample_id"].astype(str).str.replace(".fna","",regex=False)
        mob = mob.set_index("rep")
        def mget(r,c): return mob.loc[r,c] if r in mob.index and c in mob.columns else ""
        ann["predicted_mobility"] = ann["representative"].map(lambda r: mget(r,"predicted_mobility"))
        ann["relaxase_MOB_type"] = ann["representative"].map(lambda r: mget(r,"relaxase_type(s)"))
        ann["mpf_type"] = ann["representative"].map(lambda r: mget(r,"mpf_type"))
        ann["replicon_type"] = ann["representative"].map(lambda r: mget(r,"rep_type(s)"))

    # --- AMR / stress genes (AMRFinderPlus, nucleotide) ---
    amr_f = f"{GP}/kp14_amrfinder.tsv"
    if os.path.exists(amr_f) and os.path.getsize(amr_f) > 0:
        amr = pd.read_csv(amr_f, sep="\t")
        ctg = "Contig id" if "Contig id" in amr.columns else amr.columns[2]
        et = "Element type"; gs = "Gene symbol"
        amr_g = amr[amr[et]=="AMR"].groupby(ctg)[gs].apply(lambda s: ";".join(sorted(set(s))))
        str_g = amr[amr[et]=="STRESS"].groupby(ctg)[gs].apply(lambda s: ";".join(sorted(set(s))))
        vir_g = amr[amr[et]=="VIRULENCE"].groupby(ctg)[gs].apply(lambda s: ";".join(sorted(set(s)))) if (amr[et]=="VIRULENCE").any() else {}
        ann["amr_genes"] = ann["representative"].map(amr_g).fillna("")
        ann["stress_genes"] = ann["representative"].map(str_g).fillna("")
        ann["amrfinder_virulence_genes"] = ann["representative"].map(vir_g).fillna("") if len(vir_g) else ""
        ann["n_amr"] = ann["amr_genes"].apply(lambda s: 0 if s=="" else len(s.split(";")))
        ann["n_stress"] = ann["stress_genes"].apply(lambda s: 0 if s=="" else len(s.split(";")))
        # broaden virulence flag to include AMRFinder VIRULENCE + VFDB
        ann["has_virulence"] = ann["has_virulence"] | (ann.get("amrfinder_virulence_genes","")!="")
    ann.to_csv(f"{FIG}/kp14_annotation.tsv", sep="\t", index=False)
    ann.to_csv(f"{GP}/kp14_annotation.tsv", sep="\t", index=False)

    print(f"[matrix] eggNOG-annotated ORFs on reps: {len(em)}")
    print(f"[matrix] COG matrix {cog_mat.shape}, family matrix {fam.shape}")
    print(f"[matrix] has_toxin: {int(ann['has_toxin'].sum())}/14 ; has_virulence(VFDB): {int(ann['has_virulence'].sum())}/14")
    print("\n[matrix] curated family totals across reps:")
    print(fam.sum().sort_values(ascending=False).to_string())


if __name__ == "__main__":
    main()
