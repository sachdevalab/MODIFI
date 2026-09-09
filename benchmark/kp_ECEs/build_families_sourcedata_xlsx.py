#!/usr/bin/env python3
"""Excel source data for fig_kp14_gene_profile_families.pdf.
Sheet1 'genes'          : one row per predicted gene on each representative -> contig, gene_ID,
                          gene_category (curated functional family) + the eggNOG annotation used to
                          classify it (the 'why' at gene level). This underlies the heatmap cell counts.
Sheet2 'contig_annotation': one row per (contig, annotation category) -> contig, contig_type, category,
                          value (yes/no or the call), and why (the evidence). This underlies the side bars.
"""
import re, csv
import pandas as pd

GP = "/home/shuaiw/borg/revision/kp_eces/gene_profile"
FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
EMAP = f"{GP}/kp14.emapper.annotations"
OUT = f"{FIG}/fig_kp14_gene_profile_families_sourcedata.xlsx"

# ---- curated families (identical to kp14_gene_matrix.py) ----
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

def classify(text, cog):
    for fam, pat in FAMILIES:
        if re.search(pat, text, flags=re.IGNORECASE):
            return fam
    if "X" in str(cog):
        return "Transposon/integrase"
    return "Other/hypothetical"

reps = pd.read_csv(f"{GP}/kp14_representatives.tsv", sep="\t")
order = list(reps.sort_values("rep_length", ascending=False)["representative"])
rep_set = set(reps["representative"])

# ---- read eggNOG annotations ----
def read_emapper(path):
    rows, header = [], None
    for line in open(path):
        if line.startswith("##"): continue
        if line.startswith("#query"): header = line.lstrip("#").rstrip("\n").split("\t"); continue
        if line.startswith("#") or not line.strip() or header is None: continue
        rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows

# ---- Sheet1: genes ----
gene_rows = []
for r in read_emapper(EMAP):
    q = r.get("query", "")
    contig = q.rsplit("_", 1)[0]
    if contig not in rep_set:
        continue
    desc = r.get("Description", "-"); pref = r.get("Preferred_name", "-"); pf = r.get("PFAMs", "-")
    text = f"{desc} | {pref} | {pf} | {r.get('KEGG_ko','')}"
    fam = classify(text, r.get("COG_category", ""))
    ann = "; ".join(x for x in [pref if pref not in ("-", "") else "",
                                desc if desc not in ("-", "") else "",
                                ("PFAM:" + pf) if pf not in ("-", "") else ""] if x) or "hypothetical protein"
    gene_rows.append(dict(contig=contig, gene_ID=q, gene_category=fam, eggNOG_annotation=ann))
genes = pd.DataFrame(gene_rows)
genes["contig"] = pd.Categorical(genes["contig"], categories=order, ordered=True)
genes = genes.sort_values(["contig", "gene_ID"]).reset_index(drop=True)

# ---- Sheet2: per-contig annotation calls ----
ann = pd.read_csv(f"{FIG}/kp14_annotation.tsv", sep="\t").set_index("representative")
kn = pd.read_csv(f"{GP}/kp14_ani_known.tsv", sep="\t").set_index("representative")
bsum = pd.read_csv(f"{FIG}/kp14_antismash_summary.tsv", sep="\t").set_index("representative")
try:
    bdet = pd.read_csv(f"{FIG}/kp14_antismash_bgc.tsv", sep="\t")
    bmap = {r["representative"]: r.get("best_known_cluster", "") for _, r in bdet.iterrows()}
except Exception:
    bmap = {}

def g(df, rep, col, default=""):
    try:
        v = df.loc[rep, col]
        return "" if pd.isna(v) else v
    except Exception:
        return default

def ynstr(v):
    return "yes" if str(v).lower() in ("true", "1", "yes") else "no"

rows2 = []
for rep in order:
    ctype = g(reps.set_index("representative"), rep, "rep_type", "")
    # Mobility
    mob = g(ann, rep, "predicted_mobility"); rel = g(ann, rep, "relaxase_MOB_type"); mpf = g(ann, rep, "mpf_type"); repl = g(ann, rep, "replicon_type")
    rows2.append(dict(contig=rep, contig_type=ctype, category="Mobility", value=mob or "unknown",
                      why=f"relaxase={rel or 'none'}; MPF/T4SS={mpf or 'none'}; replicon={repl or 'none'} (MOB-suite)"))
    # Circular
    circ = ynstr(g(reps.set_index("representative"), rep, "rep_circular"))
    rows2.append(dict(contig=rep, contig_type=ctype, category="Circular", value=circ,
                      why="circular contig (hifiasm 'C' / geNomad DTR)" if circ == "yes" else "linear contig"))
    # Known (ANI)
    known = ynstr(g(kn, rep, "known_ANI")); bref = g(kn, rep, "best_ref"); dbn = g(kn, rep, "db")
    ani = g(kn, rep, "ANI"); afs = g(kn, rep, "aln_frac_smaller"); rdesc = g(kn, rep, "ref_description")
    if known == "yes":
        why = f"matches {bref} ({rdesc}) in {dbn}: ANI {ani}%, aln frac(smaller) {afs}% (>=95% ANI & >=85%)"
    else:
        why = f"no reference at >=95% ANI & >=85% cov (best {bref or 'NA'} in {dbn}: ANI {ani}%, cov {afs}%)"
    rows2.append(dict(contig=rep, contig_type=ctype, category="Known (ANI)", value=known, why=why))
    # AMR
    amr = g(ann, rep, "amr_genes"); namr = g(ann, rep, "n_amr")
    rows2.append(dict(contig=rep, contig_type=ctype, category="AMR", value=("yes" if str(namr) not in ("", "0", "nan") else "no"),
                      why=f"AMRFinderPlus AMR genes: {amr}" if amr else "no acquired antibiotic-resistance gene"))
    # Toxin
    tox = ynstr(g(ann, rep, "has_toxin")); ntx = g(ann, rep, "n_toxin_genes")
    rows2.append(dict(contig=rep, contig_type=ctype, category="Toxin", value=tox,
                      why=f"{ntx} toxin-antitoxin gene(s) (eggNOG)" if tox == "yes" else "no toxin-antitoxin gene"))
    # Virulence
    vir = ynstr(g(ann, rep, "has_virulence")); vfdb = g(ann, rep, "vfdb_genes"); avir = g(ann, rep, "amrfinder_virulence_genes")
    vwhy = "; ".join(x for x in [f"VFDB:{vfdb}" if vfdb else "", f"AMRFinder:{avir}" if avir else ""] if x)
    rows2.append(dict(contig=rep, contig_type=ctype, category="Virulence", value=vir,
                      why=vwhy if vir == "yes" else "no virulence factor (VFDB/AMRFinder)"))
    # BGC
    nbgc = g(bsum, rep, "n_BGC"); btypes = g(bsum, rep, "BGC_types")
    has_bgc = "yes" if str(nbgc) not in ("", "0", "nan", "0.0") else "no"
    bwhy = f"antiSMASH: {btypes}" + (f"; best MIBiG: {bmap.get(rep,'')}" if bmap.get(rep, "") else "") if has_bgc == "yes" else "no biosynthetic gene cluster (antiSMASH)"
    rows2.append(dict(contig=rep, contig_type=ctype, category="BGC", value=has_bgc, why=bwhy))

sheet2 = pd.DataFrame(rows2, columns=["contig", "contig_type", "category", "value", "why"])

with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    genes[["contig", "gene_ID", "gene_category", "eggNOG_annotation"]].to_excel(xw, sheet_name="genes", index=False)
    sheet2.to_excel(xw, sheet_name="contig_annotation", index=False)

print(f"wrote {OUT}")
print(f"sheet1 genes: {len(genes)} rows, {genes['contig'].nunique()} contigs, "
      f"{genes['gene_category'].nunique()} categories")
print(f"sheet2 contig_annotation: {len(sheet2)} rows ({sheet2['contig'].nunique()} contigs x "
      f"{sheet2['category'].nunique()} categories)")
print("\ncategory gene totals:")
print(genes['gene_category'].value_counts().to_string())
