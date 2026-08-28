#!/usr/bin/env python3
"""Combined 9-panel profile of the FINAL strict ECE set (3,880), plasmid vs virus:
 a Venn (5 callers)  b length dist  c depth dist  d environment  e linkage-score dist
 f GC similarity (MGE vs host)  g depth similarity  h tetranucleotide cosine  i CRISPR consistent.
Reads expanded/filterpass_FINAL.csv + final_profile/{ece_table,linkage_table,crispr_by_sample}.csv."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from venn import venn

A = "/home/shuaiw/borg/revision/ece_anno"
FP = f"{A}/expanded/final_profile"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
CP, CV = "#F8766D", "#00BFC4"   # plasmid (salmon), virus (teal) - match manuscript
pal = {"plasmid": CP, "virus": CV}
TYPES = ["plasmid", "virus"]

ece = pd.read_csv(f"{FP}/ece_table.csv")
for c in ("length", "depth", "best_score"): ece[c] = pd.to_numeric(ece[c], errors="coerce")
lk = pd.read_csv(f"{FP}/linkage_table.csv")
for c in ("MGE_gc", "host_gc", "cos_sim", "MGE_cov", "host_cov"): lk[c] = pd.to_numeric(lk[c], errors="coerce")
fp = pd.read_csv(f"{A}/expanded/filterpass_FINAL.csv")
fp["key"] = fp["sample"] + "|" + fp["MGE"].astype(str)

plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold"})
fig, axes = plt.subplots(4, 3, figsize=(18, 21))
def title(ax, l, t): ax.set_title(f"{l}. {t}", loc="left")

# a) Venn
m = fp["methods"].fillna("")
sets = {"geNomad": set(fp.loc[m.str.contains("genomad"),"key"]),
        "VirSorter1": set(fp.loc[m.str.contains("virsorter1"),"key"]),
        "VirSorter2": set(fp.loc[m.str.contains("virsorter2"),"key"]),
        "VIBRANT": set(fp.loc[m.str.contains("vibrant"),"key"]),
        "PlasX": set(fp.loc[m.str.contains("plasx"),"key"])}
venn(sets, cmap=["#4C78A8","#F58518","#54A24B","#B279A2","#E45756"], fontsize=7, legend_loc="upper left", ax=axes[0,0])
title(axes[0,0], "a", f"Callers (n={len(fp)})")

# b) length
ax=axes[0,1]; title(ax,"b","ECE length")
for t in TYPES:
    v=ece.loc[ece.type==t,"length"].dropna()
    ax.hist(v, bins=np.logspace(np.log10(5000), np.log10(max(v.max(),1e6)), 40), histtype="step", lw=2, color=pal[t], label=f"{t} (n={len(v)})")
ax.set_xscale("log"); ax.set_xlabel("length (bp)"); ax.set_ylabel("ECEs"); ax.legend()

# c) depth
ax=axes[0,2]; title(ax,"c","ECE depth")
for t in TYPES:
    v=ece.loc[ece.type==t,"depth"].dropna(); v=v[v>0]
    ax.hist(v, bins=np.logspace(np.log10(max(v.min(),1)), np.log10(v.max()), 40), histtype="step", lw=2, color=pal[t], label=t)
ax.set_xscale("log"); ax.set_xlabel("mean depth (x)"); ax.set_ylabel("ECEs"); ax.legend()

# d) environment
ax=axes[1,0]; title(ax,"d","Environment")
ct=ece.groupby(["environment","type"]).size().unstack(fill_value=0).reindex(columns=TYPES, fill_value=0)
ct["s"]=ct.sum(1); ct=ct.sort_values("s",ascending=False).drop(columns="s")
x=np.arange(len(ct)); w=0.4
ax.bar(x-w/2, ct["plasmid"], w, color=CP, label="plasmid"); ax.bar(x+w/2, ct["virus"], w, color=CV, label="virus")
ax.set_xticks(x); ax.set_xticklabels(ct.index, rotation=45, ha="right"); ax.set_ylabel("ECEs"); ax.legend()

# e) linkage score distribution
ax=axes[1,1]; title(ax,"e","Best linkage score")
for t in TYPES:
    v=ece.loc[ece.type==t,"best_score"].dropna()
    ax.hist(v, bins=np.linspace(0,1,26), histtype="step", lw=2, color=pal[t], label=t)
ax.axvline(0.8, ls="--", c="grey", lw=1); ax.set_yscale("log"); ax.set_xlabel("final_score (best host)"); ax.set_ylabel("ECEs"); ax.legend()

# f) GC similarity
def scat(ax,l,t,x,y,logg):
    title(ax,l,t); d=lk.dropna(subset=[x,y]).copy()
    if logg: d=d[(d[x]>0)&(d[y]>0)]
    for tp in TYPES:
        s=d[d.type==tp]; ax.scatter(s[x],s[y],s=22,c=pal[tp],alpha=0.7,marker=('o' if tp=='plasmid' else 'x'),linewidths=(0 if tp=='plasmid' else 1),label=tp)
    xv,yv=d[x].values,d[y].values
    if logg:
        r,p=pearsonr(np.log10(xv),np.log10(yv)); ax.set_xscale("log"); ax.set_yscale("log")
        b,a0=np.polyfit(np.log10(xv),np.log10(yv),1); xs=np.linspace(np.log10(xv.min()),np.log10(xv.max()),50); ax.plot(10**xs,10**(a0+b*xs),"k-",lw=1.3)
    else:
        r,p=pearsonr(xv,yv); b,a0=np.polyfit(xv,yv,1); xs=np.linspace(xv.min(),xv.max(),50); ax.plot(xs,a0+b*xs,"k-",lw=1.3)
    ax.text(0.05,0.93,f"r={r:.2f}\np={p:.1e}",transform=ax.transAxes,va="top",fontsize=9); ax.legend(fontsize=8)
scat(axes[1,2],"f","GC similarity (ECE vs host)","MGE_gc","host_gc",False); axes[1,2].set_xlabel("ECE GC"); axes[1,2].set_ylabel("host GC")
scat(axes[2,0],"g","Depth similarity (ECE vs host)","MGE_cov","host_cov",True); axes[2,0].set_xlabel("ECE cov (log)"); axes[2,0].set_ylabel("host cov (log)")

# h) cosine similarity by habitat (ordered by median), colored by habitat
ax=axes[2,1]; title(ax,"h","Tetranucleotide cosine by habitat")
d=lk.dropna(subset=["cos_sim"])
order=d.groupby("environment")["cos_sim"].median().sort_values(ascending=False).index.tolist()
hcmap=plt.get_cmap("tab10"); hcol={h:hcmap(i%10) for i,h in enumerate(order)}
data=[d[d.environment==h]["cos_sim"].values for h in order]
bp=ax.boxplot(data, patch_artist=True, showfliers=True, flierprops=dict(marker=".",markersize=3))
for patch,h in zip(bp["boxes"],order): patch.set_facecolor(hcol[h])
for med in bp["medians"]: med.set_color("black")
ax.set_xticks(range(1,len(order)+1)); ax.set_xticklabels(order, rotation=45, ha="right"); ax.set_ylabel("Cosine similarity")

# i) CRISPR
ax=axes[2,2]; title(ax,"i","CRISPR-consistent linkages")
cf=f"{FP}/crispr_by_sample.csv"
if os.path.exists(cf) and os.path.getsize(cf)>30:
    c=pd.read_csv(cf); piv=c.pivot_table(index="Sample",columns="MGE Type",values="Consistent Linkages",aggfunc="sum",fill_value=0)
    for t in TYPES:
        if t not in piv.columns: piv[t]=0
    piv["s"]=piv["plasmid"]+piv["virus"]; piv=piv.sort_values("s",ascending=False).drop(columns="s")
    x=np.arange(len(piv)); ax.bar(x,piv["plasmid"],color=CP,label="plasmid"); ax.bar(x,piv["virus"],bottom=piv["plasmid"],color=CV,label="virus")
    ax.set_xticks(x); ax.set_xticklabels(piv.index,rotation=90,fontsize=7); ax.set_ylabel("consistent linkages"); ax.legend()

# j) IMG-catalogue host validation (full-width banner spanning the bottom row; KNOWN = ANI>=95 & qcov>=85)
for _c in range(3): axes[3,_c].remove()
axj=fig.add_subplot(axes[0,0].get_gridspec()[3, :])
title(axj,"j","IMG-catalogue host validation")
vsum=f"{OUT}/ece_validation_summary_strict.csv"
if os.path.exists(vsum):
    vs=pd.read_csv(vsum)
    catlev=["host-supported (agrees)","host mismatch","no comparable reference"]
    vpal={"host-supported (agrees)":"#238b45","host mismatch":"#d7301f","no comparable reference":"#e0e0e0"}
    vtypes=["Plasmids (IMG/PR)","Viruses (IMG/VR)"]
    piv=vs.pivot_table(index="type",columns="category",values="n",fill_value=0).reindex(vtypes)
    # grey = no strict IMG hit OR a strict hit whose reference lacks a genus-level host
    piv["no comparable reference"]=piv.get("no strict match",0)+piv.get("known (strict), host not genus-resolved",0)
    y=np.arange(len(vtypes)); left=np.zeros(len(vtypes))
    for c in catlev:
        vals=piv[c].values if c in piv.columns else np.zeros(len(vtypes))
        axj.barh(y,vals,left=left,color=vpal[c],label=c,edgecolor="white",linewidth=0.8,height=0.62)
        for yi,(v,l) in enumerate(zip(vals,left)):
            if v>=3: axj.text(l+v/2,yi,int(v),ha="center",va="center",fontsize=10,color="#262626")
        left=left+vals
    labs=[]
    for t in vtypes:
        row=vs[vs.type==t]
        total=row.n.sum(); known=row[row.category!="no strict match"].n.sum()
        resolved=row[row.category.isin(["host-supported (agrees)","host mismatch"])].n.sum()
        supp=row[row.category=="host-supported (agrees)"].n.sum()
        labs.append(t)
    axj.set_yticks(y); axj.set_yticklabels(labs, fontsize=11)
    axj.set_xlabel("linked ECEs", fontsize=11); axj.set_ylim(-0.65, len(vtypes)-0.35)
    axj.margins(x=0); axj.set_xlim(0, left.max()*1.02)
    for sp in ("top","right","left"): axj.spines[sp].set_visible(False)
    axj.tick_params(axis="y", length=0)
    axj.legend(fontsize=9, loc="upper center", bbox_to_anchor=(0.5,-0.18), ncol=4, frameon=False,
               handlelength=1.2, columnspacing=1.6)

fig.suptitle(f"Metagenome ECE set (strict criteria): {len(fp)} ECEs ({(fp.MGE_type=='plasmid').sum()} plasmid, {(fp.MGE_type=='virus').sum()} virus)", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.98])
# give panel j (full-width banner) a real left margin so its 3-line y-labels are not clipped
try:
    _p=axj.get_position(); axj.set_position([0.16, _p.y0, 0.80, _p.height])
except Exception: pass
for ext in ("pdf","png"): fig.savefig(f"{OUT}/ece_profile_final.{ext}", bbox_inches="tight", dpi=200)
# sourcedata
ece.to_csv(f"{OUT}/ece_profile_final_sourcedata_ece.csv", index=False)
lk.to_csv(f"{OUT}/ece_profile_final_sourcedata_linkage.csv", index=False)
print("wrote ece_profile_final.pdf/.png + sourcedata")
