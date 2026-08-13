#!/usr/bin/env python3
"""
evaluate_c1.py — assess MODIFI ECE->host linkage on the C1 mock communities.

Ground truth = SRA-accession prefix on each contig (ERR10042285_3_L -> ERR10042285).
For each ECE, MODIFI's `host_summary.csv` gives the best predicted host; we label the
prediction correct iff the ECE and the predicted host share the SRA prefix.

Adds, on top of benchmark/linkage/evaluate_linkage.py:
  - `ece_class` planted/orphan labeling (orphans = ECEs whose true host is ABSENT;
    any confident assignment of an orphan is a false positive by construction)
  - orphan-aware FPR / TNR, ROC (AUROC) and PR (AUPRC) using orphans as the negative class
  - recall / precision / FDR vs community size
  - strain / species / isolate resolution confusion (needs Cdb.csv + lineage)

Outputs CSVs under simu_meta_dir/C1/eval/ for the plotting scripts.

Usage:
  python evaluate_c1.py                          # all ladder_* communities found
  python evaluate_c1.py --labels ladder_25,ladder_40
  python evaluate_c1.py --orphans orphans.txt    # file of orphan ECE contig names
"""
import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd

sys.path.append("/home/shuaiw/MODIFI/benchmark/linkage")
from evaluate_linkage import get_sra_id, build_step_curve  # noqa: E402

OUT_ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
EVAL_DIR = os.path.join(OUT_ROOT, "eval")
ISOLATION_SUMMARY = "/home/shuaiw/borg/paper/isolation/GTDB_tree/anno/isolation_sample_summary.tsv"
CDB_CSV = "/home/shuaiw/borg/paper/specificity/iso_99_out/data_tables/Cdb.csv"

# MODIFI's shipped default operating point (for the "confirm default" comparison)
DEFAULT_SCORE = 0.5
DEFAULT_SPEC = 0.01


# ------------------------------------------------------------------ taxonomy maps
def load_taxonomy():
    """SRA -> (species, strain_cluster)."""
    df = pd.read_csv(ISOLATION_SUMMARY, sep="\t")
    df["Species"] = df["Lineage"].apply(
        lambda s: next((p[3:] for p in s.split(";") if p.startswith("s__")), None)
        if isinstance(s, str) else None
    )
    sp = df.set_index("Sample")["Species"].to_dict()

    cdb = pd.read_csv(CDB_CSV)
    cdb["Sample"] = cdb["genome"].str.extract(r"^([A-Za-z]+\d+)")
    strain = cdb.drop_duplicates("Sample").set_index("Sample")["secondary_cluster"].to_dict()
    return sp, strain


# ------------------------------------------------------------------ load one community
def host_summary_path(label):
    return os.path.join(OUT_ROOT, label, "modifi", label, "host_summary.csv")


def community_size(label):
    """Genome count from the manifest (fallback: parse trailing int in label)."""
    man = os.path.join(OUT_ROOT, label, f"{label}.manifest.csv")
    if os.path.exists(man):
        return len(pd.read_csv(man))
    digits = "".join(c for c in label if c.isdigit())
    return int(digits) if digits else -1


def mge_list_path(label):
    return os.path.join(OUT_ROOT, label, f"{label}.mge_list.tsv")


def load_community(label, orphan_ids, sp_map, strain_map):
    """Return one row per ECE in the community's MGE list (the correct denominator).

    ECEs absent from host_summary.csv are UNASSIGNED (host=NaN, final_score=NaN) and
    count as false negatives — not silently dropped.
    """
    hs = host_summary_path(label)
    mlp = mge_list_path(label)
    if not os.path.exists(mlp):
        print(f"  [skip] no mge_list for {label} ({mlp})")
        return None
    ece = pd.read_csv(mlp, sep="\t").rename(columns={"seq_name": "MGE"})
    ece = ece[["MGE", "mge_type", "length"]].drop_duplicates("MGE")

    hs_df = pd.read_csv(hs) if os.path.exists(hs) else pd.DataFrame(columns=["MGE"])
    # left-join predictions onto the full ECE list -> unassigned ECEs kept as NaN
    df = ece.merge(hs_df, on="MGE", how="left")

    df["community"] = label
    df["size"] = community_size(label)
    df["assigned"] = df["host"].notna() if "host" in df.columns else False
    df["mge_sra"] = df["MGE"].apply(get_sra_id)
    df["host_sra"] = df["host"].apply(lambda h: get_sra_id(h) if isinstance(h, str) else None)
    df["ece_class"] = np.where(df["MGE"].isin(orphan_ids), "orphan", "planted")
    # planted correctness = assigned to same isolate; orphans never correct (host absent)
    df["is_correct"] = df["assigned"] & (df["mge_sra"] == df["host_sra"]) & (df["ece_class"] == "planted")
    # taxonomic resolution of the (predicted) call, for the confusion breakdown
    df["mge_sp"] = df["mge_sra"].map(sp_map)
    df["host_sp"] = df["host_sra"].map(sp_map)
    df["mge_strain"] = df["mge_sra"].map(strain_map)
    df["host_strain"] = df["host_sra"].map(strain_map)
    # replicate parsing: "ladder_58_rep3" -> setting "ladder_58", rep 3; bare label -> rep 1
    import re as _re
    m = _re.search(r"_rep(\d+)$", label)
    df["rep"] = int(m.group(1)) if m else 1
    df["setting"] = _re.sub(r"_rep\d+$", "", label)
    return df


def resolution(row):
    if row["ece_class"] == "orphan":
        return "orphan"
    if not row.get("assigned", False):
        return "unassigned"
    if row["mge_sra"] == row["host_sra"]:
        return "correct_isolate"
    if pd.notna(row["mge_strain"]) and row["mge_strain"] == row["host_strain"]:
        return "correct_strain"
    if pd.notna(row["mge_sp"]) and row["mge_sp"] == row["host_sp"]:
        return "correct_species_wrong_strain"
    return "wrong_species"


# ------------------------------------------------------------------ metrics
def confusion_at(df, min_score, min_spec):
    """Apply an operating point; return counts split by ece_class."""
    called = df[(df["final_score"] > min_score) & (df["specificity"] < min_spec)]
    planted = df[df["ece_class"] == "planted"]
    orphan = df[df["ece_class"] == "orphan"]
    called_planted = called[called["ece_class"] == "planted"]
    called_orphan = called[called["ece_class"] == "orphan"]

    tp = int(called_planted["is_correct"].sum())
    fp_mislink = int((~called_planted["is_correct"]).sum())
    fp_orphan = len(called_orphan)
    n_planted = len(planted)
    n_orphan = len(orphan)

    recall = tp / n_planted if n_planted else np.nan
    precision = tp / (tp + fp_mislink + fp_orphan) if (tp + fp_mislink + fp_orphan) else np.nan
    orphan_fpr = fp_orphan / n_orphan if n_orphan else np.nan
    fdr = (fp_mislink + fp_orphan) / (tp + fp_mislink + fp_orphan) if (tp + fp_mislink + fp_orphan) else np.nan
    return dict(min_score=min_score, min_spec=min_spec, tp=tp, fp_mislink=fp_mislink,
                fp_orphan=fp_orphan, n_planted=n_planted, n_orphan=n_orphan,
                recall=recall, precision=precision, orphan_fpr=orphan_fpr, fdr=fdr)


def roc_pr(df):
    """ROC/PR sweeping final_score, positives = correct planted, negatives = everything
    else (mislinks + orphans). Returns arrays and AUCs (only if both classes present)."""
    d = df.dropna(subset=["final_score"]).sort_values("final_score", ascending=False)
    y = d["is_correct"].astype(int).values
    P = int(y.sum())
    N = int((y == 0).sum())
    if P == 0 or N == 0:
        return None
    tp = fp = 0
    rows = []
    for yi, s in zip(y, d["final_score"].values):
        if yi:
            tp += 1
        else:
            fp += 1
        rows.append((s, tp / P, fp / N, tp / (tp + fp)))  # score, tpr(=recall), fpr, precision
    arr = pd.DataFrame(rows, columns=["score", "tpr", "fpr", "precision"])
    trapz = getattr(np, "trapezoid", None) or np.trapz  # np.trapz removed in numpy 2.0
    auroc = float(trapz(arr["tpr"], arr["fpr"]))
    # AUPRC via recall-sorted precision
    ap = arr.sort_values("tpr")
    auprc = float(trapz(ap["precision"], ap["tpr"]))
    return arr, auroc, auprc


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default=None, help="comma-separated community labels")
    ap.add_argument("--orphans", default=None, help="file with orphan ECE contig names (one/line)")
    args = ap.parse_args()

    os.makedirs(EVAL_DIR, exist_ok=True)
    sp_map, strain_map = load_taxonomy()

    orphan_ids = set()
    if args.orphans and os.path.exists(args.orphans):
        orphan_ids = {l.strip() for l in open(args.orphans) if l.strip()}
        print(f"[orphans] {len(orphan_ids)} orphan ECE ids loaded")

    if args.labels:
        labels = args.labels.split(",")
    else:
        labels = sorted(
            os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(p))))
            for p in glob.glob(os.path.join(OUT_ROOT, "*", "modifi", "*", "host_summary.csv"))
        )
    print(f"[eval] communities: {labels}")

    all_pred, curve_rows, size_rows, auc_rows = [], [], [], []
    for label in labels:
        df = load_community(label, orphan_ids, sp_map, strain_map)
        if df is None or df.empty:
            continue
        df["resolution"] = df.apply(resolution, axis=1)
        all_pred.append(df)

        planted = df[df["ece_class"] == "planted"]
        planted_assigned = planted[planted["assigned"]]
        total_planted = len(planted)                       # recall denominator = ALL planted ECEs
        n = df["size"].iloc[0]

        # step curve over final_score; iterate only assigned rows, but recall is out of
        # ALL planted ECEs (unassigned = never called at any threshold)
        if len(planted_assigned):
            fps, recalls, precs, fprs, scores = build_step_curve(planted_assigned, max(total_planted, 1))
            for fp, rec, pr, fpr, sc in zip(fps, recalls, precs, fprs, scores):
                curve_rows.append(dict(community=label, size=n, score_threshold=sc,
                                       fp=fp, fpr=fpr, recall=rec, precision=pr))

        # default operating-point summary
        c = confusion_at(df, DEFAULT_SCORE, DEFAULT_SPEC)
        c.update(community=label, setting=df["setting"].iloc[0], rep=int(df["rep"].iloc[0]),
                 size=n,
                 n_ece=len(df), n_planted=total_planted,
                 n_assigned=int(planted["assigned"].sum()),
                 n_unassigned=int((~planted["assigned"]).sum()),
                 res_correct_isolate=int((df["resolution"] == "correct_isolate").sum()),
                 res_correct_strain=int((df["resolution"] == "correct_strain").sum()),
                 res_correct_sp_wrong_strain=int((df["resolution"] == "correct_species_wrong_strain").sum()),
                 res_wrong_species=int((df["resolution"] == "wrong_species").sum()),
                 res_unassigned=int((df["resolution"] == "unassigned").sum()))
        size_rows.append(c)

        rp = roc_pr(df)
        if rp is not None:
            _, auroc, auprc = rp
            auc_rows.append(dict(community=label, size=n, auroc=auroc, auprc=auprc))

        print(f"  {label}: ECEs={len(df)} planted={len(planted)} orphan={int((df.ece_class=='orphan').sum())} "
              f"TP@default={c['tp']} mislinkFP={c['fp_mislink']} orphanFP={c['fp_orphan']} "
              f"recall={c['recall']:.3f} precision={c['precision']:.3f}")

    if all_pred:
        pd.concat(all_pred, ignore_index=True).to_csv(os.path.join(EVAL_DIR, "predictions.csv"), index=False)
    pd.DataFrame(curve_rows).to_csv(os.path.join(EVAL_DIR, "step_curves.csv"), index=False)
    size_df = pd.DataFrame(size_rows)
    size_df.to_csv(os.path.join(EVAL_DIR, "size_summary.csv"), index=False)
    pd.DataFrame(auc_rows).to_csv(os.path.join(EVAL_DIR, "auc.csv"), index=False)

    # ---- aggregate across replicates: mean +/- 95% CI per setting ----
    if len(size_df):
        agg = aggregate_replicates(size_df)
        agg.to_csv(os.path.join(EVAL_DIR, "size_summary_agg.csv"), index=False)
        print("\n[eval] per-setting mean +/- 95% CI across replicates:")
        show = agg[["setting", "size", "n_rep", "recall_mean", "recall_ci",
                    "precision_mean", "precision_ci"]]
        print(show.to_string(index=False))
    print(f"[eval] wrote CSVs -> {EVAL_DIR}")


def aggregate_replicates(size_df):
    """Mean +/- 95% CI across replicates for each setting (t-based; falls back to
    half-range for n<2)."""
    from scipy import stats
    metrics = ["recall", "precision", "orphan_fpr", "fdr", "tp", "fp_mislink",
               "fp_orphan", "n_assigned", "n_unassigned"]
    out = []
    for setting, g in size_df.groupby("setting"):
        row = {"setting": setting, "size": int(g["size"].iloc[0]), "n_rep": len(g)}
        for m in metrics:
            vals = g[m].dropna().values
            mean = float(np.mean(vals)) if len(vals) else np.nan
            if len(vals) >= 2:
                sem = stats.sem(vals)
                ci = float(sem * stats.t.ppf(0.975, len(vals) - 1))
            else:
                ci = np.nan
            row[f"{m}_mean"] = round(mean, 4) if mean == mean else mean
            row[f"{m}_ci"] = round(ci, 4) if ci == ci else ci
        out.append(row)
    return pd.DataFrame(out).sort_values(["setting", "size"]).reset_index(drop=True)


if __name__ == "__main__":
    main()
