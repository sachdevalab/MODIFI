import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

def get_sra_id(name):
    """Extract SRA ID (e.g. ERR10042297) from contig name like ERR10042297_2_L."""
    return name.split("_")[0]

def evaluate_at_cutoff(df, min_score, total_tp):
    subset = df[df["final_score"] > min_score]
    tp = subset["is_correct"].sum()
    fp = len(subset) - tp
    recall = tp / total_tp if total_tp > 0 else 0
    precision = tp / len(subset) if len(subset) > 0 else 0
    return int(tp), int(fp), recall, precision

def build_step_curve(df, total_tp):
    """Sort by final_score desc, accumulate TP/FP to get a monotone step function."""
    subset = df.sort_values("final_score", ascending=False)
    fps, recalls, precisions, scores = [], [], [], []
    tp = fp = 0
    for _, row in subset.iterrows():
        if row["is_correct"]:
            tp += 1
        else:
            fp += 1
        fps.append(fp)
        recalls.append(tp / total_tp)
        precisions.append(tp / (tp + fp))
        scores.append(row["final_score"])
    # prepend origin
    return (
        [0] + fps,
        [0.0] + recalls,
        [1.0] + precisions,
        [None] + scores,
    )

def main(csv_path):
    df = pd.read_csv(csv_path)
    df["mge_sra"] = df["MGE"].apply(get_sra_id)
    df["host_sra"] = df["host"].apply(get_sra_id)
    df["is_correct"] = df["mge_sra"] == df["host_sra"]

    total_tp = int(df["is_correct"].sum())
    total = len(df)

    print(f"Total predictions: {total}")
    print(f"True pairs (correct):    {total_tp}")
    print(f"False pairs (incorrect): {total - total_tp}")
    print()

    fp_rows = df[~df["is_correct"]]
    if not fp_rows.empty:
        print("False positive linkages:")
        print(fp_rows[["MGE", "host", "final_score", "specificity"]].to_string(index=False))
        print()

    # Named cutoffs (score only)
    named_cutoffs = [
        (0.68, "score > 0.68\n(0 FP)"),   # above FP1 → no FPs, recall=44%
        (0.40, "score > 0.40\n(1 FP)"),   # between FP1 and FP2 → 1 FP, recall=92%
        (0.0,  "score > 0\n(all)"),        # all predictions → 2 FPs, recall=100%
    ]

    print(f"{'Cutoff':<28} {'TP':>4} {'FP':>4} {'Recall':>8} {'Precision':>10}")
    print("-" * 58)
    for min_score, label in named_cutoffs:
        tp, fp, recall, prec = evaluate_at_cutoff(df, min_score, total_tp)
        print(f"{label.replace(chr(10), ' '):<28} {tp:>4} {fp:>4} {recall:>8.3f} {prec:>10.3f}")

    # ── Save CSV ──────────────────────────────────────────────────────────
    out_dir = "/home/shuaiw/MODIFI/tmp/figures/link_accuracy"
    os.makedirs(out_dir, exist_ok=True)

    fps, recalls, precs, scores = build_step_curve(df, total_tp)
    curve_csv = pd.DataFrame({
        "score_threshold": scores,
        "fp": fps,
        "recall": recalls,
        "precision": precs,
    })
    curve_csv_path = os.path.join(out_dir, "linkage_curve.csv")
    curve_csv.to_csv(curve_csv_path, index=False)
    print(f"Curve data saved to {curve_csv_path}")

    cutoff_rows = []
    for min_score, label in named_cutoffs:
        tp, fp, recall, prec = evaluate_at_cutoff(df, min_score, total_tp)
        cutoff_rows.append({
            "label": label.replace("\n", " "),
            "min_score": min_score,
            "tp": tp,
            "fp": fp,
            "recall": recall,
            "precision": prec,
        })
    cutoff_csv = pd.DataFrame(cutoff_rows)
    cutoff_csv_path = os.path.join(out_dir, "linkage_cutoffs.csv")
    cutoff_csv.to_csv(cutoff_csv_path, index=False)
    print(f"Cutoff data saved to {cutoff_csv_path}")

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else \
        "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/mix_16/modifi/mix_16/host_summary.csv"
    main(path)
