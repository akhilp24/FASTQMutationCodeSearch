#!/usr/bin/env python3
"""Generate plots from combined telomere CSV.

This script mirrors the notebook plotting flow for:
- Histograms
- Spearman-by-feature scatter plots
- Pairwise mutation/age Spearman heatmap
- 2x2 trendline + Spearman panels
"""

from __future__ import annotations

import argparse
import os
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr


def plot_histograms_from_csv(csv_path: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    data = pd.read_csv(csv_path)
    numeric_cols = [c for c in data.select_dtypes(include=[np.number]).columns if c != "Age"]
    for col in numeric_cols:
        plt.figure(figsize=(8, 5))
        sns.histplot(data[col].dropna(), kde=True)
        plt.title(f"Histogram: {col}")
        plt.tight_layout()
        safe = col.replace("/", "_").replace(" ", "_").replace(">", "to")
        plt.savefig(os.path.join(output_dir, f"{safe}_hist.png"), dpi=220)
        plt.close()


def plot_spearman_with_age(csv_path: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(csv_path)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if "Age" not in numeric_cols:
        print("No numeric Age column found; skipping spearman scatter plots.")
        return
    results = []
    for col in [c for c in numeric_cols if c != "Age"]:
        sub_df = df.dropna(subset=["Age", col])
        if sub_df.shape[0] < 2:
            continue
        corr, pval = stats.spearmanr(sub_df["Age"], sub_df[col])
        results.append({"Column": col, "Spearman_r": corr, "p_value": pval})
        plt.figure(figsize=(8, 6))
        ax = sns.scatterplot(x=sub_df["Age"], y=sub_df[col], alpha=0.65)
        sns.regplot(x=sub_df["Age"], y=sub_df[col], scatter=False, ci=None, line_kws={"color": "red", "linestyle": "--"}, ax=ax)
        ax.set_xlabel("Age (years)")
        ax.set_ylabel(col)
        ax.set_title(f"{col} vs Age | Spearman r={corr:.2f}, p={pval:.2g}")
        plt.tight_layout()
        safe = col.replace("/", "_").replace(" ", "_").replace(">", "to")
        plt.savefig(os.path.join(output_dir, f"{safe}_vs_age.png"), dpi=220)
        plt.close()
    pd.DataFrame(results).to_csv(os.path.join(output_dir, "spearman_results.csv"), index=False)


def plot_pairwise_r_heatmap(csv_path: str, output_dir: str) -> None:
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(csv_path)
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    mutation_cols = [c for c in numeric_cols if "per_1k" in c or "per1k" in c]
    if "Age" in numeric_cols:
        mutation_cols.append("Age")
    mutation_cols = sorted(set(mutation_cols))
    if len(mutation_cols) < 2:
        print("Insufficient numeric mutation columns for pairwise heatmap.")
        return
    corr = df[mutation_cols].corr(method="spearman")
    plt.figure(figsize=(max(10, len(mutation_cols) * 0.5), max(8, len(mutation_cols) * 0.4)))
    sns.heatmap(corr, cmap="coolwarm", center=0, square=True)
    plt.title("Pairwise Spearman Correlation Heatmap")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pairwise_spearman_heatmap.png"), dpi=220)
    plt.close()


def plot_trendlines(data: pd.DataFrame, output_path: str, variables: list[str], titles: list[str]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    for i, (var, title) in enumerate(zip(variables, titles)):
        ax = axes[i]
        plot_data = data.dropna(subset=["Age", var])
        if len(plot_data) > 1:
            corr, pval = spearmanr(plot_data["Age"], plot_data[var])
            sns.scatterplot(data=plot_data, x="Age", y=var, ax=ax, alpha=0.6)
            sns.regplot(data=plot_data, x="Age", y=var, ax=ax, scatter=False, line_kws={"color": "red", "linestyle": "--"})
            ax.set_title(f"{title}\nSpearman r = {corr:.2f}, p = {pval:.2g}", fontweight="bold")
        else:
            ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Age")
        ax.set_ylabel(var)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot all telomere analysis charts from combined CSV.")
    parser.add_argument("--combined-csv", required=True, help="Path to final combined CSV")
    parser.add_argument("--output-dir", default="telomere_plots", help="Directory for plot outputs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    df = pd.read_csv(args.combined_csv)

    plot_histograms_from_csv(args.combined_csv, os.path.join(args.output_dir, "histograms"))
    plot_spearman_with_age(args.combined_csv, os.path.join(args.output_dir, "spearman_plots"))
    plot_pairwise_r_heatmap(args.combined_csv, os.path.join(args.output_dir, "heatmaps"))

    trend_vars = [
        "total_mutations_over_total_g_strand_per_1k",
        "g_strand_T>C_sum_per_1k",
        "g_strand_G>T_sum_per_1k",
        "g_strand_T>G_sum_per_1k",
    ]
    trend_titles = [
        "Total Mutations Normalized",
        "T > C Mutation Rate Normalized",
        "G > T Mutation Rate Normalized",
        "T > G Mutation Rate Normalized",
    ]
    existing_vars, existing_titles = [], []
    for v, t in zip(trend_vars, trend_titles):
        if v in df.columns:
            existing_vars.append(v)
            existing_titles.append(t)
    if len(existing_vars) == 4 and "Age" in df.columns:
        plot_trendlines(
            df,
            os.path.join(args.output_dir, f"trendlines_{ts}.png"),
            existing_vars,
            existing_titles,
        )
    else:
        print("Skipping 2x2 trendline panel: required columns not found.")

    print(f"Plot generation complete: {args.output_dir}")


if __name__ == "__main__":
    main()
