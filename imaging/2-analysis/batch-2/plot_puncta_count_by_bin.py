#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BIN_ORDER = ["0_40", "40_80", "80_120", "120_160"]
BIN_LABELS = {
    "0_40": "0-40 um",
    "40_80": "40-80 um",
    "80_120": "80-120 um",
    "120_160": "120-160 um",
}


def build_threshold_summary(df: pd.DataFrame, max_threshold: int) -> pd.DataFrame:
    rows: list[dict[str, int | str | float]] = []
    for bin_name in BIN_ORDER:
        subset = df[df["bin"] == bin_name].copy()
        n_images = int(len(subset))
        for threshold in range(1, max_threshold + 1):
            kept = int((subset["total_puncta"] >= threshold).sum())
            rows.append(
                {
                    "bin": bin_name,
                    "min_total_puncta": threshold,
                    "images_kept": kept,
                    "fraction_kept": kept / n_images if n_images else np.nan,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observed", default="results_250nm/observed_per_image.csv")
    parser.add_argument("--outdir", default="results_250nm")
    parser.add_argument("--max-threshold", type=int, default=15)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    observed = pd.read_csv(args.observed)
    observed = observed[observed["bin"].isin(BIN_ORDER)].copy()
    observed["bin"] = pd.Categorical(observed["bin"], categories=BIN_ORDER, ordered=True)
    observed["bin_label"] = observed["bin"].map(BIN_LABELS)
    observed["actin_fraction"] = observed["actin_n"] / observed["total_puncta"].replace(0, np.nan)
    observed["arc_fraction"] = observed["arc_n"] / observed["total_puncta"].replace(0, np.nan)

    threshold_summary = build_threshold_summary(observed, args.max_threshold)
    threshold_summary["bin_label"] = threshold_summary["bin"].map(BIN_LABELS)
    threshold_summary.to_csv(outdir / "puncta_count_threshold_summary.csv", index=False)

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.linewidth": 0.8,
        }
    )
    sns.set_theme(style="white", rc={"axes.grid": False})

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), gridspec_kw={"width_ratios": [1.15, 1.0]})

    ax = axes[0]
    sns.boxplot(
        data=observed,
        x="bin_label",
        y="total_puncta",
        order=[BIN_LABELS[b] for b in BIN_ORDER],
        color="#efefef",
        width=0.55,
        showcaps=True,
        fliersize=0,
        linewidth=0.9,
        ax=ax,
    )
    sns.stripplot(
        data=observed,
        x="bin_label",
        y="total_puncta",
        order=[BIN_LABELS[b] for b in BIN_ORDER],
        color="#4d4d4d",
        size=4,
        jitter=0.23,
        alpha=0.75,
        ax=ax,
    )
    medians = observed.groupby("bin", observed=False)["total_puncta"].median().reindex(BIN_ORDER)
    for idx, (bin_name, median_value) in enumerate(medians.items()):
        ax.text(
            idx,
            median_value + 1.0,
            f"med={int(median_value)}",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#4d4d4d",
        )
    ax.set_xlabel("Distance from soma")
    ax.set_ylabel("Total puncta per image")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.8)

    ax = axes[1]
    sns.lineplot(
        data=threshold_summary,
        x="min_total_puncta",
        y="fraction_kept",
        hue="bin_label",
        hue_order=[BIN_LABELS[b] for b in BIN_ORDER],
        palette=["#1b9e77", "#7570b3", "#d95f02", "#e7298a"],
        marker="o",
        linewidth=1.6,
        markersize=4.5,
        ax=ax,
    )
    ax.set_xlabel("Minimum total puncta required in bin")
    ax.set_ylabel("Fraction of images retained")
    ax.set_ylim(-0.02, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.8)
    ax.legend(frameon=False, title="", fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(outdir / "puncta_count_by_bin.png", dpi=350, bbox_inches="tight")
    fig.savefig(outdir / "puncta_count_by_bin.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
