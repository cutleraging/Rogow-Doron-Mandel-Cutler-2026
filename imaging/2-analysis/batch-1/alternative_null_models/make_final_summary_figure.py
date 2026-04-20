#!/usr/bin/env python3

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


LABEL_MAP = {
    "whole_image": "Whole image",
    "0_30": "0-30 um",
    "30_60": "30-60 um",
    "60_90": "60-90 um",
    "90_120": "90-120 um",
}
ORDER = ["whole_image", "0_30", "30_60", "60_90", "90_120"]


def main() -> None:
    base = Path("expanded_results")
    observed = pd.read_csv(base / "observed_per_image.csv")
    null_summary = pd.read_csv(base / "expanded_null_model_summary.csv")
    bauer = pd.read_csv(base / "bauer_inspired_per_image.csv")
    bauer_summary = pd.read_csv(base / "bauer_inspired_summary.csv")
    cbc = pd.read_csv(base / "cbc_per_image.csv")
    cbc_summary = pd.read_csv(base / "cbc_summary.csv")

    null_panel = null_summary[null_summary["null_model"] == "local_masked_randomization"].copy()
    null_panel["bin_label"] = null_panel["bin"].map(LABEL_MAP)
    observed["bin_label"] = observed["bin"].map(LABEL_MAP)
    bauer["bin_label"] = bauer["bin"].map(LABEL_MAP)
    bauer_summary["bin_label"] = bauer_summary["bin"].map(LABEL_MAP)
    cbc["bin_label"] = cbc["bin"].map(LABEL_MAP)
    cbc_summary["bin_label"] = cbc_summary["bin"].map(LABEL_MAP)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.3))
    x_positions = np.arange(len(ORDER))
    x_labels = [LABEL_MAP[name] for name in ORDER]

    ax = axes[0]
    sns.stripplot(
        data=observed,
        x="bin_label",
        y="coassoc_freq",
        order=x_labels,
        color="#355070",
        alpha=0.5,
        jitter=0.22,
        size=4.5,
        ax=ax,
    )
    ax.errorbar(
        x_positions,
        null_panel["observed_mean"],
        yerr=null_panel["observed_sem"],
        fmt="D",
        color="#b80c09",
        markersize=6,
        capsize=4,
        linewidth=1.4,
        zorder=5,
        label="Observed mean ± SEM",
    )
    ax.errorbar(
        x_positions,
        null_panel["null_mean"],
        yerr=[
            null_panel["null_mean"] - null_panel["null_ci_low"],
            null_panel["null_ci_high"] - null_panel["null_mean"],
        ],
        fmt="s",
        color="#6c757d",
        markersize=5,
        capsize=4,
        linewidth=1.4,
        zorder=4,
        label="Local-mask null ± 95% interval",
    )
    for x_pos, (_, row) in zip(x_positions, null_panel.iterrows()):
        ax.text(
            x_pos,
            max(row["observed_mean"], row["null_ci_high"]) + 0.03,
            f"p={row['p_value_enrichment']:.3g}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="#343a40",
        )
    ax.set_title("A. Conservative null model")
    ax.set_xlabel("Distance from soma")
    ax.set_ylabel("Co-associated puncta / total puncta")
    ax.set_ylim(bottom=-0.01)

    ax = axes[1]
    sns.stripplot(
        data=bauer,
        x="bin_label",
        y="bauer_specific_freq_bin",
        order=x_labels,
        color="#355070",
        alpha=0.5,
        jitter=0.22,
        size=4.5,
        ax=ax,
    )
    ax.errorbar(
        x_positions,
        bauer_summary["mean"],
        yerr=bauer_summary["sem"],
        fmt="D",
        color="#b80c09",
        markersize=6,
        capsize=4,
        linewidth=1.4,
        zorder=5,
    )
    ax.axhline(0, color="#adb5bd", linewidth=1, linestyle="--")
    ax.set_title("B. Bauer 2021-inspired")
    ax.set_xlabel("Distance from soma")
    ax.set_ylabel("Specific puncta / total puncta")

    ax = axes[2]
    sns.stripplot(
        data=cbc,
        x="bin_label",
        y="cbc_mean_score",
        order=x_labels,
        color="#355070",
        alpha=0.5,
        jitter=0.22,
        size=4.5,
        ax=ax,
    )
    ax.errorbar(
        x_positions,
        cbc_summary["mean"],
        yerr=cbc_summary["sem"],
        fmt="D",
        color="#b80c09",
        markersize=6,
        capsize=4,
        linewidth=1.4,
        zorder=5,
    )
    ax.axhline(0, color="#adb5bd", linewidth=1, linestyle="--")
    ax.set_title("C. Malkusch 2012-inspired CBC")
    ax.set_xlabel("Distance from soma")
    ax.set_ylabel("Mean CBC score")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.03))
    fig.suptitle("Final summary of puncta co-association under conservative and literature-inspired analyses", y=1.08, fontsize=15)
    fig.tight_layout()
    fig.savefig(base / "final_summary_figure.png", dpi=320, bbox_inches="tight")
    fig.savefig(base / "final_summary_figure.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
