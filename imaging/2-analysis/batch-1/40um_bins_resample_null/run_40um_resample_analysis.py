#!/usr/bin/env python3

import argparse
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BINS = [
    ("whole_image", None, None),
    ("0_40", 0.0, 40.0),
    ("40_80", 40.0, 80.0),
    ("80_120", 80.0, 120.0),
    ("120_160", 120.0, 160.0),
]


@dataclass(frozen=True)
class BinSpec:
    name: str
    low: float | None
    high: float | None


def canonicalize_source_file(name: str) -> str:
    return re.sub(r"^C\d+-", "", name)


def load_puncta_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(
        columns={
            "Source File": "source_file",
            "X (um)": "x_um",
            "Y (um)": "y_um",
        }
    )
    df["image_id"] = df["source_file"].map(canonicalize_source_file)
    return df[["image_id", "source_file", "x_um", "y_um"]]


def subset_points(points: np.ndarray, bin_spec: BinSpec) -> np.ndarray:
    if bin_spec.low is None:
        return points
    mask = (points[:, 0] >= bin_spec.low) & (points[:, 0] < bin_spec.high)
    return points[mask]


def overlap_matrix(actin_points: np.ndarray, arc_points: np.ndarray, threshold_um: float) -> np.ndarray:
    if len(actin_points) == 0 or len(arc_points) == 0:
        return np.zeros((len(actin_points), len(arc_points)), dtype=bool)
    deltas = actin_points[:, None, :] - arc_points[None, :, :]
    return np.sum(deltas * deltas, axis=2) <= (threshold_um ** 2)


def compute_coassociation(actin_points: np.ndarray, arc_points: np.ndarray, threshold_um: float) -> dict[str, float]:
    within_threshold = overlap_matrix(actin_points, arc_points, threshold_um)
    actin_n = int(len(actin_points))
    arc_n = int(len(arc_points))
    actin_assoc_n = int(within_threshold.any(axis=1).sum()) if actin_n else 0
    arc_assoc_n = int(within_threshold.any(axis=0).sum()) if arc_n else 0
    total_puncta = actin_n + arc_n
    total_assoc = actin_assoc_n + arc_assoc_n
    return {
        "actin_n": actin_n,
        "arc_n": arc_n,
        "total_puncta": total_puncta,
        "actin_assoc_n": actin_assoc_n,
        "arc_assoc_n": arc_assoc_n,
        "total_assoc_puncta": total_assoc,
        "pair_count": int(within_threshold.sum()),
        "coassoc_freq": (total_assoc / total_puncta) if total_puncta else 0.0,
    }


def resample_with_replacement(points: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    indices = rng.integers(0, len(points), size=len(points))
    return points[indices]


def summarize_with_null(observed: np.ndarray, null_distribution: np.ndarray) -> dict[str, float]:
    observed_mean = float(np.mean(observed)) if len(observed) else np.nan
    observed_sd = float(np.std(observed, ddof=1)) if len(observed) > 1 else np.nan
    observed_sem = observed_sd / np.sqrt(len(observed)) if len(observed) > 1 else np.nan
    if len(null_distribution) == 0:
        return {
            "observed_mean": observed_mean,
            "observed_sd": observed_sd,
            "observed_sem": observed_sem,
            "null_mean": np.nan,
            "null_sd": np.nan,
            "delta_vs_null": np.nan,
            "z_score": np.nan,
            "p_value_enrichment": np.nan,
            "null_ci_low": np.nan,
            "null_ci_high": np.nan,
            "n_images": int(len(observed)),
        }

    null_mean = float(np.mean(null_distribution))
    null_sd = float(np.std(null_distribution, ddof=1)) if len(null_distribution) > 1 else 0.0
    p_value = float((1 + np.sum(null_distribution >= observed_mean)) / (len(null_distribution) + 1))
    z_score = (observed_mean - null_mean) / null_sd if null_sd > 0 else np.nan
    return {
        "observed_mean": observed_mean,
        "observed_sd": observed_sd,
        "observed_sem": observed_sem,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "delta_vs_null": observed_mean - null_mean,
        "z_score": z_score,
        "p_value_enrichment": p_value,
        "null_ci_low": float(np.quantile(null_distribution, 0.025)),
        "null_ci_high": float(np.quantile(null_distribution, 0.975)),
        "n_images": int(len(observed)),
    }


def summarize_replicates(values: np.ndarray) -> dict[str, float]:
    mean = float(np.mean(values)) if len(values) else np.nan
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
    sem = sd / np.sqrt(len(values)) if len(values) > 1 else np.nan
    return {"mean": mean, "sd": sd, "sem": sem, "n_images": int(len(values))}


def significance_stars(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "ns"
    if p_value <= 0.0001:
        return "****"
    if p_value <= 0.001:
        return "***"
    if p_value <= 0.01:
        return "**"
    if p_value <= 0.05:
        return "*"
    return "ns"


def plot_publication_summary(
    observed_df: pd.DataFrame,
    null_mean_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_png: Path,
    output_pdf: Path,
) -> None:
    label_map = {
        "whole_image": "Whole image",
        "0_40": "0-40 um",
        "40_80": "40-80 um",
        "80_120": "80-120 um",
        "120_160": "120-160 um",
    }
    order = [name for name, _, _ in BINS if name != "whole_image"]
    x_labels = [label_map[name] for name in order]

    observed_df = observed_df.copy()
    null_mean_df = null_mean_df.copy()
    summary_df = summary_df[summary_df["bin"].isin(order)].copy()
    observed_df = observed_df[observed_df["bin"].isin(order)].copy()
    null_mean_df = null_mean_df[null_mean_df["bin"].isin(order)].copy()
    observed_df["plot_value"] = observed_df["coassoc_freq"] * 100.0
    null_mean_df["plot_value"] = null_mean_df["null_mean_coassoc_freq"] * 100.0

    observed_stats = []
    null_stats = []
    for bin_name in order:
        observed_values = observed_df.loc[observed_df["bin"] == bin_name, "coassoc_freq"].to_numpy(dtype=float)
        null_values = null_mean_df.loc[null_mean_df["bin"] == bin_name, "null_mean_coassoc_freq"].to_numpy(dtype=float)
        observed_stats.append({"bin": bin_name, **summarize_replicates(observed_values)})
        null_stats.append({"bin": bin_name, **summarize_replicates(null_values)})
    observed_stats_df = pd.DataFrame(observed_stats)
    null_stats_df = pd.DataFrame(null_stats)

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
    fig, ax = plt.subplots(figsize=(7.2, 4.4))

    x_positions = np.arange(len(order))
    width = 0.34
    observed_x = x_positions - width / 2
    null_x = x_positions + width / 2

    ax.bar(
        observed_x,
        observed_stats_df["mean"] * 100.0,
        width=width,
        color="#bfbfbf",
        edgecolor="#4d4d4d",
        linewidth=0.8,
        zorder=2,
        label="Observed",
    )
    ax.bar(
        null_x,
        null_stats_df["mean"] * 100.0,
        width=width,
        color="#ededed",
        edgecolor="#7a7a7a",
        linewidth=0.8,
        zorder=2,
        label="Null",
    )

    ax.errorbar(
        observed_x,
        observed_stats_df["mean"] * 100.0,
        yerr=observed_stats_df["sem"] * 100.0,
        fmt="none",
        ecolor="#4d4d4d",
        elinewidth=0.9,
        capsize=2.5,
        zorder=3,
    )
    ax.errorbar(
        null_x,
        null_stats_df["mean"] * 100.0,
        yerr=null_stats_df["sem"] * 100.0,
        fmt="none",
        ecolor="#7a7a7a",
        elinewidth=0.9,
        capsize=2.5,
        zorder=3,
    )

    rng = np.random.default_rng(20260320)
    for idx, bin_name in enumerate(order):
        observed_points = observed_df.loc[observed_df["bin"] == bin_name, "plot_value"].to_numpy(dtype=float)
        null_points = null_mean_df.loc[null_mean_df["bin"] == bin_name, "plot_value"].to_numpy(dtype=float)
        obs_jitter = rng.uniform(-width * 0.22, width * 0.22, size=len(observed_points))
        null_jitter = rng.uniform(-width * 0.22, width * 0.22, size=len(null_points))
        ax.scatter(
            np.full(len(observed_points), observed_x[idx]) + obs_jitter,
            observed_points,
            s=12,
            color="#4d4d4d",
            alpha=0.75,
            zorder=4,
        )
        ax.scatter(
            np.full(len(null_points), null_x[idx]) + null_jitter,
            null_points,
            s=12,
            color="#8c8c8c",
            alpha=0.75,
            zorder=4,
        )

    bar_tops = np.maximum(
        observed_stats_df["mean"].to_numpy(dtype=float) * 100.0
        + np.nan_to_num(observed_stats_df["sem"].to_numpy(dtype=float) * 100.0, nan=0.0),
        null_stats_df["mean"].to_numpy(dtype=float) * 100.0
        + np.nan_to_num(null_stats_df["sem"].to_numpy(dtype=float) * 100.0, nan=0.0),
    )
    point_tops = np.maximum(
        observed_df.groupby("bin")["plot_value"].max().reindex(order).fillna(0.0).to_numpy(dtype=float),
        null_mean_df.groupby("bin")["plot_value"].max().reindex(order).fillna(0.0).to_numpy(dtype=float),
    )
    annotation_base = np.maximum(bar_tops, point_tops)
    for idx, (_, row) in enumerate(summary_df.iterrows()):
        y = annotation_base[idx] + 1.6
        h = 0.7
        ax.plot(
            [observed_x[idx], observed_x[idx], null_x[idx], null_x[idx]],
            [y, y + h, y + h, y],
            color="#4d4d4d",
            linewidth=0.8,
            zorder=5,
        )
        ax.text(
            x_positions[idx],
            y + h + 0.25,
            significance_stars(float(row["p_value_enrichment"])),
            ha="center",
            va="bottom",
            fontsize=9,
            color="#4d4d4d",
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel("Distance from soma")
    ax.set_ylabel("Co-association freq (%)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.8)
    ax.legend(frameon=False, loc="upper right", fontsize=9, handlelength=1.2)
    ymax = float(np.max(annotation_base + 4.0)) if len(annotation_base) else 8.0
    ax.set_ylim(0, max(8.0, ymax))
    fig.tight_layout()
    fig.savefig(output_png, dpi=400, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actin", default="../actin puncta locations.csv")
    parser.add_argument("--arc", default="../arc puncta locations.csv")
    parser.add_argument("--threshold-um", type=float, default=0.1)
    parser.add_argument("--iterations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260318)
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.plot_only:
        plot_publication_summary(
            pd.read_csv(outdir / "observed_per_image.csv"),
            pd.read_csv(outdir / "null_mean_per_image.csv"),
            pd.read_csv(outdir / "summary.csv"),
            outdir / "coassociation_40um_bins_publication.png",
            outdir / "coassociation_40um_bins_publication.pdf",
        )
        return

    actin_df = load_puncta_table(args.actin)
    arc_df = load_puncta_table(args.arc)
    shared_images = sorted(set(actin_df["image_id"]) & set(arc_df["image_id"]))
    if not shared_images:
        raise SystemExit("No shared images after canonicalizing source file names.")

    rng = np.random.default_rng(args.seed)
    bin_specs = [BinSpec(*bin_spec) for bin_spec in BINS]

    actin_by_image = {
        image_id: actin_df.loc[actin_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }
    arc_by_image = {
        image_id: arc_df.loc[arc_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }

    observed_rows: list[dict[str, float | str | int]] = []
    null_mean_rows: list[dict[str, float | str | int]] = []
    observed_metric_values = {bin_spec.name: [] for bin_spec in bin_specs}
    null_metric_values = {bin_spec.name: [] for bin_spec in bin_specs}

    for image_id in shared_images:
        actin_points = actin_by_image[image_id]
        arc_points = arc_by_image[image_id]

        for bin_spec in bin_specs:
            observed = compute_coassociation(
                subset_points(actin_points, bin_spec),
                subset_points(arc_points, bin_spec),
                args.threshold_um,
            )
            observed_rows.append({"image_id": image_id, "bin": bin_spec.name, **observed})
            observed_metric_values[bin_spec.name].append(observed["coassoc_freq"])

        image_null = {bin_spec.name: [] for bin_spec in bin_specs}
        for _ in range(args.iterations):
            null_actin = resample_with_replacement(actin_points, rng)
            null_arc = resample_with_replacement(arc_points, rng)
            for bin_spec in bin_specs:
                null_metrics = compute_coassociation(
                    subset_points(null_actin, bin_spec),
                    subset_points(null_arc, bin_spec),
                    args.threshold_um,
                )
                image_null[bin_spec.name].append(null_metrics["coassoc_freq"])

        for bin_name, values in image_null.items():
            null_metric_values[bin_name].append(values)
            null_mean_rows.append(
                {
                    "image_id": image_id,
                    "bin": bin_name,
                    "null_mean_coassoc_freq": float(np.mean(values)) if len(values) else np.nan,
                }
            )

    observed_df = pd.DataFrame(observed_rows)
    null_mean_df = pd.DataFrame(null_mean_rows)
    observed_df.to_csv(outdir / "observed_per_image.csv", index=False)
    null_mean_df.to_csv(outdir / "null_mean_per_image.csv", index=False)

    summary_rows = []
    for bin_spec in bin_specs:
        observed = np.asarray(observed_metric_values[bin_spec.name], dtype=float)
        per_image_null = np.asarray(null_metric_values[bin_spec.name], dtype=float)
        group_null = np.mean(per_image_null, axis=0) if per_image_null.size else np.asarray([])
        summary_rows.append({"bin": bin_spec.name, **summarize_with_null(observed, group_null)})

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "summary.csv", index=False)
    plot_publication_summary(
        observed_df,
        null_mean_df,
        summary_df,
        outdir / "coassociation_40um_bins_publication.png",
        outdir / "coassociation_40um_bins_publication.pdf",
    )


if __name__ == "__main__":
    main()
