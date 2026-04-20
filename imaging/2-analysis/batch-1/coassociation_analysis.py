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
    ("0_30", 0.0, 30.0),
    ("30_60", 30.0, 60.0),
    ("60_90", 60.0, 90.0),
    ("90_120", 90.0, 120.0),
]


@dataclass(frozen=True)
class BinSpec:
    name: str
    low: float | None
    high: float | None


def canonicalize_source_file(name: str) -> str:
    return re.sub(r"^C\d+-", "", name)


def load_puncta_table(path: str, puncta_type: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(
        columns={
            "Source File": "source_file",
            "X (um)": "x_um",
            "Y (um)": "y_um",
        }
    )
    df["puncta_type"] = puncta_type
    df["image_id"] = df["source_file"].map(canonicalize_source_file)
    return df[["image_id", "source_file", "x_um", "y_um", "puncta_type"]]


def subset_points(points: np.ndarray, bin_spec: BinSpec) -> tuple[np.ndarray, np.ndarray]:
    if bin_spec.low is None:
        indices = np.arange(len(points))
        return points, indices
    mask = (points[:, 0] >= bin_spec.low) & (points[:, 0] < bin_spec.high)
    indices = np.flatnonzero(mask)
    return points[mask], indices


def overlap_matrix(
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    threshold_um: float,
) -> np.ndarray:
    if len(actin_points) == 0 or len(arc_points) == 0:
        return np.zeros((len(actin_points), len(arc_points)), dtype=bool)
    deltas = actin_points[:, None, :] - arc_points[None, :, :]
    return np.sum(deltas * deltas, axis=2) <= (threshold_um ** 2)


def compute_overlap_metrics(
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    threshold_um: float,
) -> dict[str, float]:
    within_threshold = overlap_matrix(actin_points, arc_points, threshold_um)
    actin_n = int(actin_points.shape[0])
    arc_n = int(arc_points.shape[0])
    actin_assoc_n = int(within_threshold.any(axis=1).sum()) if actin_n else 0
    arc_assoc_n = int(within_threshold.any(axis=0).sum()) if arc_n else 0
    total_puncta = actin_n + arc_n
    total_assoc_puncta = actin_assoc_n + arc_assoc_n
    pair_count = int(within_threshold.sum())

    return {
        "actin_n": actin_n,
        "arc_n": arc_n,
        "total_puncta": total_puncta,
        "actin_assoc_n": actin_assoc_n,
        "arc_assoc_n": arc_assoc_n,
        "total_assoc_puncta": total_assoc_puncta,
        "pair_count": pair_count,
        "coassoc_freq": (total_assoc_puncta / total_puncta) if total_puncta else np.nan,
    }


def permute_points(points: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(points) <= 1:
        return points.copy()
    return np.column_stack((rng.permutation(points[:, 0]), rng.permutation(points[:, 1])))


def summarize_with_null(observed: np.ndarray, null_distribution: np.ndarray) -> dict[str, float]:
    valid_observed = observed[~np.isnan(observed)]
    observed_mean = float(np.mean(valid_observed)) if valid_observed.size else np.nan

    if null_distribution.size == 0:
        return {
            "observed_mean": observed_mean,
            "observed_sd": float(np.std(valid_observed, ddof=1)) if len(valid_observed) > 1 else np.nan,
            "observed_sem": (
                float(np.std(valid_observed, ddof=1) / np.sqrt(len(valid_observed)))
                if len(valid_observed) > 1
                else np.nan
            ),
            "null_mean": np.nan,
            "null_sd": np.nan,
            "delta_vs_null": np.nan,
            "z_score": np.nan,
            "p_value_enrichment": np.nan,
            "null_ci_low": np.nan,
            "null_ci_high": np.nan,
            "n_images": int(valid_observed.size),
        }

    null_mean = float(np.mean(null_distribution))
    null_sd = float(np.std(null_distribution, ddof=1)) if len(null_distribution) > 1 else 0.0
    observed_sd = float(np.std(valid_observed, ddof=1)) if len(valid_observed) > 1 else np.nan
    observed_sem = observed_sd / np.sqrt(len(valid_observed)) if len(valid_observed) > 1 else np.nan
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
        "n_images": int(valid_observed.size),
    }


def sanitize_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def plot_bin_summary(
    per_image_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_path: Path,
) -> None:
    order = [bin_name for bin_name, _, _ in BINS]
    label_map = {
        "whole_image": "Whole image",
        "0_30": "0-30 um",
        "30_60": "30-60 um",
        "60_90": "60-90 um",
        "90_120": "90-120 um",
    }
    plot_df = per_image_df.dropna(subset=["coassoc_freq"]).copy()
    plot_df["bin_label"] = plot_df["bin"].map(label_map)
    summary_df = summary_df.copy()
    summary_df["bin_label"] = summary_df["bin"].map(label_map)

    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(11, 6.5))
    sns.stripplot(
        data=plot_df,
        x="bin_label",
        y="coassoc_freq",
        order=[label_map[name] for name in order],
        color="#274c77",
        alpha=0.55,
        jitter=0.22,
        size=5,
        ax=ax,
    )

    x_positions = np.arange(len(order))
    ax.errorbar(
        x=x_positions,
        y=summary_df["observed_mean"],
        yerr=summary_df["observed_sem"],
        fmt="D",
        color="#c1121f",
        markersize=7,
        capsize=4,
        linewidth=1.5,
        label="Observed mean ± SEM",
        zorder=5,
    )

    ax.errorbar(
        x=x_positions,
        y=summary_df["null_mean"],
        yerr=[
            summary_df["null_mean"] - summary_df["null_ci_low"],
            summary_df["null_ci_high"] - summary_df["null_mean"],
        ],
        fmt="s",
        color="#6c757d",
        markersize=6,
        capsize=4,
        linewidth=1.5,
        label="Null mean ± 95% interval",
        zorder=4,
    )

    for x_pos, (_, row) in zip(x_positions, summary_df.iterrows()):
        p_text = f"p={row['p_value_enrichment']:.3g}" if pd.notna(row["p_value_enrichment"]) else "p=NA"
        n_text = f"n={int(row['n_images'])}"
        ax.text(
            x_pos,
            max(row["observed_mean"], row["null_ci_high"]) + 0.035,
            f"{p_text}\n{n_text}",
            ha="center",
            va="bottom",
            fontsize=9,
            color="#343a40",
        )

    ax.set_xlabel("Distance from soma")
    ax.set_ylabel("Co-associated puncta / total puncta")
    ax.set_title("Actin-Arc co-association by image and distance bin")
    ax.set_ylim(bottom=-0.01)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_image_map(
    image_id: str,
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    threshold_um: float,
    output_dir: Path,
) -> None:
    within_threshold = overlap_matrix(actin_points, arc_points, threshold_um)
    actin_assoc = within_threshold.any(axis=1) if len(actin_points) else np.array([], dtype=bool)
    arc_assoc = within_threshold.any(axis=0) if len(arc_points) else np.array([], dtype=bool)

    fig, ax = plt.subplots(figsize=(9, 3.4))
    if len(actin_points):
        ax.scatter(
            actin_points[~actin_assoc, 0],
            actin_points[~actin_assoc, 1],
            s=28,
            c="#7aa6c2",
            alpha=0.55,
            label="Actin",
        )
        ax.scatter(
            actin_points[actin_assoc, 0],
            actin_points[actin_assoc, 1],
            s=62,
            c="#0b5d7a",
            edgecolors="white",
            linewidths=0.7,
            label="Actin co-associated",
            zorder=3,
        )
    if len(arc_points):
        ax.scatter(
            arc_points[~arc_assoc, 0],
            arc_points[~arc_assoc, 1],
            s=28,
            c="#f4a261",
            alpha=0.55,
            marker="^",
            label="Arc",
        )
        ax.scatter(
            arc_points[arc_assoc, 0],
            arc_points[arc_assoc, 1],
            s=72,
            c="#c75100",
            edgecolors="white",
            linewidths=0.7,
            marker="^",
            label="Arc co-associated",
            zorder=3,
        )

    pair_indices = np.argwhere(within_threshold)
    for actin_idx, arc_idx in pair_indices:
        ax.plot(
            [actin_points[actin_idx, 0], arc_points[arc_idx, 0]],
            [actin_points[actin_idx, 1], arc_points[arc_idx, 1]],
            color="#adb5bd",
            alpha=0.35,
            linewidth=0.8,
            zorder=1,
        )

    ax.set_title(f"{image_id}\nCo-associated puncta within {threshold_um:.1f} um")
    ax.set_xlabel("Distance from soma (um)")
    ax.set_ylabel("Y position (um)")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_xlim(left=-1)
    ax.set_ylim(bottom=-0.2)
    fig.tight_layout()
    fig.savefig(output_dir / f"{sanitize_filename(image_id)}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actin", default="actin puncta locations.csv")
    parser.add_argument("--arc", default="arc puncta locations.csv")
    parser.add_argument("--threshold-um", type=float, default=0.1)
    parser.add_argument("--iterations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260313)
    parser.add_argument("--summary-out", default="coassociation_summary.csv")
    parser.add_argument("--per-image-out", default="coassociation_per_image.csv")
    parser.add_argument("--summary-plot-out", default="coassociation_by_bin.png")
    parser.add_argument("--image-plot-dir", default="coassociation_image_plots")
    args = parser.parse_args()

    actin_df = load_puncta_table(args.actin, "actin")
    arc_df = load_puncta_table(args.arc, "arc")

    shared_images = sorted(set(actin_df["image_id"]) & set(arc_df["image_id"]))
    if not shared_images:
        raise SystemExit("No shared images after canonicalizing source file names.")

    bin_specs = [BinSpec(*bin_spec) for bin_spec in BINS]
    rng = np.random.default_rng(args.seed)
    image_plot_dir = Path(args.image_plot_dir)
    image_plot_dir.mkdir(parents=True, exist_ok=True)

    actin_by_image = {
        image_id: actin_df.loc[actin_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }
    arc_by_image = {
        image_id: arc_df.loc[arc_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }

    per_image_rows: list[dict[str, float | str | int]] = []
    observed_metric_values: dict[str, list[float]] = {}
    null_metric_values: dict[str, list[list[float]]] = {}
    for bin_spec in bin_specs:
        observed_metric_values[bin_spec.name] = []
        null_metric_values[bin_spec.name] = []

    for image_id in shared_images:
        actin_points = actin_by_image[image_id]
        arc_points = arc_by_image[image_id]

        plot_image_map(
            image_id=image_id,
            actin_points=actin_points,
            arc_points=arc_points,
            threshold_um=args.threshold_um,
            output_dir=image_plot_dir,
        )

        image_null_rows = {bin_spec.name: [] for bin_spec in bin_specs}
        for bin_spec in bin_specs:
            actin_subset, _ = subset_points(actin_points, bin_spec)
            arc_subset, _ = subset_points(arc_points, bin_spec)
            observed = compute_overlap_metrics(actin_subset, arc_subset, args.threshold_um)
            per_image_rows.append(
                {
                    "image_id": image_id,
                    "bin": bin_spec.name,
                    **observed,
                }
            )
            observed_metric_values[bin_spec.name].append(observed["coassoc_freq"])

        for _ in range(args.iterations):
            permuted_actin = permute_points(actin_points, rng)
            permuted_arc = permute_points(arc_points, rng)
            for bin_spec in bin_specs:
                permuted_actin_subset, _ = subset_points(permuted_actin, bin_spec)
                permuted_arc_subset, _ = subset_points(permuted_arc, bin_spec)
                permuted = compute_overlap_metrics(permuted_actin_subset, permuted_arc_subset, args.threshold_um)
                image_null_rows[bin_spec.name].append(permuted["coassoc_freq"])

        for bin_name, values in image_null_rows.items():
            null_metric_values[bin_name].append(values)

    summary_rows: list[dict[str, float | str | int]] = []
    for bin_spec in bin_specs:
        observed = np.asarray(observed_metric_values[bin_spec.name], dtype=float)
        per_image_null = np.asarray(null_metric_values[bin_spec.name], dtype=float)

        valid_mask = ~np.isnan(observed)
        valid_observed = observed[valid_mask]
        valid_null = per_image_null[valid_mask, :]
        group_null = np.mean(valid_null, axis=0) if valid_null.size else np.asarray([])

        summary_rows.append(
            {
                "bin": bin_spec.name,
                **summarize_with_null(valid_observed, group_null),
            }
        )

    per_image_df = pd.DataFrame(per_image_rows)
    summary_df = pd.DataFrame(summary_rows)
    per_image_df.to_csv(args.per_image_out, index=False)
    summary_df.to_csv(args.summary_out, index=False)
    plot_bin_summary(per_image_df, summary_df, Path(args.summary_plot_out))


if __name__ == "__main__":
    main()
