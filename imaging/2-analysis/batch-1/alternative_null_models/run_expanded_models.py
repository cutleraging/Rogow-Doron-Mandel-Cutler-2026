#!/usr/bin/env python3

import argparse
import math
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

NULL_MODELS = {
    "permute_xy": "Permute x/y reference",
    "uniform_random": "Uniform random placement",
    "resample_with_replacement": "Resample puncta with replacement",
    "local_jitter": "Local jitter window",
    "axial_shift_arc": "Axial shift of arc channel",
    "local_masked_randomization": "Local x-density / y-mask randomization",
}


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


def pairwise_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if len(a) == 0 or len(b) == 0:
        return np.zeros((len(a), len(b)), dtype=float)
    deltas = a[:, None, :] - b[None, :, :]
    return np.sqrt(np.sum(deltas * deltas, axis=2))


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


def summarize_observed(values: np.ndarray) -> dict[str, float]:
    mean = float(np.mean(values)) if len(values) else np.nan
    sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
    sem = sd / np.sqrt(len(values)) if len(values) > 1 else np.nan
    return {
        "mean": mean,
        "sd": sd,
        "sem": sem,
        "n_images": int(len(values)),
    }


def permute_xy(points: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(points) <= 1:
        return points.copy()
    return np.column_stack((rng.permutation(points[:, 0]), rng.permutation(points[:, 1])))


def uniform_random(points: np.ndarray, bounds: tuple[float, float, float, float], rng: np.random.Generator) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    x_min, x_max, y_min, y_max = bounds
    return np.column_stack(
        (
            rng.uniform(x_min, x_max, size=len(points)),
            rng.uniform(y_min, y_max, size=len(points)),
        )
    )


def resample_with_replacement(points: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    indices = rng.integers(0, len(points), size=len(points))
    return points[indices]


def local_jitter(
    points: np.ndarray,
    bounds: tuple[float, float, float, float],
    rng: np.random.Generator,
    x_half_window: float,
    y_half_window: float,
) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    x_min, x_max, y_min, y_max = bounds
    jittered = points.copy()
    jittered[:, 0] += rng.uniform(-x_half_window, x_half_window, size=len(points))
    jittered[:, 1] += rng.uniform(-y_half_window, y_half_window, size=len(points))
    jittered[:, 0] = np.clip(jittered[:, 0], x_min, x_max)
    jittered[:, 1] = np.clip(jittered[:, 1], y_min, y_max)
    return jittered


def axial_shift_arc(points: np.ndarray, x_extent: float, rng: np.random.Generator, min_shift_um: float) -> np.ndarray:
    if len(points) == 0 or x_extent <= 0:
        return points.copy()
    max_shift = max(x_extent - min_shift_um, min_shift_um)
    shift = rng.uniform(min_shift_um, max_shift) if max_shift > min_shift_um else min_shift_um
    shifted = points.copy()
    shifted[:, 0] = np.mod(shifted[:, 0] + shift, x_extent)
    return shifted


def local_masked_randomization(
    points: np.ndarray,
    reference_points: np.ndarray,
    bounds: tuple[float, float, float, float],
    rng: np.random.Generator,
    x_window_um: float,
) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    x_min, x_max, y_min_global, y_max_global = bounds
    ref_x = reference_points[:, 0]
    ref_y = reference_points[:, 1]
    randomized = np.empty_like(points)
    for idx, (x0, _) in enumerate(points):
        local_mask = np.abs(ref_x - x0) <= x_window_um
        if np.any(local_mask):
            local_x = ref_x[local_mask]
            local_y = ref_y[local_mask]
            new_x = float(rng.choice(local_x))
            y_low = float(np.min(local_y))
            y_high = float(np.max(local_y))
        else:
            new_x = float(rng.uniform(x_min, x_max))
            y_low = y_min_global
            y_high = y_max_global
        if y_high <= y_low:
            new_y = y_low
        else:
            new_y = float(rng.uniform(y_low, y_high))
        randomized[idx, 0] = np.clip(new_x, x_min, x_max)
        randomized[idx, 1] = np.clip(new_y, y_min_global, y_max_global)
    return randomized


def generate_null_points(
    mode: str,
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    combined_points: np.ndarray,
    bounds: tuple[float, float, float, float],
    rng: np.random.Generator,
    jitter_x: float,
    jitter_y: float,
    local_x_window: float,
    min_axial_shift: float,
) -> tuple[np.ndarray, np.ndarray]:
    x_extent = bounds[1] - bounds[0]
    if mode == "permute_xy":
        return permute_xy(actin_points, rng), permute_xy(arc_points, rng)
    if mode == "uniform_random":
        return uniform_random(actin_points, bounds, rng), uniform_random(arc_points, bounds, rng)
    if mode == "resample_with_replacement":
        return resample_with_replacement(actin_points, rng), resample_with_replacement(arc_points, rng)
    if mode == "local_jitter":
        return (
            local_jitter(actin_points, bounds, rng, jitter_x, jitter_y),
            local_jitter(arc_points, bounds, rng, jitter_x, jitter_y),
        )
    if mode == "axial_shift_arc":
        return actin_points.copy(), axial_shift_arc(arc_points, x_extent, rng, min_axial_shift)
    if mode == "local_masked_randomization":
        return (
            local_masked_randomization(actin_points, combined_points, bounds, rng, local_x_window),
            local_masked_randomization(arc_points, combined_points, bounds, rng, local_x_window),
        )
    raise ValueError(f"Unsupported null model: {mode}")


def mutual_nearest_neighbor_pairs(actin_points: np.ndarray, arc_points: np.ndarray) -> list[tuple[int, int, float]]:
    if len(actin_points) == 0 or len(arc_points) == 0:
        return []
    distances = pairwise_distances(actin_points, arc_points)
    nearest_arc = np.argmin(distances, axis=1)
    nearest_actin = np.argmin(distances, axis=0)
    pairs: list[tuple[int, int, float]] = []
    for actin_idx, arc_idx in enumerate(nearest_arc):
        if nearest_actin[arc_idx] == actin_idx:
            pairs.append((actin_idx, int(arc_idx), float(distances[actin_idx, arc_idx])))
    return pairs


def local_monte_carlo_pair_pvalue(
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    pair_midpoint: np.ndarray,
    observed_distance: float,
    image_bounds: tuple[float, float, float, float],
    rng: np.random.Generator,
    iterations: int,
    x_half_window: float,
    y_half_window: float,
) -> float:
    x_min = max(image_bounds[0], pair_midpoint[0] - x_half_window)
    x_max = min(image_bounds[1], pair_midpoint[0] + x_half_window)
    y_min = max(image_bounds[2], pair_midpoint[1] - y_half_window)
    y_max = min(image_bounds[3], pair_midpoint[1] + y_half_window)
    if x_max <= x_min:
        x_min, x_max = image_bounds[0], image_bounds[1]
    if y_max <= y_min:
        y_min, y_max = image_bounds[2], image_bounds[3]

    actin_mask = (
        (actin_points[:, 0] >= x_min)
        & (actin_points[:, 0] <= x_max)
        & (actin_points[:, 1] >= y_min)
        & (actin_points[:, 1] <= y_max)
    )
    arc_mask = (
        (arc_points[:, 0] >= x_min)
        & (arc_points[:, 0] <= x_max)
        & (arc_points[:, 1] >= y_min)
        & (arc_points[:, 1] <= y_max)
    )
    local_actin_n = int(np.sum(actin_mask))
    local_arc_n = int(np.sum(arc_mask))
    if local_actin_n == 0 or local_arc_n == 0:
        return 1.0

    min_distances = np.empty(iterations, dtype=float)
    for i in range(iterations):
        sim_actin = np.column_stack(
            (
                rng.uniform(x_min, x_max, size=local_actin_n),
                rng.uniform(y_min, y_max, size=local_actin_n),
            )
        )
        sim_arc = np.column_stack(
            (
                rng.uniform(x_min, x_max, size=local_arc_n),
                rng.uniform(y_min, y_max, size=local_arc_n),
            )
        )
        sim_distances = pairwise_distances(sim_actin, sim_arc)
        min_distances[i] = float(np.min(sim_distances))

    return float((1 + np.sum(min_distances <= observed_distance)) / (iterations + 1))


def bauer_inspired_metrics(
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    threshold_um: float,
    image_bounds: tuple[float, float, float, float],
    rng: np.random.Generator,
    iterations: int,
    x_half_window: float,
    y_half_window: float,
) -> tuple[dict[str, float], list[dict[str, float]]]:
    pairs = mutual_nearest_neighbor_pairs(actin_points, arc_points)
    total_puncta = len(actin_points) + len(arc_points)
    bin_rows = []
    for bin_name, _, _ in BINS:
        bin_rows.append(
            {
                "bin": bin_name,
                "specific_pair_count": 0,
                "specific_puncta_count": 0,
                "pair_count": 0,
                "specific_freq": 0.0,
                "pair_midpoint_x_mean": np.nan,
            }
        )
    if total_puncta == 0 or not pairs:
        return {"whole_image_specific_freq": 0.0, "whole_image_specific_pairs": 0}, bin_rows

    bin_lookup = {name: idx for idx, (name, _, _) in enumerate(BINS)}
    pair_midpoints_by_bin: dict[str, list[float]] = {name: [] for name, _, _ in BINS}
    specific_pair_total = 0
    pair_total = 0
    for actin_idx, arc_idx, distance in pairs:
        pair_total += 1
        actin_point = actin_points[actin_idx]
        arc_point = arc_points[arc_idx]
        midpoint = (actin_point + arc_point) / 2.0
        pair_bin = "whole_image"
        for name, low, high in BINS[1:]:
            if low <= midpoint[0] < high:
                pair_bin = name
                break
        p_local = local_monte_carlo_pair_pvalue(
            actin_points,
            arc_points,
            midpoint,
            distance,
            image_bounds,
            rng,
            iterations,
            x_half_window,
            y_half_window,
        )
        if distance <= threshold_um and p_local < 0.05:
            specific_pair_total += 1
            for target_bin in ("whole_image", pair_bin):
                row = bin_rows[bin_lookup[target_bin]]
                row["specific_pair_count"] += 1
                row["specific_puncta_count"] += 2
                pair_midpoints_by_bin[target_bin].append(float(midpoint[0]))
        for target_bin in ("whole_image", pair_bin):
            bin_rows[bin_lookup[target_bin]]["pair_count"] += 1

    for row in bin_rows:
        row["specific_freq"] = row["specific_puncta_count"] / total_puncta if total_puncta else 0.0
        bin_midpoints = pair_midpoints_by_bin[row["bin"]]
        row["pair_midpoint_x_mean"] = float(np.mean(bin_midpoints)) if bin_midpoints else np.nan

    return {
        "whole_image_specific_freq": (2 * specific_pair_total / total_puncta) if total_puncta else 0.0,
        "whole_image_specific_pairs": specific_pair_total,
        "whole_image_mnn_pairs": pair_total,
    }, bin_rows


def rankdata(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average").to_numpy(dtype=float)


def spearman_from_vectors(a: np.ndarray, b: np.ndarray) -> float:
    if np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return 0.0
    a_rank = rankdata(a)
    b_rank = rankdata(b)
    corr = np.corrcoef(a_rank, b_rank)[0, 1]
    return float(corr) if np.isfinite(corr) else 0.0


def cbc_scores_for_channel(
    source_points: np.ndarray,
    same_points: np.ndarray,
    other_points: np.ndarray,
    radii: np.ndarray,
) -> np.ndarray:
    if len(source_points) == 0:
        return np.zeros(0, dtype=float)
    same_dist = pairwise_distances(source_points, same_points)
    other_dist = pairwise_distances(source_points, other_points)
    scores = np.zeros(len(source_points), dtype=float)
    r_max = float(np.max(radii))

    for idx in range(len(source_points)):
        same_counts = []
        other_counts = []
        for radius in radii:
            if len(same_points):
                counts_same = np.sum(same_dist[idx] <= radius)
                if same_points is source_points:
                    counts_same -= 1
            else:
                counts_same = 0
            counts_other = int(np.sum(other_dist[idx] <= radius)) if len(other_points) else 0
            same_counts.append(counts_same)
            other_counts.append(counts_other)
        rho = spearman_from_vectors(np.asarray(same_counts, dtype=float), np.asarray(other_counts, dtype=float))
        nearest_other = float(np.min(other_dist[idx])) if len(other_points) else np.inf
        penalty = math.exp(-nearest_other / r_max) if np.isfinite(nearest_other) else 0.0
        scores[idx] = rho * penalty
    return scores


def cbc_metrics(
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    actin_scores = cbc_scores_for_channel(actin_points, actin_points, arc_points, radii)
    arc_scores = cbc_scores_for_channel(arc_points, arc_points, actin_points, radii)
    return actin_scores, arc_scores


def plot_null_model_comparison(observed_df: pd.DataFrame, summary_df: pd.DataFrame, output_path: Path) -> None:
    sns.set_theme(style="whitegrid")
    label_map = {
        "whole_image": "Whole image",
        "0_30": "0-30 um",
        "30_60": "30-60 um",
        "60_90": "60-90 um",
        "90_120": "90-120 um",
    }
    order = [name for name, _, _ in BINS]
    x_labels = [label_map[name] for name in order]
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=True)
    axes = axes.ravel()
    observed_plot = observed_df.copy()
    observed_plot["bin_label"] = observed_plot["bin"].map(label_map)

    for ax, model in zip(axes, NULL_MODELS):
        model_summary = summary_df[summary_df["null_model"] == model].copy()
        model_summary["bin_label"] = model_summary["bin"].map(label_map)
        sns.stripplot(
            data=observed_plot,
            x="bin_label",
            y="coassoc_freq",
            order=x_labels,
            color="#355070",
            alpha=0.45,
            jitter=0.22,
            size=4,
            ax=ax,
        )
        x_positions = np.arange(len(order))
        ax.errorbar(
            x_positions,
            model_summary["observed_mean"],
            yerr=model_summary["observed_sem"],
            fmt="D",
            color="#b80c09",
            markersize=6,
            capsize=4,
            linewidth=1.4,
            zorder=5,
        )
        ax.errorbar(
            x_positions,
            model_summary["null_mean"],
            yerr=[
                model_summary["null_mean"] - model_summary["null_ci_low"],
                model_summary["null_ci_high"] - model_summary["null_mean"],
            ],
            fmt="s",
            color="#6c757d",
            markersize=5,
            capsize=4,
            linewidth=1.4,
            zorder=4,
        )
        for x_pos, (_, row) in zip(x_positions, model_summary.iterrows()):
            ax.text(
                x_pos,
                max(row["observed_mean"], row["null_ci_high"]) + 0.03,
                f"p={row['p_value_enrichment']:.3g}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#343a40",
            )
        ax.set_title(NULL_MODELS[model])
        ax.set_xlabel("Distance from soma")
        ax.set_ylim(bottom=-0.01)
        if ax in (axes[0], axes[3]):
            ax.set_ylabel("Co-associated puncta / total puncta")
        else:
            ax.set_ylabel("")
    fig.suptitle("Observed co-association against expanded null models", y=1.02, fontsize=15)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_observed_summary(df: pd.DataFrame, value_column: str, title: str, ylabel: str, output_path: Path) -> None:
    sns.set_theme(style="whitegrid")
    label_map = {
        "whole_image": "Whole image",
        "0_30": "0-30 um",
        "30_60": "30-60 um",
        "60_90": "60-90 um",
        "90_120": "90-120 um",
    }
    order = [name for name, _, _ in BINS]
    plot_df = df.copy()
    plot_df["bin_label"] = plot_df["bin"].map(label_map)
    summary_rows = []
    for bin_name in order:
        values = plot_df.loc[plot_df["bin"] == bin_name, value_column].to_numpy(dtype=float)
        stats = summarize_observed(values)
        summary_rows.append({"bin": bin_name, **stats})
    summary_df = pd.DataFrame(summary_rows)

    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    sns.stripplot(
        data=plot_df,
        x="bin_label",
        y=value_column,
        order=[label_map[name] for name in order],
        color="#355070",
        alpha=0.5,
        jitter=0.22,
        size=4.5,
        ax=ax,
    )
    x_positions = np.arange(len(order))
    ax.errorbar(
        x_positions,
        summary_df["mean"],
        yerr=summary_df["sem"],
        fmt="D",
        color="#b80c09",
        markersize=6,
        capsize=4,
        linewidth=1.4,
        zorder=5,
    )
    ax.axhline(0, color="#adb5bd", linewidth=1, linestyle="--")
    ax.set_title(title)
    ax.set_xlabel("Distance from soma")
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actin", default="../actin puncta locations.csv")
    parser.add_argument("--arc", default="../arc puncta locations.csv")
    parser.add_argument("--threshold-um", type=float, default=0.1)
    parser.add_argument("--null-iterations", type=int, default=3000)
    parser.add_argument("--bauer-iterations", type=int, default=250)
    parser.add_argument("--seed", type=int, default=20260316)
    parser.add_argument("--jitter-x-um", type=float, default=1.0)
    parser.add_argument("--jitter-y-um", type=float, default=0.5)
    parser.add_argument("--local-x-window-um", type=float, default=5.0)
    parser.add_argument("--min-axial-shift-um", type=float, default=5.0)
    parser.add_argument("--bauer-local-x-half-window-um", type=float, default=5.0)
    parser.add_argument("--bauer-local-y-half-window-um", type=float, default=1.0)
    parser.add_argument("--cbc-max-radius-um", type=float, default=2.0)
    parser.add_argument("--cbc-radius-steps", type=int, default=10)
    parser.add_argument("--outdir", default="expanded_results")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    bin_specs = [BinSpec(*bin_spec) for bin_spec in BINS]

    actin_df = load_puncta_table(args.actin)
    arc_df = load_puncta_table(args.arc)
    shared_images = sorted(set(actin_df["image_id"]) & set(arc_df["image_id"]))
    if not shared_images:
        raise SystemExit("No shared images after canonicalizing source file names.")

    actin_by_image = {
        image_id: actin_df.loc[actin_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }
    arc_by_image = {
        image_id: arc_df.loc[arc_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }
    bounds_by_image = {}
    combined_by_image = {}
    for image_id in shared_images:
        combined = np.vstack((actin_by_image[image_id], arc_by_image[image_id]))
        combined_by_image[image_id] = combined
        x_max = float(np.max(combined[:, 0])) if len(combined) else 0.0
        y_max = float(np.max(combined[:, 1])) if len(combined) else 0.0
        bounds_by_image[image_id] = (0.0, x_max, 0.0, y_max)

    observed_rows: list[dict[str, float | str | int]] = []
    observed_metric_values = {bin_spec.name: [] for bin_spec in bin_specs}
    null_metric_values = {(model, bin_spec.name): [] for model in NULL_MODELS for bin_spec in bin_specs}

    bauer_rows: list[dict[str, float | str | int]] = []
    cbc_rows: list[dict[str, float | str | int]] = []
    cbc_point_rows: list[dict[str, float | str | int]] = []
    cbc_radii = np.linspace(args.cbc_max_radius_um / args.cbc_radius_steps, args.cbc_max_radius_um, args.cbc_radius_steps)

    for image_id in shared_images:
        actin_points = actin_by_image[image_id]
        arc_points = arc_by_image[image_id]
        combined_points = combined_by_image[image_id]
        bounds = bounds_by_image[image_id]

        for bin_spec in bin_specs:
            observed = compute_coassociation(
                subset_points(actin_points, bin_spec),
                subset_points(arc_points, bin_spec),
                args.threshold_um,
            )
            observed_rows.append({"image_id": image_id, "bin": bin_spec.name, **observed})
            observed_metric_values[bin_spec.name].append(observed["coassoc_freq"])

        per_model_null = {(model, bin_spec.name): [] for model in NULL_MODELS for bin_spec in bin_specs}
        for _ in range(args.null_iterations):
            for model in NULL_MODELS:
                null_actin, null_arc = generate_null_points(
                    mode=model,
                    actin_points=actin_points,
                    arc_points=arc_points,
                    combined_points=combined_points,
                    bounds=bounds,
                    rng=rng,
                    jitter_x=args.jitter_x_um,
                    jitter_y=args.jitter_y_um,
                    local_x_window=args.local_x_window_um,
                    min_axial_shift=args.min_axial_shift_um,
                )
                for bin_spec in bin_specs:
                    null_metrics = compute_coassociation(
                        subset_points(null_actin, bin_spec),
                        subset_points(null_arc, bin_spec),
                        args.threshold_um,
                    )
                    per_model_null[(model, bin_spec.name)].append(null_metrics["coassoc_freq"])
        for key, values in per_model_null.items():
            null_metric_values[key].append(values)

        _, bauer_bin_rows = bauer_inspired_metrics(
            actin_points=actin_points,
            arc_points=arc_points,
            threshold_um=args.threshold_um,
            image_bounds=bounds,
            rng=rng,
            iterations=args.bauer_iterations,
            x_half_window=args.bauer_local_x_half_window_um,
            y_half_window=args.bauer_local_y_half_window_um,
        )
        total_puncta_by_bin = {}
        for bin_spec in bin_specs:
            total_puncta_by_bin[bin_spec.name] = len(subset_points(actin_points, bin_spec)) + len(subset_points(arc_points, bin_spec))
        for row in bauer_bin_rows:
            denominator = total_puncta_by_bin[row["bin"]]
            row["bauer_specific_freq_bin"] = row["specific_puncta_count"] / denominator if denominator else 0.0
            bauer_rows.append({"image_id": image_id, **row})

        actin_cbc, arc_cbc = cbc_metrics(actin_points, arc_points, cbc_radii)
        for channel_name, points, scores in (
            ("actin", actin_points, actin_cbc),
            ("arc", arc_points, arc_cbc),
        ):
            for point, score in zip(points, scores):
                cbc_point_rows.append(
                    {
                        "image_id": image_id,
                        "channel": channel_name,
                        "x_um": float(point[0]),
                        "y_um": float(point[1]),
                        "cbc_score": float(score),
                    }
                )
        combined_scores = np.concatenate((actin_cbc, arc_cbc))
        combined_x = np.concatenate((actin_points[:, 0], arc_points[:, 0]))
        for bin_spec in bin_specs:
            if bin_spec.low is None:
                bin_scores = combined_scores
            else:
                mask = (combined_x >= bin_spec.low) & (combined_x < bin_spec.high)
                bin_scores = combined_scores[mask]
            cbc_rows.append(
                {
                    "image_id": image_id,
                    "bin": bin_spec.name,
                    "cbc_mean_score": float(np.mean(bin_scores)) if len(bin_scores) else 0.0,
                    "cbc_positive_fraction": float(np.mean(bin_scores > 0.5)) if len(bin_scores) else 0.0,
                    "n_points": int(len(bin_scores)),
                }
            )

    observed_df = pd.DataFrame(observed_rows)
    observed_df.to_csv(outdir / "observed_per_image.csv", index=False)

    null_summary_rows = []
    for model in NULL_MODELS:
        model_dir = outdir / model
        model_dir.mkdir(exist_ok=True)
        model_rows = []
        for bin_spec in bin_specs:
            observed = np.asarray(observed_metric_values[bin_spec.name], dtype=float)
            per_image_null = np.asarray(null_metric_values[(model, bin_spec.name)], dtype=float)
            group_null = np.mean(per_image_null, axis=0) if per_image_null.size else np.asarray([])
            row = {
                "null_model": model,
                "null_model_label": NULL_MODELS[model],
                "bin": bin_spec.name,
                **summarize_with_null(observed, group_null),
            }
            null_summary_rows.append(row)
            model_rows.append(row)
        pd.DataFrame(model_rows).to_csv(model_dir / "summary.csv", index=False)
    null_summary_df = pd.DataFrame(null_summary_rows)
    null_summary_df.to_csv(outdir / "expanded_null_model_summary.csv", index=False)
    plot_null_model_comparison(observed_df, null_summary_df, outdir / "expanded_null_model_comparison.png")

    bauer_df = pd.DataFrame(bauer_rows)
    bauer_df.to_csv(outdir / "bauer_inspired_per_image.csv", index=False)
    bauer_summary_rows = []
    for bin_spec in bin_specs:
        values = bauer_df.loc[bauer_df["bin"] == bin_spec.name, "bauer_specific_freq_bin"].to_numpy(dtype=float)
        bauer_summary_rows.append({"bin": bin_spec.name, **summarize_observed(values)})
    bauer_summary_df = pd.DataFrame(bauer_summary_rows)
    bauer_summary_df.to_csv(outdir / "bauer_inspired_summary.csv", index=False)
    plot_observed_summary(
        bauer_df,
        value_column="bauer_specific_freq_bin",
        title="Bauer 2021-inspired local density corrected co-association",
        ylabel="Specific puncta / total puncta",
        output_path=outdir / "bauer_inspired_plot.png",
    )

    cbc_df = pd.DataFrame(cbc_rows)
    cbc_df.to_csv(outdir / "cbc_per_image.csv", index=False)
    pd.DataFrame(cbc_point_rows).to_csv(outdir / "cbc_per_punctum.csv", index=False)
    cbc_summary_rows = []
    for bin_spec in bin_specs:
        values = cbc_df.loc[cbc_df["bin"] == bin_spec.name, "cbc_mean_score"].to_numpy(dtype=float)
        cbc_summary_rows.append({"bin": bin_spec.name, **summarize_observed(values)})
    cbc_summary_df = pd.DataFrame(cbc_summary_rows)
    cbc_summary_df.to_csv(outdir / "cbc_summary.csv", index=False)
    plot_observed_summary(
        cbc_df,
        value_column="cbc_mean_score",
        title="Malkusch 2012-inspired coordinate-based colocalization score",
        ylabel="Mean CBC score",
        output_path=outdir / "cbc_plot.png",
    )


if __name__ == "__main__":
    main()
