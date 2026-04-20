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

NULL_MODELS = {
    "permute_xy": "Permute x/y reference",
    "uniform_random": "Uniform random placement",
    "resample_with_replacement": "Resample puncta with replacement",
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


def compute_coassociation(actin_points: np.ndarray, arc_points: np.ndarray, threshold_um: float) -> dict[str, float]:
    within_threshold = overlap_matrix(actin_points, arc_points, threshold_um)
    actin_n = int(len(actin_points))
    arc_n = int(len(arc_points))
    actin_assoc_n = int(within_threshold.any(axis=1).sum()) if actin_n else 0
    arc_assoc_n = int(within_threshold.any(axis=0).sum()) if arc_n else 0
    total_puncta = actin_n + arc_n
    total_assoc = actin_assoc_n + arc_assoc_n
    pair_count = int(within_threshold.sum())
    return {
        "actin_n": actin_n,
        "arc_n": arc_n,
        "total_puncta": total_puncta,
        "actin_assoc_n": actin_assoc_n,
        "arc_assoc_n": arc_assoc_n,
        "total_assoc_puncta": total_assoc,
        "pair_count": pair_count,
        "coassoc_freq": (total_assoc / total_puncta) if total_puncta else 0.0,
    }


def permute_xy(points: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(points) <= 1:
        return points.copy()
    return np.column_stack((rng.permutation(points[:, 0]), rng.permutation(points[:, 1])))


def uniform_random(points: np.ndarray, bounds: tuple[float, float, float, float], rng: np.random.Generator) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    x_min, x_max, y_min, y_max = bounds
    x = rng.uniform(x_min, x_max, size=len(points))
    y = rng.uniform(y_min, y_max, size=len(points))
    return np.column_stack((x, y))


def resample_with_replacement(points: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    if len(points) == 0:
        return points.copy()
    indices = rng.integers(0, len(points), size=len(points))
    return points[indices]


def generate_null_points(
    mode: str,
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    bounds: tuple[float, float, float, float],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    if mode == "permute_xy":
        return permute_xy(actin_points, rng), permute_xy(arc_points, rng)
    if mode == "uniform_random":
        return uniform_random(actin_points, bounds, rng), uniform_random(arc_points, bounds, rng)
    if mode == "resample_with_replacement":
        return resample_with_replacement(actin_points, rng), resample_with_replacement(arc_points, rng)
    raise ValueError(f"Unsupported null model: {mode}")


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


def plot_null_model_comparison(
    observed_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    output_path: Path,
) -> None:
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

    fig, axes = plt.subplots(1, len(NULL_MODELS), figsize=(18, 6.5), sharey=True)
    if len(NULL_MODELS) == 1:
        axes = [axes]

    observed_plot = observed_df.dropna(subset=["coassoc_freq"]).copy()
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
            alpha=0.5,
            jitter=0.22,
            size=4.5,
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
            label="Observed mean ± SEM",
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
            label="Null mean ± 95% interval",
        )

        for x_pos, (_, row) in zip(x_positions, model_summary.iterrows()):
            p_text = f"p={row['p_value_enrichment']:.3g}" if pd.notna(row["p_value_enrichment"]) else "p=NA"
            ax.text(
                x_pos,
                max(row["observed_mean"], row["null_ci_high"]) + 0.03,
                f"{p_text}\nn={int(row['n_images'])}",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color="#343a40",
            )

        ax.set_title(NULL_MODELS[model])
        ax.set_xlabel("Distance from soma")
        ax.set_ylim(bottom=-0.01)
        if ax is axes[0]:
            ax.set_ylabel("Co-associated puncta / total puncta")
        else:
            ax.set_ylabel("")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Observed co-association compared against alternative null models", y=1.08, fontsize=15)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actin", default="../actin puncta locations.csv")
    parser.add_argument("--arc", default="../arc puncta locations.csv")
    parser.add_argument("--threshold-um", type=float, default=0.1)
    parser.add_argument("--iterations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260316)
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    actin_df = load_puncta_table(args.actin)
    arc_df = load_puncta_table(args.arc)
    shared_images = sorted(set(actin_df["image_id"]) & set(arc_df["image_id"]))
    if not shared_images:
        raise SystemExit("No shared images after canonicalizing source file names.")

    bin_specs = [BinSpec(*bin_spec) for bin_spec in BINS]
    rng = np.random.default_rng(args.seed)

    actin_by_image = {
        image_id: actin_df.loc[actin_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }
    arc_by_image = {
        image_id: arc_df.loc[arc_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        for image_id in shared_images
    }
    bounds_by_image = {}
    for image_id in shared_images:
        combined = np.vstack((actin_by_image[image_id], arc_by_image[image_id]))
        x_max = float(np.max(combined[:, 0])) if len(combined) else 0.0
        y_max = float(np.max(combined[:, 1])) if len(combined) else 0.0
        bounds_by_image[image_id] = (0.0, x_max, 0.0, y_max)

    observed_rows: list[dict[str, float | str | int]] = []
    observed_metric_values: dict[str, list[float]] = {bin_spec.name: [] for bin_spec in bin_specs}
    null_metric_values: dict[tuple[str, str], list[list[float]]] = {
        (model, bin_spec.name): [] for model in NULL_MODELS for bin_spec in bin_specs
    }

    for image_id in shared_images:
        actin_points = actin_by_image[image_id]
        arc_points = arc_by_image[image_id]
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
        for _ in range(args.iterations):
            for model in NULL_MODELS:
                null_actin, null_arc = generate_null_points(model, actin_points, arc_points, bounds, rng)
                for bin_spec in bin_specs:
                    null_metrics = compute_coassociation(
                        subset_points(null_actin, bin_spec),
                        subset_points(null_arc, bin_spec),
                        args.threshold_um,
                    )
                    per_model_null[(model, bin_spec.name)].append(null_metrics["coassoc_freq"])

        for key, values in per_model_null.items():
            null_metric_values[key].append(values)

    observed_df = pd.DataFrame(observed_rows)
    observed_df.to_csv(outdir / "observed_per_image.csv", index=False)

    summary_rows: list[dict[str, float | str | int]] = []
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
            summary_rows.append(row)
            model_rows.append(row)
        pd.DataFrame(model_rows).to_csv(model_dir / "summary.csv", index=False)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "null_model_summary.csv", index=False)
    plot_null_model_comparison(observed_df, summary_df, outdir / "null_model_comparison.png")


if __name__ == "__main__":
    main()
