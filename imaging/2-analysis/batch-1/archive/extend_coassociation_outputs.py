from __future__ import annotations

import argparse
import re
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SOURCE_COLUMN = "Source File"
COORD_COLUMNS = ["X (um)", "Y (um)", "Z (um)"]
BIN_ORDER = ["0_30", "30_60", "60_90", "90_120"]
REGION_LABELS = {"0_30": "0-30", "30_60": "30-60", "60_90": "60-90", "90_120": "90-120"}
GROUP_SPECS = [
    (
        "arc",
        "Arc puncta",
        "arc_coassociation_fraction",
        "null_arc_coassociation_fraction_mean",
        "null_arc_coassociation_fraction_ci_low",
        "null_arc_coassociation_fraction_ci_high",
        "coassociated_arc_observed",
        "arc_total",
    ),
    (
        "actin",
        "Actin puncta",
        "actin_coassociation_fraction",
        "null_actin_coassociation_fraction_mean",
        "null_actin_coassociation_fraction_ci_low",
        "null_actin_coassociation_fraction_ci_high",
        "coassociated_actin_observed",
        "actin_total",
    ),
    (
        "all_puncta",
        "All puncta",
        "all_puncta_coassociation_fraction",
        "null_all_puncta_coassociation_fraction_mean",
        "null_all_puncta_coassociation_fraction_ci_low",
        "null_all_puncta_coassociation_fraction_ci_high",
        "coassociated_all_puncta_observed",
        "all_puncta_total",
    ),
]


def normalize_source_name(value: str) -> str:
    return re.sub(r"^C\d-", "", value)


def hopcroft_karp_pairs(adjacency: list[list[int]], n_right: int) -> list[tuple[int, int]]:
    n_left = len(adjacency)
    pair_left = [-1] * n_left
    pair_right = [-1] * n_right
    dist = [0] * n_left
    infinity = n_left + n_right + 1

    def bfs() -> bool:
        queue: deque[int] = deque()
        found_augmenting_path = False
        for left_index in range(n_left):
            if pair_left[left_index] == -1:
                dist[left_index] = 0
                queue.append(left_index)
            else:
                dist[left_index] = infinity

        while queue:
            left_index = queue.popleft()
            for right_index in adjacency[left_index]:
                partner = pair_right[right_index]
                if partner == -1:
                    found_augmenting_path = True
                elif dist[partner] == infinity:
                    dist[partner] = dist[left_index] + 1
                    queue.append(partner)
        return found_augmenting_path

    def dfs(left_index: int) -> bool:
        for right_index in adjacency[left_index]:
            partner = pair_right[right_index]
            if partner == -1 or (dist[partner] == dist[left_index] + 1 and dfs(partner)):
                pair_left[left_index] = right_index
                pair_right[right_index] = left_index
                return True
        dist[left_index] = infinity
        return False

    while bfs():
        for left_index in range(n_left):
            if pair_left[left_index] == -1:
                dfs(left_index)

    return [(left_index, right_index) for left_index, right_index in enumerate(pair_left) if right_index != -1]


def match_pairs(actin_coords: np.ndarray, arc_coords: np.ndarray, threshold_um: float) -> list[tuple[int, int]]:
    if len(actin_coords) == 0 or len(arc_coords) == 0:
        return []
    delta = actin_coords[:, None, :] - arc_coords[None, :, :]
    distance_sq = np.sum(delta * delta, axis=2)
    threshold_sq = threshold_um * threshold_um
    adjacency: list[list[int]] = []
    for actin_index in range(distance_sq.shape[0]):
        neighbors = np.flatnonzero(distance_sq[actin_index] <= threshold_sq)
        adjacency.append(neighbors.tolist())
    return hopcroft_karp_pairs(adjacency, len(arc_coords))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate extended co-association plots and coordinate tables for actin/arc/all puncta."
    )
    parser.add_argument("--summary-csv", type=Path, required=True, help="Path to puncta_colocalization_summary.csv")
    parser.add_argument("--actin-csv", type=Path, required=True, help="Path to actin puncta CSV")
    parser.add_argument("--arc-csv", type=Path, required=True, help="Path to arc puncta CSV")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for plots and tables")
    parser.add_argument("--threshold-um", type=float, default=0.1, help="Co-association threshold in um")
    parser.add_argument("--null-model", type=str, default="fix_x_permute_yz", help="Null model label to filter")
    return parser.parse_args()


def load_summary(summary_path: Path, threshold_um: float, null_model: str) -> pd.DataFrame:
    summary = pd.read_csv(summary_path)
    summary = summary.loc[np.isclose(summary["threshold_um"], threshold_um)].copy()
    if "null_model" in summary.columns:
        summary = summary.loc[summary["null_model"] == null_model].copy()
    if summary.empty:
        raise ValueError(
            f"No summary rows found for threshold={threshold_um} and null_model={null_model} in {summary_path}."
        )
    return summary


def plot_whole_image_relative(summary: pd.DataFrame, output_path: Path) -> None:
    whole = summary.loc[summary["region"] == "whole_image"].iloc[0]
    x = np.arange(len(GROUP_SPECS))
    observed = np.array([100.0 * whole[spec[2]] for spec in GROUP_SPECS], dtype=float)
    null_mean = np.array([100.0 * whole[spec[3]] for spec in GROUP_SPECS], dtype=float)
    null_low = np.array([100.0 * whole[spec[4]] for spec in GROUP_SPECS], dtype=float)
    null_high = np.array([100.0 * whole[spec[5]] for spec in GROUP_SPECS], dtype=float)

    fig, ax = plt.subplots(figsize=(9, 5.2), dpi=220)
    width = 0.36
    ax.bar(
        x - width / 2,
        observed,
        width=width,
        color="#0f766e",
        edgecolor="#222222",
        linewidth=1.0,
        label="Observed",
    )
    ax.bar(
        x + width / 2,
        null_mean,
        width=width,
        color="#cbd5e1",
        edgecolor="#222222",
        linewidth=1.0,
        label="Null mean",
    )
    yerr_low = np.maximum(null_mean - null_low, 0.0)
    yerr_high = np.maximum(null_high - null_mean, 0.0)
    ax.errorbar(
        x + width / 2,
        null_mean,
        yerr=np.vstack((yerr_low, yerr_high)),
        fmt="none",
        ecolor="#111111",
        elinewidth=1.2,
        capsize=3,
        capthick=1.2,
        zorder=3,
    )

    ax.set_xticks(x, [spec[1] for spec in GROUP_SPECS])
    ax.set_ylabel("Co-associated fraction (%)")
    ax.set_title("Whole Image Co-association (100 nm)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)
    ax.set_ylim(0, max(np.max(observed), np.max(null_high)) * 1.25 + 1.0)

    ax.text(
        0.02,
        0.98,
        f"Enrichment vs null (pair-level): {whole['enrichment_vs_null']:.2f}x\nEmpirical p: {whole['empirical_p_enrichment']:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#dddddd"},
    )

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def plot_binned_relative(summary: pd.DataFrame, output_path: Path) -> None:
    bins = summary.loc[summary["region"].isin(BIN_ORDER)].copy()
    bins["region"] = pd.Categorical(bins["region"], categories=BIN_ORDER, ordered=True)
    bins = bins.sort_values("region")

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), dpi=220, sharey=True)
    x = np.arange(len(BIN_ORDER))
    width = 0.36

    for ax, spec in zip(axes, GROUP_SPECS):
        label = spec[1]
        observed = 100.0 * bins[spec[2]].to_numpy(dtype=float)
        null_mean = 100.0 * bins[spec[3]].to_numpy(dtype=float)
        null_low = 100.0 * bins[spec[4]].to_numpy(dtype=float)
        null_high = 100.0 * bins[spec[5]].to_numpy(dtype=float)
        coassociated_counts = bins[spec[6]].to_numpy(dtype=int)
        totals = bins[spec[7]].to_numpy(dtype=int)

        ax.bar(x - width / 2, observed, width=width, color="#0f766e", edgecolor="#222222", linewidth=1.0, label="Observed")
        ax.bar(x + width / 2, null_mean, width=width, color="#cbd5e1", edgecolor="#222222", linewidth=1.0, label="Null mean")
        yerr_low = np.maximum(null_mean - null_low, 0.0)
        yerr_high = np.maximum(null_high - null_mean, 0.0)
        ax.errorbar(
            x + width / 2,
            null_mean,
            yerr=np.vstack((yerr_low, yerr_high)),
            fmt="none",
            ecolor="#111111",
            elinewidth=1.1,
            capsize=3,
            capthick=1.1,
            zorder=3,
        )

        for idx in range(len(x)):
            ax.text(
                x[idx] - width / 2,
                observed[idx] + 0.6,
                f"{coassociated_counts[idx]}/{totals[idx]}",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=7.5,
            )

        ax.set_title(label)
        ax.set_xticks(x, [REGION_LABELS[region] for region in BIN_ORDER])
        ax.set_xlabel("Distance from soma (um)")
        ax.grid(axis="y", alpha=0.25)

    axes[0].set_ylabel("Co-associated fraction (%)")
    axes[0].legend(frameon=False, loc="upper right")
    axes[0].set_ylim(
        0,
        max(
            np.max(100.0 * bins["arc_coassociation_fraction"]),
            np.max(100.0 * bins["actin_coassociation_fraction"]),
            np.max(100.0 * bins["all_puncta_coassociation_fraction"]),
        )
        * 1.30
        + 1.0,
    )

    fig.suptitle("Distance-binned Co-association at 100 nm", y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def export_coassociated_coordinate_tables(
    actin_df: pd.DataFrame,
    arc_df: pd.DataFrame,
    threshold_um: float,
    pairs_output_path: Path,
    puncta_output_path: Path,
) -> tuple[int, int]:
    actin_df = actin_df.copy()
    arc_df = arc_df.copy()
    actin_df["image_key"] = actin_df[SOURCE_COLUMN].map(normalize_source_name)
    arc_df["image_key"] = arc_df[SOURCE_COLUMN].map(normalize_source_name)

    pair_rows: list[dict[str, float | int | str]] = []
    puncta_rows: list[dict[str, float | int | str]] = []
    pair_id = 0
    image_keys = sorted(set(actin_df["image_key"]) & set(arc_df["image_key"]))
    for image_key in image_keys:
        actin_rows = actin_df.loc[actin_df["image_key"] == image_key].reset_index(drop=True)
        arc_rows = arc_df.loc[arc_df["image_key"] == image_key].reset_index(drop=True)
        actin_xyz = actin_rows[COORD_COLUMNS].to_numpy(float)
        arc_xyz = arc_rows[COORD_COLUMNS].to_numpy(float)
        pairs = match_pairs(actin_xyz, arc_xyz, threshold_um)

        for actin_idx, arc_idx in pairs:
            actin_xyz_point = actin_xyz[actin_idx]
            arc_xyz_point = arc_xyz[arc_idx]
            distance_um = float(np.linalg.norm(actin_xyz_point - arc_xyz_point))
            current_pair_id = pair_id
            pair_id += 1

            pair_rows.append(
                {
                    "pair_id": current_pair_id,
                    "image_key": image_key,
                    "actin_source_file": str(actin_rows.loc[actin_idx, SOURCE_COLUMN]),
                    "arc_source_file": str(arc_rows.loc[arc_idx, SOURCE_COLUMN]),
                    "actin_index_within_image": int(actin_idx),
                    "arc_index_within_image": int(arc_idx),
                    "actin_x_um": float(actin_xyz_point[0]),
                    "actin_y_um": float(actin_xyz_point[1]),
                    "actin_z_um": float(actin_xyz_point[2]),
                    "arc_x_um": float(arc_xyz_point[0]),
                    "arc_y_um": float(arc_xyz_point[1]),
                    "arc_z_um": float(arc_xyz_point[2]),
                    "distance_um": distance_um,
                }
            )

            puncta_rows.append(
                {
                    "pair_id": current_pair_id,
                    "image_key": image_key,
                    "channel": "actin",
                    "source_file": str(actin_rows.loc[actin_idx, SOURCE_COLUMN]),
                    "index_within_image": int(actin_idx),
                    "x_um": float(actin_xyz_point[0]),
                    "y_um": float(actin_xyz_point[1]),
                    "z_um": float(actin_xyz_point[2]),
                    "distance_um_to_partner": distance_um,
                }
            )
            puncta_rows.append(
                {
                    "pair_id": current_pair_id,
                    "image_key": image_key,
                    "channel": "arc",
                    "source_file": str(arc_rows.loc[arc_idx, SOURCE_COLUMN]),
                    "index_within_image": int(arc_idx),
                    "x_um": float(arc_xyz_point[0]),
                    "y_um": float(arc_xyz_point[1]),
                    "z_um": float(arc_xyz_point[2]),
                    "distance_um_to_partner": distance_um,
                }
            )

    pairs_df = pd.DataFrame.from_records(pair_rows).sort_values(["image_key", "pair_id"]).reset_index(drop=True)
    puncta_df = pd.DataFrame.from_records(puncta_rows).sort_values(["image_key", "pair_id", "channel"]).reset_index(drop=True)
    pairs_df.to_csv(pairs_output_path, index=False)
    puncta_df.to_csv(puncta_output_path, index=False)
    return len(pairs_df), len(puncta_df)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = load_summary(args.summary_csv, args.threshold_um, args.null_model)
    plot_whole_image_relative(summary, args.output_dir / "coassociation_whole_image_relative_100nm.png")
    plot_binned_relative(summary, args.output_dir / "coassociation_bins_relative_100nm.png")

    actin_df = pd.read_csv(args.actin_csv)
    arc_df = pd.read_csv(args.arc_csv)
    pair_count, puncta_count = export_coassociated_coordinate_tables(
        actin_df=actin_df,
        arc_df=arc_df,
        threshold_um=args.threshold_um,
        pairs_output_path=args.output_dir / "coassociated_pairs_100nm.csv",
        puncta_output_path=args.output_dir / "coassociated_puncta_coordinates_100nm.csv",
    )

    print(f"Wrote extended plots and tables to {args.output_dir}")
    print(f"Co-associated pairs: {pair_count}")
    print(f"Co-associated puncta rows: {puncta_count}")


if __name__ == "__main__":
    main()
