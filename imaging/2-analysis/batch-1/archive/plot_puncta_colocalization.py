from __future__ import annotations

import argparse
import re
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

SOURCE_COLUMN = "Source File"
COORD_COLUMNS = ["X (um)", "Y (um)", "Z (um)"]
BIN_ORDER = ["0_30", "30_60", "60_90", "90_120"]


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
    parser = argparse.ArgumentParser(description="Create puncta co-association figures.")
    parser.add_argument("--summary-csv", type=Path, required=True, help="Path to analysis summary CSV.")
    parser.add_argument("--actin-csv", type=Path, required=True, help="Path to actin puncta CSV.")
    parser.add_argument("--arc-csv", type=Path, required=True, help="Path to arc puncta CSV.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory for figures.")
    parser.add_argument("--threshold-um", type=float, default=0.1, help="Co-association threshold in um.")
    parser.add_argument("--null-model", type=str, default="fix_x_permute_yz", help="Null model to display.")
    parser.add_argument(
        "--example-image-key",
        type=str,
        default=None,
        help="Normalized image key to use for the example puncta map. Defaults to the image with most matches.",
    )
    return parser.parse_args()


def load_summary(summary_path: Path, threshold_um: float, null_model: str) -> pd.DataFrame:
    summary = pd.read_csv(summary_path)
    if "null_model" in summary.columns:
        summary = summary.loc[summary["null_model"] == null_model]
    summary = summary.loc[np.isclose(summary["threshold_um"], threshold_um)]
    if summary.empty:
        raise ValueError(
            f"No rows found in summary for threshold={threshold_um} and null_model={null_model}."
        )
    return summary.copy()


def plot_whole_image(summary: pd.DataFrame, output_path: Path) -> None:
    row = summary.loc[summary["region"] == "whole_image"].iloc[0]

    observed_pct = 100.0 * row["arc_match_fraction"]
    null_pct = 100.0 * row["null_arc_match_fraction_mean"]
    null_low = 100.0 * row["null_pairs_ci_low"] / row["arc_total"]
    null_high = 100.0 * row["null_pairs_ci_high"] / row["arc_total"]

    fig, ax = plt.subplots(figsize=(6, 5), dpi=200)
    x = np.arange(2)
    bars = ax.bar(
        x,
        [observed_pct, null_pct],
        color=["#0b7285", "#bfc6d0"],
        edgecolor="#222222",
        width=0.65,
        linewidth=1.0,
    )
    yerr_low = max(null_pct - null_low, 0.0)
    yerr_high = max(null_high - null_pct, 0.0)
    ax.errorbar(
        x[1],
        null_pct,
        yerr=np.array([[yerr_low], [yerr_high]]),
        fmt="none",
        ecolor="#222222",
        elinewidth=1.2,
        capsize=4,
        capthick=1.2,
        zorder=3,
    )

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height + 1.2, f"{height:.1f}%", ha="center", va="bottom", fontsize=10)

    ax.set_xticks(x, ["Observed", "Null mean"])
    ax.set_ylabel("Arc puncta co-associated within 100 nm (%)")
    ax.set_title("Whole Image Co-association")
    ax.set_ylim(0, max(observed_pct, null_high) * 1.25 + 2)
    ax.grid(axis="y", alpha=0.3)

    annotation = (
        f"Observed pairs: {int(row['observed_pairs'])}/{int(row['arc_total'])} arc\n"
        f"Enrichment: {row['enrichment_vs_null']:.2f}x\n"
        f"Empirical p: {row['empirical_p_enrichment']:.4f}"
    )
    ax.text(
        0.03,
        0.97,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#dddddd"},
    )

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def plot_bins(summary: pd.DataFrame, output_path: Path) -> None:
    bins = summary.loc[summary["region"].isin(BIN_ORDER)].copy()
    bins["region"] = pd.Categorical(bins["region"], categories=BIN_ORDER, ordered=True)
    bins = bins.sort_values("region")

    observed_pct = 100.0 * bins["arc_match_fraction"].to_numpy()
    null_pct = 100.0 * bins["null_arc_match_fraction_mean"].to_numpy()
    null_low = 100.0 * (bins["null_pairs_ci_low"] / bins["arc_total"]).to_numpy()
    null_high = 100.0 * (bins["null_pairs_ci_high"] / bins["arc_total"]).to_numpy()

    fig, ax = plt.subplots(figsize=(8, 5), dpi=200)
    x = np.arange(len(BIN_ORDER))
    width = 0.36

    ax.bar(x - width / 2, observed_pct, width=width, color="#006d77", edgecolor="#222222", linewidth=1.0, label="Observed")
    ax.bar(x + width / 2, null_pct, width=width, color="#bfc6d0", edgecolor="#222222", linewidth=1.0, label="Null mean")
    yerr_low = np.maximum(null_pct - null_low, 0.0)
    yerr_high = np.maximum(null_high - null_pct, 0.0)
    ax.errorbar(
        x + width / 2,
        null_pct,
        yerr=np.vstack((yerr_low, yerr_high)),
        fmt="none",
        ecolor="#222222",
        elinewidth=1.2,
        capsize=3,
        capthick=1.2,
        zorder=3,
    )

    for idx, (_, row) in enumerate(bins.iterrows()):
        ax.text(
            x[idx] - width / 2,
            observed_pct[idx] + 0.8,
            f"{int(row['observed_pairs'])}/{int(row['arc_total'])}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )

    ax.set_xticks(x, ["0-30", "30-60", "60-90", "90-120"])
    ax.set_xlabel("Distance from soma (um)")
    ax.set_ylabel("Arc puncta co-associated within 100 nm (%)")
    ax.set_title("Distance-Binned Co-association")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(frameon=False, loc="upper right")
    ax.set_ylim(0, max(np.max(observed_pct), np.max(null_high)) * 1.3 + 2)

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)


def choose_example_image(
    actin_df: pd.DataFrame,
    arc_df: pd.DataFrame,
    threshold_um: float,
    example_image_key: str | None,
) -> str:
    image_keys = sorted(set(actin_df["image_key"]) & set(arc_df["image_key"]))
    if example_image_key is not None:
        if example_image_key not in image_keys:
            raise ValueError(f"Requested --example-image-key not found: {example_image_key}")
        return example_image_key

    scored: list[tuple[float, int, int, str]] = []
    for key in image_keys:
        actin_coords = actin_df.loc[actin_df["image_key"] == key, COORD_COLUMNS].to_numpy(float)
        arc_coords = arc_df.loc[arc_df["image_key"] == key, COORD_COLUMNS].to_numpy(float)
        match_count = len(match_pairs(actin_coords, arc_coords, threshold_um))
        total_points = len(actin_coords) + len(arc_coords)
        arc_fraction = match_count / len(arc_coords) if len(arc_coords) else 0.0
        scored.append((arc_fraction, match_count, -total_points, key))

    readable = [item for item in scored if -item[2] <= 80 and item[1] >= 5]
    candidates = readable if readable else scored
    candidates.sort(reverse=True)
    return candidates[0][3]


def plot_example_image(
    actin_df: pd.DataFrame,
    arc_df: pd.DataFrame,
    threshold_um: float,
    image_key: str,
    output_path: Path,
    pairs_csv_path: Path,
) -> None:
    actin = actin_df.loc[actin_df["image_key"] == image_key, COORD_COLUMNS].to_numpy(float)
    arc = arc_df.loc[arc_df["image_key"] == image_key, COORD_COLUMNS].to_numpy(float)
    pairs = match_pairs(actin, arc, threshold_um)

    actin_match_idx = np.array([pair[0] for pair in pairs], dtype=int) if pairs else np.array([], dtype=int)
    arc_match_idx = np.array([pair[1] for pair in pairs], dtype=int) if pairs else np.array([], dtype=int)
    ref_pair = None
    if pairs:
        ref_pair = min(
            pairs,
            key=lambda pair: float(np.linalg.norm(actin[pair[0]] - arc[pair[1]])),
        )

    fig = plt.figure(figsize=(10.5, 7.5), dpi=220)
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(actin[:, 0], actin[:, 1], actin[:, 2], s=28, color="#1d4e89", alpha=0.60, label="Actin puncta")
    ax.scatter(arc[:, 0], arc[:, 1], arc[:, 2], s=38, color="#d62828", alpha=0.70, marker="^", label="Arc puncta")

    if len(actin_match_idx):
        ax.scatter(
            actin[actin_match_idx, 0],
            actin[actin_match_idx, 1],
            actin[actin_match_idx, 2],
            s=85,
            facecolor="#72efdd",
            edgecolor="#111111",
            linewidth=0.8,
            label="Actin in matched pairs",
            zorder=4,
        )
    if len(arc_match_idx):
        ax.scatter(
            arc[arc_match_idx, 0],
            arc[arc_match_idx, 1],
            arc[arc_match_idx, 2],
            s=95,
            facecolor="#ffd166",
            edgecolor="#111111",
            linewidth=0.8,
            marker="^",
            label="Arc in matched pairs",
            zorder=4,
        )

    for left_idx, right_idx in pairs:
        ax.plot(
            [actin[left_idx, 0], arc[right_idx, 0]],
            [actin[left_idx, 1], arc[right_idx, 1]],
            [actin[left_idx, 2], arc[right_idx, 2]],
            color="#2a9d8f",
            linewidth=1.8,
            alpha=0.95,
            zorder=3,
        )

    xmin = float(min(np.min(actin[:, 0]), np.min(arc[:, 0])))
    xmax = float(max(np.max(actin[:, 0]), np.max(arc[:, 0])))
    ymin = float(min(np.min(actin[:, 1]), np.min(arc[:, 1])))
    ymax = float(max(np.max(actin[:, 1]), np.max(arc[:, 1])))
    zmin = float(min(np.min(actin[:, 2]), np.min(arc[:, 2])))
    zmax = float(max(np.max(actin[:, 2]), np.max(arc[:, 2])))
    xspan = xmax - xmin if xmax > xmin else 1.0
    yspan = ymax - ymin if ymax > ymin else 1.0
    zspan = zmax - zmin if zmax > zmin else 1.0

    ax.set_xlim(xmin - 0.04 * xspan, xmax + 0.04 * xspan)
    ax.set_ylim(ymin - 0.08 * yspan, ymax + 0.08 * yspan)
    ax.set_zlim(zmin - 0.12 * zspan, zmax + 0.12 * zspan)
    # Use a display aspect that keeps the long x-axis visible but still shows y/z depth.
    ax.set_box_aspect((2.4, 1.2, 1.0))
    ax.set_xlabel("X position from soma (um)")
    ax.set_ylabel("Y position (um)")
    ax.set_zlabel("Z position (um)")
    ax.view_init(elev=27, azim=-58)
    ax.set_title(
        f"Example Image: {image_key}\n"
        f"3D matched pairs at <=100 nm: {len(pairs)} of {len(arc)} arc puncta"
    )
    ax.grid(alpha=0.2)

    if ref_pair is not None:
        ref_left, ref_right = ref_pair
        center_x, center_y, center_z = arc[ref_right]
        u = np.linspace(0, 2 * np.pi, 32)
        v = np.linspace(0, np.pi, 16)
        sphere_x = center_x + threshold_um * np.outer(np.cos(u), np.sin(v))
        sphere_y = center_y + threshold_um * np.outer(np.sin(u), np.sin(v))
        sphere_z = center_z + threshold_um * np.outer(np.ones_like(u), np.cos(v))
        ax.plot_wireframe(
            sphere_x,
            sphere_y,
            sphere_z,
            rstride=2,
            cstride=2,
            color="#111111",
            linewidth=0.6,
            alpha=0.35,
        )
        distance_nm = float(np.linalg.norm(actin[ref_left] - arc[ref_right]) * 1000.0)
        ax.text(
            center_x,
            center_y,
            center_z + 1.25 * threshold_um,
            f"100 nm sphere\nref pair: {distance_nm:.1f} nm",
            fontsize=8,
            ha="center",
            va="bottom",
            color="#111111",
        )

    bar_x0 = xmin + 0.06 * xspan
    bar_x1 = bar_x0 + threshold_um
    bar_y = ymin + 0.05 * yspan
    bar_z = zmin + 0.05 * zspan
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y], [bar_z, bar_z], color="#111111", linewidth=3.0)
    ax.text((bar_x0 + bar_x1) / 2, bar_y, bar_z + 0.04 * zspan, "100 nm", fontsize=9, ha="center")
    ax.legend(frameon=True, framealpha=0.95, fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(output_path)
    plt.close(fig)

    pair_rows = []
    for left_idx, right_idx in pairs:
        pair_rows.append(
            {
                "image_key": image_key,
                "actin_index": int(left_idx),
                "arc_index": int(right_idx),
                "actin_x_um": float(actin[left_idx, 0]),
                "actin_y_um": float(actin[left_idx, 1]),
                "actin_z_um": float(actin[left_idx, 2]),
                "arc_x_um": float(arc[right_idx, 0]),
                "arc_y_um": float(arc[right_idx, 1]),
                "arc_z_um": float(arc[right_idx, 2]),
                "distance_um": float(np.linalg.norm(actin[left_idx] - arc[right_idx])),
            }
        )
    pd.DataFrame(pair_rows).to_csv(pairs_csv_path, index=False)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = load_summary(args.summary_csv, args.threshold_um, args.null_model)
    plot_whole_image(summary, args.output_dir / "coassociation_whole_image_100nm.png")
    plot_bins(summary, args.output_dir / "coassociation_bins_100nm.png")

    actin_df = pd.read_csv(args.actin_csv)
    arc_df = pd.read_csv(args.arc_csv)
    actin_df["image_key"] = actin_df[SOURCE_COLUMN].map(normalize_source_name)
    arc_df["image_key"] = arc_df[SOURCE_COLUMN].map(normalize_source_name)

    selected_key = choose_example_image(
        actin_df=actin_df,
        arc_df=arc_df,
        threshold_um=args.threshold_um,
        example_image_key=args.example_image_key,
    )
    plot_example_image(
        actin_df=actin_df,
        arc_df=arc_df,
        threshold_um=args.threshold_um,
        image_key=selected_key,
        output_path=args.output_dir / "example_puncta_map_100nm.png",
        pairs_csv_path=args.output_dir / "example_puncta_pairs_100nm.csv",
    )

    print(f"Wrote figures to {args.output_dir}")
    print(f"Example image key: {selected_key}")


if __name__ == "__main__":
    main()
