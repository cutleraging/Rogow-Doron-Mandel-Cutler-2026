#!/usr/bin/env python3

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BINS = [
    ("0_40", 0.0, 40.0),
    ("40_80", 40.0, 80.0),
    ("80_120", 80.0, 120.0),
    ("120_160", 120.0, 160.0),
]


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


def overlap_matrix(actin_points: np.ndarray, arc_points: np.ndarray, threshold_um: float) -> np.ndarray:
    if len(actin_points) == 0 or len(arc_points) == 0:
        return np.zeros((len(actin_points), len(arc_points)), dtype=bool)
    deltas = actin_points[:, None, :] - arc_points[None, :, :]
    return np.sum(deltas * deltas, axis=2) <= (threshold_um ** 2)


def bin_label(x_um: float) -> str:
    for name, low, high in BINS:
        if low <= x_um < high:
            return name
    return "outside_bins"


def compute_pair_details(actin_points: np.ndarray, arc_points: np.ndarray, threshold_um: float) -> dict[str, object]:
    within = overlap_matrix(actin_points, arc_points, threshold_um)
    actin_assoc = within.any(axis=1) if len(actin_points) else np.zeros(0, dtype=bool)
    arc_assoc = within.any(axis=0) if len(arc_points) else np.zeros(0, dtype=bool)
    pair_idx = np.argwhere(within)
    pair_rows: list[dict[str, object]] = []
    for actin_i, arc_i in pair_idx:
        ax, ay = actin_points[actin_i]
        rx, ry = arc_points[arc_i]
        pair_rows.append(
            {
                "actin_index": int(actin_i),
                "arc_index": int(arc_i),
                "actin_x_um": float(ax),
                "actin_y_um": float(ay),
                "arc_x_um": float(rx),
                "arc_y_um": float(ry),
                "distance_um": float(np.hypot(ax - rx, ay - ry)),
                "actin_bin": bin_label(float(ax)),
                "arc_bin": bin_label(float(rx)),
            }
        )
    return {
        "within": within,
        "actin_assoc": actin_assoc,
        "arc_assoc": arc_assoc,
        "pairs": pd.DataFrame(pair_rows),
    }


def sanitize_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def per_bin_frequency(points_a: np.ndarray, points_b: np.ndarray, threshold_um: float) -> dict[str, float]:
    output: dict[str, float] = {}
    for name, low, high in BINS:
        actin_subset = points_a[(points_a[:, 0] >= low) & (points_a[:, 0] < high)]
        arc_subset = points_b[(points_b[:, 0] >= low) & (points_b[:, 0] < high)]
        within = overlap_matrix(actin_subset, arc_subset, threshold_um)
        assoc_total = int(within.any(axis=1).sum()) + int(within.any(axis=0).sum())
        total_puncta = len(actin_subset) + len(arc_subset)
        output[name] = assoc_total / total_puncta if total_puncta else 0.0
    return output


def plot_image(
    image_id: str,
    actin_points: np.ndarray,
    arc_points: np.ndarray,
    threshold_um: float,
    output_path: Path,
) -> dict[str, object]:
    pair_details = compute_pair_details(actin_points, arc_points, threshold_um)
    pairs = pair_details["pairs"]
    actin_assoc = pair_details["actin_assoc"]
    arc_assoc = pair_details["arc_assoc"]

    bin_freqs = per_bin_frequency(actin_points, arc_points, threshold_um)
    total_assoc = int(actin_assoc.sum()) + int(arc_assoc.sum())
    total_puncta = len(actin_points) + len(arc_points)
    whole_freq = total_assoc / total_puncta if total_puncta else 0.0

    max_x = 160.0
    if len(actin_points):
        max_x = max(max_x, float(np.max(actin_points[:, 0])) + 2.0)
    if len(arc_points):
        max_x = max(max_x, float(np.max(arc_points[:, 0])) + 2.0)
    all_y = np.concatenate(
        [arr[:, 1] for arr in (actin_points, arc_points) if len(arr)],
        axis=0,
    ) if (len(actin_points) or len(arc_points)) else np.array([0.0])
    y_min = float(np.min(all_y) - 1.0)
    y_max = float(np.max(all_y) + 1.0)
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.linewidth": 0.8,
        }
    )
    fig, ax = plt.subplots(figsize=(10.5, 4.8))

    for boundary in (40, 80, 120, 160):
        ax.axvline(boundary, color="#d9d9d9", linestyle="--", linewidth=0.8, zorder=0)

    if len(pairs):
        for row in pairs.itertuples(index=False):
            ax.plot(
                [row.actin_x_um, row.arc_x_um],
                [row.actin_y_um, row.arc_y_um],
                color="#9a9a9a",
                linewidth=0.8,
                alpha=0.8,
                zorder=1,
            )

    if len(actin_points):
        ax.scatter(
            actin_points[:, 0],
            actin_points[:, 1],
            s=28,
            facecolors="none",
            edgecolors="#1f77b4",
            linewidths=1.0,
            alpha=0.9,
            zorder=2,
            label=f"Actin (n={len(actin_points)})",
        )
        if actin_assoc.any():
            ax.scatter(
                actin_points[actin_assoc, 0],
                actin_points[actin_assoc, 1],
                s=52,
                facecolors="#1f77b4",
                edgecolors="white",
                linewidths=0.5,
                zorder=4,
                label=f"Actin assoc (n={int(actin_assoc.sum())})",
            )

    if len(arc_points):
        ax.scatter(
            arc_points[:, 0],
            arc_points[:, 1],
            s=28,
            marker="s",
            facecolors="none",
            edgecolors="#d62728",
            linewidths=1.0,
            alpha=0.9,
            zorder=2,
            label=f"Arc (n={len(arc_points)})",
        )
        if arc_assoc.any():
            ax.scatter(
                arc_points[arc_assoc, 0],
                arc_points[arc_assoc, 1],
                s=52,
                marker="s",
                facecolors="#d62728",
                edgecolors="white",
                linewidths=0.5,
                zorder=4,
                label=f"Arc assoc (n={int(arc_assoc.sum())})",
            )

    ax.set_xlim(-1, max_x)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("Distance from soma (um)")
    ax.set_ylabel("Y position (um)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.8)
    title = (
        f"{image_id}\n"
        f"Whole: {whole_freq * 100:.1f}% ({total_assoc}/{total_puncta}) | "
        f"0-40: {bin_freqs['0_40'] * 100:.1f}% | "
        f"40-80: {bin_freqs['40_80'] * 100:.1f}% | "
        f"80-120: {bin_freqs['80_120'] * 100:.1f}% | "
        f"120-160: {bin_freqs['120_160'] * 100:.1f}%"
    )
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, loc="upper right", fontsize=8, ncol=2, handlelength=1.4, columnspacing=1.0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=250, bbox_inches="tight")
    plt.close(fig)

    return {
        "image_id": image_id,
        "whole_coassoc_freq": whole_freq,
        "80_120_coassoc_freq": bin_freqs["80_120"],
        "actin_n": int(len(actin_points)),
        "arc_n": int(len(arc_points)),
        "actin_assoc_n": int(actin_assoc.sum()),
        "arc_assoc_n": int(arc_assoc.sum()),
        "pair_count": int(len(pairs)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--actin", default="../../1-data/batch-2/actin puncta locations.csv")
    parser.add_argument("--arc", default="../../1-data/batch-2/arc puncta locations.csv")
    parser.add_argument("--threshold-um", type=float, default=0.25)
    parser.add_argument("--outdir", default="results_250nm/image_puncta_maps")
    parser.add_argument("--pair-details-csv", default="results_250nm/image_pair_details.csv")
    parser.add_argument("--image-summary-csv", default="results_250nm/image_plot_summary.csv")
    args = parser.parse_args()

    actin_df = load_puncta_table(args.actin)
    arc_df = load_puncta_table(args.arc)
    shared_images = sorted(set(actin_df["image_id"]) & set(arc_df["image_id"]))
    if not shared_images:
        raise SystemExit("No shared images after canonicalizing source file names.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    image_rows: list[dict[str, object]] = []
    pair_rows: list[pd.DataFrame] = []
    for image_id in shared_images:
        actin_points = actin_df.loc[actin_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)
        arc_points = arc_df.loc[arc_df["image_id"] == image_id, ["x_um", "y_um"]].to_numpy(float)

        image_rows.append(
            plot_image(
                image_id=image_id,
                actin_points=actin_points,
                arc_points=arc_points,
                threshold_um=args.threshold_um,
                output_path=outdir / f"{sanitize_filename(image_id)}.png",
            )
        )

        pair_df = compute_pair_details(actin_points, arc_points, args.threshold_um)["pairs"]
        if len(pair_df):
            pair_df.insert(0, "image_id", image_id)
            pair_rows.append(pair_df)

    pd.DataFrame(image_rows).sort_values("image_id").to_csv(args.image_summary_csv, index=False)
    if pair_rows:
        pd.concat(pair_rows, ignore_index=True).to_csv(args.pair_details_csv, index=False)
    else:
        pd.DataFrame(
            columns=[
                "image_id",
                "actin_index",
                "arc_index",
                "actin_x_um",
                "actin_y_um",
                "arc_x_um",
                "arc_y_um",
                "distance_um",
                "actin_bin",
                "arc_bin",
            ]
        ).to_csv(args.pair_details_csv, index=False)


if __name__ == "__main__":
    main()
