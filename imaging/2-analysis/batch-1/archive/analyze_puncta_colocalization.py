from __future__ import annotations

import argparse
import math
import re
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


COORD_COLUMNS = ["X (um)", "Y (um)", "Z (um)"]
SOURCE_COLUMN = "Source File"
DEFAULT_THRESHOLDS = (0.25, 0.5, 1.0)
NULL_MODELS = (
    "fix_x_permute_yz",
    "permute_x_fix_yz",
    "shift_x_circular",
    "shuffle_arc_between_images",
)
NULL_MODEL_DESCRIPTIONS = {
    "fix_x_permute_yz": (
        "Within each image and channel, x is fixed while y and z are independently permuted."
    ),
    "permute_x_fix_yz": (
        "Within each image and channel, x is permuted while y and z are fixed."
    ),
    "shift_x_circular": (
        "Within each image and channel, all x values are circularly shifted by a random offset."
    ),
    "shuffle_arc_between_images": (
        "Actin stays in each image while arc coordinates are reassigned by permuting image identities."
    ),
}


@dataclass(frozen=True)
class Region:
    label: str
    xmin: float
    xmax: float

    def mask(self, x_values: np.ndarray) -> np.ndarray:
        mask = x_values >= self.xmin
        if math.isfinite(self.xmax):
            mask &= x_values < self.xmax
        return mask


@dataclass
class ImageData:
    image_key: str
    actin_coords: np.ndarray
    arc_coords: np.ndarray


def normalize_source_name(value: str) -> str:
    return re.sub(r"^C\d-", "", value)


def hopcroft_karp_size(adjacency: list[list[int]], n_right: int) -> int:
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

    matching_size = 0
    while bfs():
        for left_index in range(n_left):
            if pair_left[left_index] == -1 and dfs(left_index):
                matching_size += 1
    return matching_size


def match_counts_by_threshold(
    actin_coords: np.ndarray,
    arc_coords: np.ndarray,
    thresholds: tuple[float, ...],
) -> list[int]:
    if len(actin_coords) == 0 or len(arc_coords) == 0:
        return [0] * len(thresholds)

    left_coords = actin_coords
    right_coords = arc_coords
    if len(left_coords) > len(right_coords):
        left_coords, right_coords = right_coords, left_coords

    deltas = left_coords[:, None, :] - right_coords[None, :, :]
    distance_sq = np.sum(deltas * deltas, axis=2)
    sorted_right = np.argsort(distance_sq, axis=1)
    sorted_dist_sq = np.take_along_axis(distance_sq, sorted_right, axis=1)

    counts: list[int] = []
    for threshold in thresholds:
        threshold_sq = threshold * threshold
        adjacency: list[list[int]] = []
        for row_index in range(sorted_dist_sq.shape[0]):
            edge_count = int(np.searchsorted(sorted_dist_sq[row_index], threshold_sq, side="right"))
            adjacency.append(sorted_right[row_index, :edge_count].tolist())
        counts.append(hopcroft_karp_size(adjacency, len(right_coords)))
    return counts


def load_tables(actin_path: Path, arc_path: Path) -> tuple[list[ImageData], pd.DataFrame, pd.DataFrame]:
    actin_df = pd.read_csv(actin_path)
    arc_df = pd.read_csv(arc_path)
    actin_df["image_key"] = actin_df[SOURCE_COLUMN].map(normalize_source_name)
    arc_df["image_key"] = arc_df[SOURCE_COLUMN].map(normalize_source_name)

    actin_keys = set(actin_df["image_key"])
    arc_keys = set(arc_df["image_key"])
    if actin_keys != arc_keys:
        missing_from_actin = sorted(arc_keys - actin_keys)
        missing_from_arc = sorted(actin_keys - arc_keys)
        raise ValueError(
            "Image keys do not align after normalization. "
            f"Missing from actin: {missing_from_actin[:5]}; missing from arc: {missing_from_arc[:5]}"
        )

    images: list[ImageData] = []
    for image_key in sorted(actin_keys):
        actin_coords = (
            actin_df.loc[actin_df["image_key"] == image_key, COORD_COLUMNS]
            .to_numpy(dtype=float, copy=True)
        )
        arc_coords = (
            arc_df.loc[arc_df["image_key"] == image_key, COORD_COLUMNS]
            .to_numpy(dtype=float, copy=True)
        )
        images.append(
            ImageData(
                image_key=image_key,
                actin_coords=actin_coords,
                arc_coords=arc_coords,
            )
        )
    return images, actin_df, arc_df


def circular_shift_x(coords: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    shifted = coords.copy()
    if len(shifted) < 2:
        return shifted
    xmin = float(np.min(shifted[:, 0]))
    xmax = float(np.max(shifted[:, 0]))
    span = xmax - xmin
    if span <= 0:
        return shifted
    offset = rng.uniform(0.0, span)
    shifted[:, 0] = ((shifted[:, 0] - xmin + offset) % span) + xmin
    return shifted


def simulate_null_coords_for_image(
    images: list[ImageData],
    image_index: int,
    null_model: str,
    rng: np.random.Generator,
    arc_image_permutation: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    image = images[image_index]

    if null_model == "shuffle_arc_between_images":
        if arc_image_permutation is None:
            raise ValueError("arc_image_permutation is required for shuffle_arc_between_images.")
        return image.actin_coords, images[int(arc_image_permutation[image_index])].arc_coords

    actin_coords = image.actin_coords.copy()
    arc_coords = image.arc_coords.copy()

    if null_model == "fix_x_permute_yz":
        if len(actin_coords) > 1:
            actin_coords[:, 1] = actin_coords[rng.permutation(len(actin_coords)), 1]
            actin_coords[:, 2] = actin_coords[rng.permutation(len(actin_coords)), 2]
        if len(arc_coords) > 1:
            arc_coords[:, 1] = arc_coords[rng.permutation(len(arc_coords)), 1]
            arc_coords[:, 2] = arc_coords[rng.permutation(len(arc_coords)), 2]
        return actin_coords, arc_coords

    if null_model == "permute_x_fix_yz":
        if len(actin_coords) > 1:
            actin_coords[:, 0] = actin_coords[rng.permutation(len(actin_coords)), 0]
        if len(arc_coords) > 1:
            arc_coords[:, 0] = arc_coords[rng.permutation(len(arc_coords)), 0]
        return actin_coords, arc_coords

    if null_model == "shift_x_circular":
        return circular_shift_x(actin_coords, rng), circular_shift_x(arc_coords, rng)

    raise ValueError(f"Unknown null model: {null_model}")


def summarize_counts(
    images: list[ImageData],
    regions: list[Region],
    thresholds: tuple[float, ...],
    simulations: int,
    seed: int,
    null_model: str,
) -> pd.DataFrame:
    region_labels = [region.label for region in regions]
    region_lookup = {region.label: region for region in regions}
    region_totals = {
        label: {
            "actin_total": int(sum(region.mask(image.actin_coords[:, 0]).sum() for image in images)),
            "arc_total": int(sum(region.mask(image.arc_coords[:, 0]).sum() for image in images)),
        }
        for label, region in region_lookup.items()
    }

    observed_pair_counts = {
        label: np.zeros(len(thresholds), dtype=int) for label in region_labels
    }
    for image in images:
        for label in region_labels:
            region = region_lookup[label]
            actin_mask = region.mask(image.actin_coords[:, 0])
            arc_mask = region.mask(image.arc_coords[:, 0])
            counts = match_counts_by_threshold(
                image.actin_coords[actin_mask],
                image.arc_coords[arc_mask],
                thresholds,
            )
            observed_pair_counts[label] += np.array(counts, dtype=int)

    rng = np.random.default_rng(seed)
    null_pair_counts = {
        label: np.zeros((simulations, len(thresholds)), dtype=int) for label in region_labels
    }
    for simulation_index in range(simulations):
        arc_image_permutation = None
        if null_model == "shuffle_arc_between_images":
            arc_image_permutation = rng.permutation(len(images))

        for image_index, image in enumerate(images):
            actin_coords, arc_coords = simulate_null_coords_for_image(
                images=images,
                image_index=image_index,
                null_model=null_model,
                rng=rng,
                arc_image_permutation=arc_image_permutation,
            )
            for label in region_labels:
                region = region_lookup[label]
                actin_mask = region.mask(actin_coords[:, 0])
                arc_mask = region.mask(arc_coords[:, 0])
                counts = match_counts_by_threshold(
                    actin_coords[actin_mask],
                    arc_coords[arc_mask],
                    thresholds,
                )
                null_pair_counts[label][simulation_index] += np.array(counts, dtype=int)

    records: list[dict[str, float | int | str]] = []
    for label in region_labels:
        actin_total = region_totals[label]["actin_total"]
        arc_total = region_totals[label]["arc_total"]
        null_counts = null_pair_counts[label]
        region = region_lookup[label]

        for threshold_index, threshold in enumerate(thresholds):
            observed_pairs = int(observed_pair_counts[label][threshold_index])
            null_pairs = null_counts[:, threshold_index]
            null_mean = float(np.mean(null_pairs))
            null_low, null_high = np.quantile(null_pairs, [0.025, 0.975])

            actin_fraction = observed_pairs / actin_total if actin_total else float("nan")
            arc_fraction = observed_pairs / arc_total if arc_total else float("nan")
            dice_fraction = (
                2 * observed_pairs / (actin_total + arc_total)
                if (actin_total + arc_total)
                else float("nan")
            )

            null_actin_mean = null_mean / actin_total if actin_total else float("nan")
            null_arc_mean = null_mean / arc_total if arc_total else float("nan")
            null_dice_mean = (
                2 * null_mean / (actin_total + arc_total)
                if (actin_total + arc_total)
                else float("nan")
            )
            null_actin_ci_low = null_low / actin_total if actin_total else float("nan")
            null_actin_ci_high = null_high / actin_total if actin_total else float("nan")
            null_arc_ci_low = null_low / arc_total if arc_total else float("nan")
            null_arc_ci_high = null_high / arc_total if arc_total else float("nan")
            null_all_ci_low = (
                2 * null_low / (actin_total + arc_total)
                if (actin_total + arc_total)
                else float("nan")
            )
            null_all_ci_high = (
                2 * null_high / (actin_total + arc_total)
                if (actin_total + arc_total)
                else float("nan")
            )
            all_puncta_total = actin_total + arc_total

            empirical_p = (1 + int(np.sum(null_pairs >= observed_pairs))) / (simulations + 1)
            enrichment = observed_pairs / null_mean if null_mean else float("inf")

            records.append(
                {
                    "threshold_um": threshold,
                    "null_model": null_model,
                    "region": label,
                    "xmin_um": region.xmin,
                    "xmax_um": region.xmax,
                    "actin_total": actin_total,
                    "arc_total": arc_total,
                    "all_puncta_total": all_puncta_total,
                    "observed_pairs": observed_pairs,
                    "coassociated_actin_observed": observed_pairs,
                    "coassociated_arc_observed": observed_pairs,
                    "coassociated_all_puncta_observed": 2 * observed_pairs,
                    "actin_coassociation_fraction": actin_fraction,
                    "arc_coassociation_fraction": arc_fraction,
                    "all_puncta_coassociation_fraction": dice_fraction,
                    "actin_match_fraction": actin_fraction,
                    "arc_match_fraction": arc_fraction,
                    "dice_fraction": dice_fraction,
                    "null_pairs_mean": null_mean,
                    "null_pairs_ci_low": float(null_low),
                    "null_pairs_ci_high": float(null_high),
                    "null_coassociated_actin_mean": null_mean,
                    "null_coassociated_arc_mean": null_mean,
                    "null_coassociated_all_puncta_mean": 2 * null_mean,
                    "null_coassociated_all_puncta_ci_low": 2 * float(null_low),
                    "null_coassociated_all_puncta_ci_high": 2 * float(null_high),
                    "null_actin_coassociation_fraction_mean": null_actin_mean,
                    "null_arc_coassociation_fraction_mean": null_arc_mean,
                    "null_all_puncta_coassociation_fraction_mean": null_dice_mean,
                    "null_actin_coassociation_fraction_ci_low": null_actin_ci_low,
                    "null_actin_coassociation_fraction_ci_high": null_actin_ci_high,
                    "null_arc_coassociation_fraction_ci_low": null_arc_ci_low,
                    "null_arc_coassociation_fraction_ci_high": null_arc_ci_high,
                    "null_all_puncta_coassociation_fraction_ci_low": null_all_ci_low,
                    "null_all_puncta_coassociation_fraction_ci_high": null_all_ci_high,
                    "null_actin_match_fraction_mean": null_actin_mean,
                    "null_arc_match_fraction_mean": null_arc_mean,
                    "null_dice_fraction_mean": null_dice_mean,
                    "enrichment_vs_null": enrichment,
                    "empirical_p_enrichment": empirical_p,
                }
            )
    return pd.DataFrame.from_records(records)


def build_regions() -> list[Region]:
    return [
        Region(label="whole_image", xmin=-math.inf, xmax=math.inf),
        Region(label="x_ge_50", xmin=50.0, xmax=math.inf),
        Region(label="0_30", xmin=0.0, xmax=30.0),
        Region(label="30_60", xmin=30.0, xmax=60.0),
        Region(label="60_90", xmin=60.0, xmax=90.0),
        Region(label="90_120", xmin=90.0, xmax=120.0),
    ]


def format_console_summary(summary_df: pd.DataFrame, primary_threshold: float) -> str:
    focus = summary_df.loc[np.isclose(summary_df["threshold_um"], primary_threshold)].copy()
    focus = focus.sort_values(
        by="region",
        key=lambda series: series.map(
            {
                "whole_image": 0,
                "x_ge_50": 1,
                "0_30": 2,
                "30_60": 3,
                "60_90": 4,
                "90_120": 5,
            }
        ),
    )
    columns = [
        "region",
        "observed_pairs",
        "actin_total",
        "arc_total",
        "actin_coassociation_fraction",
        "arc_coassociation_fraction",
        "all_puncta_coassociation_fraction",
        "null_pairs_mean",
        "enrichment_vs_null",
        "empirical_p_enrichment",
    ]
    return focus[columns].to_string(index=False, float_format=lambda value: f"{value:0.4f}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze actin/arc puncta co-association by distance from the soma.")
    parser.add_argument("--actin", type=Path, required=True, help="Path to the actin puncta CSV.")
    parser.add_argument("--arc", type=Path, required=True, help="Path to the arc puncta CSV.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory where summary outputs should be written.",
    )
    parser.add_argument(
        "--simulations",
        type=int,
        default=1000,
        help="Number of null-model permutations to run.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility.",
    )
    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=list(DEFAULT_THRESHOLDS),
        help="Colocalization thresholds in micrometers.",
    )
    parser.add_argument(
        "--primary-threshold",
        type=float,
        default=0.5,
        help="Threshold used for the console summary.",
    )
    parser.add_argument(
        "--null-model",
        type=str,
        default="fix_x_permute_yz",
        choices=list(NULL_MODELS),
        help="Null model used to generate background co-association counts.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    regions = build_regions()
    thresholds = tuple(sorted(args.thresholds))
    images, actin_df, arc_df = load_tables(args.actin, args.arc)
    summary_df = summarize_counts(
        images=images,
        regions=regions,
        thresholds=thresholds,
        simulations=args.simulations,
        seed=args.seed,
        null_model=args.null_model,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "puncta_colocalization_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    metadata_path = args.output_dir / "puncta_colocalization_metadata.txt"
    metadata_path.write_text(
        "\n".join(
            [
                f"actin_rows={len(actin_df)}",
                f"arc_rows={len(arc_df)}",
                f"images={len(images)}",
                f"simulations={args.simulations}",
                f"seed={args.seed}",
                f"thresholds_um={','.join(str(value) for value in thresholds)}",
                f"null_model={args.null_model}",
                f"null_model_description={NULL_MODEL_DESCRIPTIONS[args.null_model]}",
                "matching=Maximum one-to-one bipartite matching under a 3D Euclidean distance cutoff.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote summary to {summary_path}")
    print(f"Wrote metadata to {metadata_path}")
    print()
    print(f"Primary threshold: {args.primary_threshold} um")
    print(format_console_summary(summary_df, args.primary_threshold))


if __name__ == "__main__":
    main()
