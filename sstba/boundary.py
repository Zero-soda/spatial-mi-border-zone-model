"""Spatial-neighbour graph and cross-label boundary summaries."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class BoundaryResult:
    """Boundary labels, graph edges and sample-level summary."""

    spots: pd.DataFrame
    edges: pd.DataFrame
    summary: dict[str, object]


def _validate_inputs(
    coordinates: pd.DataFrame,
    labels: pd.Series,
    scores: pd.DataFrame | None,
) -> None:
    if len(coordinates) != len(labels) or not coordinates.index.equals(labels.index):
        raise ValueError("Coordinates and labels must be aligned by spot identifier")
    if coordinates.index.has_duplicates:
        raise ValueError("Duplicate spot identifiers are not allowed")
    if coordinates.shape[1] != 2:
        raise ValueError("Coordinates must contain exactly two columns")
    if coordinates.isna().any().any():
        raise ValueError("Coordinates contain missing values")
    if labels.isna().any():
        raise ValueError("Labels contain missing values")
    if scores is not None and not coordinates.index.equals(scores.index):
        raise ValueError("Coordinates and scores must be aligned by spot identifier")


def analyse_boundary(
    coordinates: pd.DataFrame,
    labels: pd.Series,
    *,
    radius: float | None = None,
    n_neighbors: int | None = None,
    knn_algorithm: str = "exact",
    scores: pd.DataFrame | None = None,
    positive_label: str | None = None,
    negative_label: str | None = None,
) -> BoundaryResult:
    """Construct a radius graph and quantify cross-label boundary structure."""

    if (radius is None) == (n_neighbors is None):
        raise ValueError("Supply exactly one of radius or n_neighbors")
    if radius is not None and (not np.isfinite(radius) or radius <= 0):
        raise ValueError("radius must be a positive finite number")
    if n_neighbors is not None and (
        not isinstance(n_neighbors, int) or n_neighbors < 1
    ):
        raise ValueError("n_neighbors must be a positive integer")
    if knn_algorithm not in {"exact", "kdtree"}:
        raise ValueError("knn_algorithm must be 'exact' or 'kdtree'")
    _validate_inputs(coordinates, labels, scores)

    coordinate_array = coordinates.astype(float).to_numpy()
    label_array = labels.astype(str).to_numpy()
    spot_ids = coordinates.index.astype(str).to_numpy()
    if radius is not None:
        edge_pairs = sorted(cKDTree(coordinate_array).query_pairs(radius))
        graph_method = f"radius_{radius:g}"
    else:
        if len(coordinate_array) < 2:
            edge_pairs = []
        else:
            k = min(int(n_neighbors), len(coordinate_array) - 1)
            edge_set: set[tuple[int, int]] = set()
            if knn_algorithm == "exact":
                for left in range(len(coordinate_array)):
                    deltas = coordinate_array - coordinate_array[left]
                    distances = np.sqrt(np.sum(deltas * deltas, axis=1))
                    # The discovery analysis used NumPy's quicksort ordering.
                    # Pinning the sort kind preserves its deterministic
                    # equal-distance behaviour on the Visium hexagonal grid.
                    order = np.argsort(distances, kind="quicksort")
                    neighbours = order[order != left][:k]
                    for right in neighbours:
                        right = int(right)
                        edge_set.add(
                            (left, right) if left < right else (right, left)
                        )
            else:
                _, neighbour_indices = cKDTree(coordinate_array).query(
                    coordinate_array,
                    k=k + 1,
                )
                if neighbour_indices.ndim == 1:
                    neighbour_indices = neighbour_indices[:, None]
                for left, neighbours in enumerate(neighbour_indices):
                    for right in neighbours[1:]:
                        right = int(right)
                        edge_set.add(
                            (left, right) if left < right else (right, left)
                        )
            edge_pairs = sorted(edge_set)
        graph_method = (
            f"symmetrized_{n_neighbors}_nearest_neighbours_{knn_algorithm}"
        )

    edge_rows: list[dict[str, object]] = []
    boundary_mask = np.zeros(len(coordinates), dtype=bool)
    degree = np.zeros(len(coordinates), dtype=int)
    cross_label_degree = np.zeros(len(coordinates), dtype=int)
    n_cross_label_edges = 0
    target_labels = (
        {str(positive_label), str(negative_label)}
        if positive_label is not None and negative_label is not None
        else None
    )
    target_edge_pairs: list[tuple[int, int]] = []

    for left, right in edge_pairs:
        is_any_cross_label = label_array[left] != label_array[right]
        is_cross_label = is_any_cross_label and (
            target_labels is None
            or {label_array[left], label_array[right]} == target_labels
        )
        degree[left] += 1
        degree[right] += 1
        if is_cross_label:
            boundary_mask[left] = True
            boundary_mask[right] = True
            cross_label_degree[left] += 1
            cross_label_degree[right] += 1
            n_cross_label_edges += 1
            target_edge_pairs.append((left, right))
        edge_rows.append(
            {
                "spot_left": spot_ids[left],
                "spot_right": spot_ids[right],
                "label_left": label_array[left],
                "label_right": label_array[right],
                "is_any_cross_label": bool(is_any_cross_label),
                "is_cross_label": bool(is_cross_label),
            }
        )

    spot_table = pd.DataFrame(
        {
            "label": label_array,
            "is_boundary": boundary_mask,
            "graph_degree": degree,
            "cross_label_degree": cross_label_degree,
        },
        index=coordinates.index,
    )
    edge_table = pd.DataFrame(
        edge_rows,
        columns=[
            "spot_left",
            "spot_right",
            "label_left",
            "label_right",
            "is_any_cross_label",
            "is_cross_label",
        ],
    )
    summary: dict[str, object] = {
        "n_spots": len(coordinates),
        "n_edges": len(edge_pairs),
        "n_cross_label_edges": n_cross_label_edges,
        "boundary_fraction": float(boundary_mask.mean()) if len(boundary_mask) else np.nan,
        "graph_method": graph_method,
    }
    for label in sorted(set(label_array)):
        label_mask = label_array == label
        summary[f"n_spots_{label}"] = int(label_mask.sum())
        summary[f"boundary_fraction_{label}"] = float(
            boundary_mask[label_mask].mean()
        )

    if scores is not None:
        numeric_scores = scores.apply(pd.to_numeric, errors="coerce")
        if positive_label is None or negative_label is None:
            raise ValueError(
                "positive_label and negative_label are required when scores are supplied"
            )
        positive_mask = label_array == str(positive_label)
        negative_mask = label_array == str(negative_label)
        if not positive_mask.any() or not negative_mask.any():
            raise ValueError("Both positive and negative labels must be present")
        for column in numeric_scores:
            summary[f"domain_gradient_{column}"] = float(
                numeric_scores.loc[positive_mask, column].mean()
                - numeric_scores.loc[negative_mask, column].mean()
            )
            edge_differences: list[float] = []
            values = numeric_scores[column].to_numpy(dtype=float)
            for left, right in target_edge_pairs:
                if label_array[left] == str(positive_label):
                    positive_index, negative_index = left, right
                else:
                    positive_index, negative_index = right, left
                edge_differences.append(
                    float(values[positive_index] - values[negative_index])
                )
            summary[f"edge_gradient_{column}"] = (
                float(np.nanmean(edge_differences))
                if edge_differences
                else float("nan")
            )
            spot_table[column] = numeric_scores[column]

    return BoundaryResult(spots=spot_table, edges=edge_table, summary=summary)
