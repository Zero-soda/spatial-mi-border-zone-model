#!/usr/bin/env python3
"""Shared helpers for spatial MI border-zone revision analyses."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def cohens_d(left: list[float] | np.ndarray, right: list[float] | np.ndarray) -> float:
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    pooled = math.sqrt(((len(x) - 1) * np.var(x, ddof=1) + (len(y) - 1) * np.var(y, ddof=1)) / (len(x) + len(y) - 2))
    if pooled == 0:
        return float("nan")
    return float((np.mean(x) - np.mean(y)) / pooled)


def bootstrap_ci(values: list[float] | np.ndarray, n_boot: int = 2000, seed: int = 614611) -> tuple[float, float, float]:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    samples = rng.choice(arr, size=(n_boot, len(arr)), replace=True)
    means = np.mean(samples, axis=1)
    return float(np.mean(arr)), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def pairwise_distances(coords: np.ndarray) -> np.ndarray:
    delta = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def morans_i(coords: np.ndarray, values: np.ndarray, radius: float) -> float:
    coords = np.asarray(coords, dtype=float)
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values) & np.all(np.isfinite(coords), axis=1)
    coords = coords[valid]
    values = values[valid]
    n = len(values)
    if n < 3:
        return float("nan")
    distances = pairwise_distances(coords)
    weights = (distances > 0) & (distances <= radius)
    w_sum = float(np.sum(weights))
    if w_sum == 0:
        return float("nan")
    centered = values - np.mean(values)
    denom = float(np.sum(centered * centered))
    if denom == 0:
        return float("nan")
    num = float(np.sum(weights * (centered[:, None] * centered[None, :])))
    return float((n / w_sum) * (num / denom))


def gearys_c(coords: np.ndarray, values: np.ndarray, radius: float) -> float:
    coords = np.asarray(coords, dtype=float)
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values) & np.all(np.isfinite(coords), axis=1)
    coords = coords[valid]
    values = values[valid]
    n = len(values)
    if n < 3:
        return float("nan")
    distances = pairwise_distances(coords)
    weights = (distances > 0) & (distances <= radius)
    w_sum = float(np.sum(weights))
    if w_sum == 0:
        return float("nan")
    denom = float(np.sum((values - np.mean(values)) ** 2))
    if denom == 0:
        return float("nan")
    diff_sq = (values[:, None] - values[None, :]) ** 2
    return float(((n - 1) / (2 * w_sum)) * (np.sum(weights * diff_sq) / denom))


def domain_label_permutation_margin(
    values: np.ndarray,
    labels: np.ndarray,
    positive_label: str,
    negative_label: str,
    n_perm: int = 5000,
    seed: int = 614611,
) -> tuple[float, float, np.ndarray]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels).astype(str)
    mask = np.isfinite(values) & np.isin(labels, [positive_label, negative_label])
    values = values[mask]
    labels = labels[mask]
    positive_mask = labels == positive_label
    negative_mask = labels == negative_label
    if not np.any(positive_mask) or not np.any(negative_mask):
        return float("nan"), float("nan"), np.array([], dtype=float)
    observed = float(np.mean(values[positive_mask]) - np.mean(values[negative_mask]))
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    for idx in range(n_perm):
        shuffled = rng.permutation(labels)
        null[idx] = np.mean(values[shuffled == positive_label]) - np.mean(values[shuffled == negative_label])
    p_upper = float((np.sum(null >= observed) + 1) / (n_perm + 1))
    return observed, p_upper, null


def spatial_block_permutation_margin(
    values: np.ndarray,
    labels: np.ndarray,
    coords: np.ndarray,
    positive_label: str,
    negative_label: str,
    n_perm: int = 5000,
    seed: int = 614611,
    grid_size: int = 4,
) -> tuple[float, float, np.ndarray, int]:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels).astype(str)
    coords = np.asarray(coords, dtype=float)
    mask = np.isfinite(values) & np.isin(labels, [positive_label, negative_label]) & np.all(np.isfinite(coords), axis=1)
    values = values[mask]
    labels = labels[mask]
    coords = coords[mask]
    observed = float(np.mean(values[labels == positive_label]) - np.mean(values[labels == negative_label]))
    if len(values) < 4:
        return observed, float("nan"), np.array([], dtype=float), 0
    x_edges = np.linspace(float(np.min(coords[:, 0])), float(np.max(coords[:, 0])), grid_size + 1)
    y_edges = np.linspace(float(np.min(coords[:, 1])), float(np.max(coords[:, 1])), grid_size + 1)
    x_bin = np.clip(np.digitize(coords[:, 0], x_edges[1:-1]), 0, grid_size - 1)
    y_bin = np.clip(np.digitize(coords[:, 1], y_edges[1:-1]), 0, grid_size - 1)
    block_ids = x_bin + grid_size * y_bin
    unique_blocks = np.unique(block_ids)
    block_labels = []
    for block_id in unique_blocks:
        block_domain_labels = labels[block_ids == block_id]
        n_pos = np.sum(block_domain_labels == positive_label)
        n_neg = np.sum(block_domain_labels == negative_label)
        block_labels.append(positive_label if n_pos >= n_neg else negative_label)
    block_labels = np.asarray(block_labels)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm, dtype=float)
    for idx in range(n_perm):
        shuffled_block_labels = rng.permutation(block_labels)
        shuffled_labels = np.empty_like(labels)
        for block_idx, block_id in enumerate(unique_blocks):
            shuffled_labels[block_ids == block_id] = shuffled_block_labels[block_idx]
        if np.any(shuffled_labels == positive_label) and np.any(shuffled_labels == negative_label):
            null[idx] = np.mean(values[shuffled_labels == positive_label]) - np.mean(values[shuffled_labels == negative_label])
        else:
            null[idx] = np.nan
    null = null[np.isfinite(null)]
    if len(null) == 0:
        return observed, float("nan"), null, len(unique_blocks)
    p_upper = float((np.sum(null >= observed) + 1) / (len(null) + 1))
    return observed, p_upper, null, len(unique_blocks)


def signed_distance_slope(distance: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    x = np.asarray(distance, dtype=float)
    y = np.asarray(score, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 3:
        return float("nan"), float("nan")
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept)
