from __future__ import annotations

import pandas as pd
import pytest

from sstba.boundary import analyse_boundary


def test_boundary_marks_cross_label_neighbours() -> None:
    coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0, 3.0], "y": [0.0, 0.0, 0.0, 0.0]},
        index=["s1", "s2", "s3", "s4"],
    )
    labels = pd.Series(["3", "3", "4", "4"], index=coordinates.index)

    result = analyse_boundary(coordinates, labels, radius=1.1)

    assert not bool(result.spots.loc["s1", "is_boundary"])
    assert bool(result.spots.loc["s2", "is_boundary"])
    assert bool(result.spots.loc["s3", "is_boundary"])
    assert not bool(result.spots.loc["s4", "is_boundary"])
    assert result.summary["n_cross_label_edges"] == 1
    assert result.summary["boundary_fraction_3"] == 0.5
    assert result.summary["boundary_fraction_4"] == 0.5


def test_boundary_score_gradient_uses_label_order() -> None:
    coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0, 3.0], "y": [0.0, 0.0, 0.0, 0.0]},
        index=["s1", "s2", "s3", "s4"],
    )
    labels = pd.Series(["3", "3", "4", "4"], index=coordinates.index)
    scores = pd.DataFrame(
        {"scar": [0.0, 1.0, 3.0, 4.0], "mechanical": [4.0, 3.0, 1.0, 0.0]},
        index=coordinates.index,
    )

    result = analyse_boundary(
        coordinates,
        labels,
        radius=1.1,
        scores=scores,
        positive_label="4",
        negative_label="3",
    )

    assert result.summary["domain_gradient_scar"] == 3.0
    assert result.summary["domain_gradient_mechanical"] == -3.0
    assert result.summary["edge_gradient_scar"] == 2.0
    assert result.summary["edge_gradient_mechanical"] == -2.0


def test_boundary_requires_coordinate_alignment() -> None:
    coordinates = pd.DataFrame({"x": [0.0], "y": [0.0]}, index=["s1"])
    labels = pd.Series(["3"], index=["other"])

    with pytest.raises(ValueError, match="aligned"):
        analyse_boundary(coordinates, labels, radius=1.0)


def test_boundary_rejects_nonpositive_radius() -> None:
    coordinates = pd.DataFrame({"x": [0.0], "y": [0.0]}, index=["s1"])
    labels = pd.Series(["3"], index=["s1"])

    with pytest.raises(ValueError, match="positive"):
        analyse_boundary(coordinates, labels, radius=0.0)


def test_knn_boundary_matches_one_nearest_neighbour_graph() -> None:
    coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 3.0], "y": [0.0, 0.0, 0.0]},
        index=["s1", "s2", "s3"],
    )
    labels = pd.Series(["3", "3", "4"], index=coordinates.index)

    result = analyse_boundary(
        coordinates,
        labels,
        n_neighbors=1,
        knn_algorithm="exact",
    )

    assert (
        result.summary["graph_method"]
        == "symmetrized_1_nearest_neighbours_exact"
    )
    assert result.summary["n_edges"] == 2
    assert result.summary["n_cross_label_edges"] == 1


def test_target_boundary_ignores_other_label_interfaces() -> None:
    coordinates = pd.DataFrame(
        {"x": [0.0, 1.0, 2.0, 3.0], "y": [0.0, 0.0, 0.0, 0.0]},
        index=["s1", "s2", "s3", "s4"],
    )
    labels = pd.Series(["3", "3", "4", "5"], index=coordinates.index)
    scores = pd.DataFrame({"scar": [0.0, 1.0, 3.0, 9.0]}, index=coordinates.index)

    result = analyse_boundary(
        coordinates,
        labels,
        radius=1.1,
        scores=scores,
        positive_label="4",
        negative_label="3",
    )

    assert result.summary["n_cross_label_edges"] == 1
    assert result.spots["is_boundary"].tolist() == [False, True, True, False]
    assert result.summary["edge_gradient_scar"] == 2.0


def test_exact_knn_breaks_equal_distance_ties_by_spot_order() -> None:
    coordinates = pd.DataFrame(
        {"x": [0.0, -1.0, 1.0], "y": [0.0, 0.0, 0.0]},
        index=["centre", "left", "right"],
    )
    labels = pd.Series(["3", "4", "5"], index=coordinates.index)

    result = analyse_boundary(
        coordinates,
        labels,
        n_neighbors=1,
        knn_algorithm="exact",
    )

    centre_edges = result.edges[
        (result.edges["spot_left"] == "centre")
        | (result.edges["spot_right"] == "centre")
    ]
    assert "left" in set(centre_edges[["spot_left", "spot_right"]].to_numpy().ravel())
