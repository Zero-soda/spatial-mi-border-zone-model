from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from sstba.cli import main


ROOT = Path(__file__).resolve().parents[1]


def test_real_demo_reproduces_released_graph_summaries(tmp_path: Path) -> None:
    output = tmp_path / "demo"
    assert (
        main(
            [
                "run",
                "--config",
                str(ROOT / "config" / "demo_boundary.json"),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    observed = pd.read_csv(output / "boundary_summary.tsv", sep="\t")
    expected = pd.read_csv(
        ROOT
        / "source_data"
        / "Source_Data_Supplementary_Figure_S6_visium_graph_boundary_analysis.tsv",
        sep="\t",
    )
    merged = observed.merge(expected, on="sample")

    assert np.array_equal(merged["n_edges"], merged["n_graph_edges"])
    assert np.array_equal(
        merged["n_cross_label_edges"],
        merged["domain3_domain4_edge_count"],
    )
    assert np.allclose(
        merged["boundary_fraction_3"],
        merged["domain3_graph_boundary_fraction"],
        atol=1e-8,
    )
    assert np.allclose(
        merged["boundary_fraction_4"],
        merged["domain4_graph_boundary_fraction"],
        atol=1e-8,
    )

    edge_expected = pd.read_csv(
        ROOT
        / "source_data"
        / "Source_Data_Supplementary_Figure_S6_visium_graph_boundary_edge_gradients.tsv",
        sep="\t",
    )
    for score, observed_column in {
        "mechanical_border": "edge_gradient_mechanical_border_score",
        "immune_fibrotic_activation": "edge_gradient_immune_fibrotic_activation_score",
        "fibroblast_scar_repair": "edge_gradient_fibroblast_scar_repair_score",
    }.items():
        score_expected = edge_expected.loc[
            edge_expected["score"] == score,
            ["sample", "mean_edge_delta"],
        ]
        comparison = observed[["sample", observed_column]].merge(
            score_expected,
            on="sample",
        )
        assert np.allclose(
            comparison[observed_column],
            comparison["mean_edge_delta"],
            atol=1e-7,
        )
