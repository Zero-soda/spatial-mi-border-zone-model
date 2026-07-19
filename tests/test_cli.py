from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from sstba.cli import main
from sstba.manifest import sha256_file


def write_demo_inputs(root: Path) -> Path:
    spots = pd.DataFrame(
        {
            "sample": ["A", "A", "A", "A"],
            "stage": ["day3", "day3", "day3", "day3"],
            "barcode": ["s1", "s2", "s3", "s4"],
            "domain": ["3", "3", "4", "4"],
            "x": [0.0, 1.0, 2.0, 3.0],
            "y": [0.0, 0.0, 0.0, 0.0],
            "mechanical": [4.0, 3.0, 1.0, 0.0],
            "scar": [0.0, 1.0, 3.0, 4.0],
        }
    )
    spots_path = root / "spots.tsv"
    spots.to_csv(spots_path, sep="\t", index=False)
    config = {
        "mode": "boundary",
        "input": str(spots_path),
        "sample_column": "sample",
        "spot_column": "barcode",
        "label_column": "domain",
        "x_column": "x",
        "y_column": "y",
        "score_columns": ["mechanical", "scar"],
        "summary_columns": ["stage"],
        "positive_label": "4",
        "negative_label": "3",
        "radius": 1.1,
        "seed": 614611,
    }
    config_path = root / "config.json"
    config_path.write_text(json.dumps(config, indent=2))
    return config_path


def test_sha256_is_content_stable(tmp_path: Path) -> None:
    path = tmp_path / "input.tsv"
    path.write_text("a\tb\n1\t2\n")

    assert sha256_file(path) == sha256_file(path)
    assert len(sha256_file(path)) == 64


def test_demo_cli_is_deterministic(tmp_path: Path) -> None:
    config = write_demo_inputs(tmp_path)
    first = tmp_path / "first"
    second = tmp_path / "second"

    assert main(["run", "--config", str(config), "--output", str(first)]) == 0
    assert main(["run", "--config", str(config), "--output", str(second)]) == 0

    assert (first / "boundary_summary.tsv").read_bytes() == (
        second / "boundary_summary.tsv"
    ).read_bytes()
    summary = pd.read_csv(first / "boundary_summary.tsv", sep="\t")
    assert summary.loc[0, "stage"] == "day3"
    manifest_first = json.loads((first / "run_manifest.json").read_text())
    manifest_second = json.loads((second / "run_manifest.json").read_text())
    for manifest in (manifest_first, manifest_second):
        manifest.pop("created_utc")
    assert manifest_first == manifest_second
    assert manifest_first["output_directory"] == "."
    assert manifest_first["command"] == "sstba run"
    assert len(manifest_first["config_sha256"]) == 64


def test_validate_reports_missing_required_field(tmp_path: Path, capsys) -> None:
    config_path = tmp_path / "invalid.json"
    config_path.write_text(json.dumps({"mode": "boundary"}))

    assert main(["validate", "--config", str(config_path)]) == 2
    assert "input" in capsys.readouterr().err


def test_validate_rejects_invalid_algorithm_equal_labels_and_missing_columns(
    tmp_path: Path,
    capsys,
) -> None:
    config_path = write_demo_inputs(tmp_path)
    config = json.loads(config_path.read_text())
    config["knn_algorithm"] = "approximate"
    config["positive_label"] = "3"
    config["negative_label"] = "3"
    config["score_columns"] = ["scar", "scar", "missing"]
    config_path.write_text(json.dumps(config))

    assert main(["validate", "--config", str(config_path)]) == 2
    message = capsys.readouterr().err
    assert (
        "knn_algorithm" in message
        or "positive_label" in message
        or "duplicate" in message
        or "missing" in message
    )


def test_score_writes_coverage_and_manifest(tmp_path: Path) -> None:
    expression = pd.DataFrame(
        {"POSTN": [1.0, 3.0], "CTHRC1": [0.0, 4.0]},
        index=["s1", "s2"],
    )
    expression_path = tmp_path / "expression.tsv"
    expression.to_csv(expression_path, sep="\t")
    signatures_path = tmp_path / "signatures.tsv"
    signatures_path.write_text("module\tgene\nscar\tPOSTN\nscar\tCTHRC1\n")
    output = tmp_path / "score-output"

    assert (
        main(
            [
                "score",
                "--expression",
                str(expression_path),
                "--signatures",
                str(signatures_path),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["command"] == "sstba score"
    assert manifest["signature_coverage"]["scar"]["detected"] == 2
    assert set(manifest["outputs"]) == {"scores.tsv", "signature_coverage.json"}
