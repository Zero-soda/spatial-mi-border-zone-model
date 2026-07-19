"""Command-line interface for the SSTBA workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

from .boundary import analyse_boundary
from .manifest import base_manifest, sha256_file
from .scoring import score_modules


BOUNDARY_REQUIRED_FIELDS = (
    "input",
    "sample_column",
    "spot_column",
    "label_column",
    "x_column",
    "y_column",
    "score_columns",
    "positive_label",
    "negative_label",
)


def _load_config(path: Path) -> dict[str, Any]:
    try:
        config = json.loads(path.read_text())
    except FileNotFoundError as error:
        raise ValueError(f"Configuration file does not exist: {path}") from error
    except json.JSONDecodeError as error:
        raise ValueError(f"Configuration is not valid JSON: {error}") from error
    if not isinstance(config, dict):
        raise ValueError("Configuration root must be a JSON object")
    return config


def _resolve_input(config_path: Path, raw_path: str) -> Path:
    input_path = Path(raw_path).expanduser()
    if not input_path.is_absolute():
        input_path = config_path.parent / input_path
    return input_path.resolve()


def _validate_boundary_config(config: dict[str, Any], config_path: Path) -> None:
    missing = [field for field in BOUNDARY_REQUIRED_FIELDS if field not in config]
    if missing:
        raise ValueError(f"Missing required configuration field(s): {', '.join(missing)}")
    if config.get("mode", "boundary") != "boundary":
        raise ValueError("Only mode 'boundary' is currently supported by run")
    input_path = _resolve_input(config_path, str(config["input"]))
    if not input_path.is_file():
        raise ValueError(f"Input file does not exist: {input_path}")
    named_fields = (
        "sample_column",
        "spot_column",
        "label_column",
        "x_column",
        "y_column",
    )
    for field in named_fields:
        if not isinstance(config[field], str) or not config[field].strip():
            raise ValueError(f"{field} must be a non-empty string")
    score_columns = config["score_columns"]
    if (
        not isinstance(score_columns, list)
        or not score_columns
        or any(not isinstance(column, str) or not column.strip() for column in score_columns)
    ):
        raise ValueError("score_columns must be a non-empty list of column names")
    if len(score_columns) != len(set(score_columns)):
        raise ValueError("score_columns contains duplicate column names")
    summary_columns = config.get("summary_columns", [])
    if (
        not isinstance(summary_columns, list)
        or any(
            not isinstance(column, str) or not column.strip()
            for column in summary_columns
        )
    ):
        raise ValueError("summary_columns must be a list of column names")
    if len(summary_columns) != len(set(summary_columns)):
        raise ValueError("summary_columns contains duplicate column names")
    if str(config["positive_label"]) == str(config["negative_label"]):
        raise ValueError("positive_label and negative_label must be different")
    algorithm = str(config.get("knn_algorithm", "exact"))
    if algorithm not in {"exact", "kdtree"}:
        raise ValueError("knn_algorithm must be 'exact' or 'kdtree'")
    has_radius = "radius" in config
    has_neighbours = "n_neighbors" in config
    if has_radius == has_neighbours:
        raise ValueError("Supply exactly one of radius or n_neighbors")
    if has_radius:
        try:
            radius = float(config["radius"])
        except (TypeError, ValueError) as error:
            raise ValueError("radius must be numeric") from error
        if radius <= 0:
            raise ValueError("radius must be positive")
    if has_neighbours:
        if (
            not isinstance(config["n_neighbors"], int)
            or config["n_neighbors"] < 1
        ):
            raise ValueError("n_neighbors must be a positive integer")
    seed = config.get("seed", 614611)
    if not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer")

    header = pd.read_csv(input_path, sep="\t", nrows=0)
    required_columns = {
        *(str(config[field]) for field in named_fields),
        *map(str, score_columns),
        *map(str, summary_columns),
    }
    missing_columns = sorted(required_columns - set(header.columns))
    if missing_columns:
        raise ValueError(f"Input table is missing column(s): {', '.join(missing_columns)}")


def _write_tsv(path: Path, table: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, lineterminator="\n")


def _run_boundary(
    config: dict[str, Any],
    config_path: Path,
    output: Path,
    *,
    command: str = "sstba run",
) -> dict[str, object]:
    _validate_boundary_config(config, config_path)
    input_path = _resolve_input(config_path, str(config["input"]))
    table = pd.read_csv(input_path, sep="\t", dtype={str(config["label_column"]): str})
    required_columns = {
        str(config["sample_column"]),
        str(config["spot_column"]),
        str(config["label_column"]),
        str(config["x_column"]),
        str(config["y_column"]),
        *map(str, config["score_columns"]),
        *map(str, config.get("summary_columns", [])),
    }
    missing_columns = sorted(required_columns - set(table.columns))
    if missing_columns:
        raise ValueError(f"Input table is missing column(s): {', '.join(missing_columns)}")

    summary_rows: list[dict[str, object]] = []
    spot_tables: list[pd.DataFrame] = []
    sample_column = str(config["sample_column"])
    spot_column = str(config["spot_column"])
    label_column = str(config["label_column"])
    score_columns = list(map(str, config["score_columns"]))
    summary_columns = list(map(str, config.get("summary_columns", [])))

    for sample, sample_table in table.groupby(sample_column, sort=True):
        sample_table = sample_table.copy()
        sample_table.index = sample_table[spot_column].astype(str)
        result = analyse_boundary(
            sample_table[[str(config["x_column"]), str(config["y_column"])]],
            sample_table[label_column].astype(str),
            radius=float(config["radius"]) if "radius" in config else None,
            n_neighbors=(
                int(config["n_neighbors"])
                if "n_neighbors" in config
                else None
            ),
            knn_algorithm=str(config.get("knn_algorithm", "exact")),
            scores=sample_table[score_columns],
            positive_label=str(config["positive_label"]),
            negative_label=str(config["negative_label"]),
        )
        summary_metadata: dict[str, object] = {}
        for column in summary_columns:
            values = sample_table[column].dropna().unique()
            if len(values) > 1:
                raise ValueError(
                    f"Summary column '{column}' varies within sample '{sample}'"
                )
            summary_metadata[column] = values[0] if len(values) else ""
        summary_rows.append(
            {sample_column: sample, **summary_metadata, **result.summary}
        )
        spots = result.spots.reset_index(names=spot_column)
        spots.insert(0, sample_column, sample)
        spot_tables.append(spots)

    summary_table = pd.DataFrame(summary_rows)
    spots_table = pd.concat(spot_tables, ignore_index=True)
    output.mkdir(parents=True, exist_ok=True)
    summary_path = output / "boundary_summary.tsv"
    spots_path = output / "boundary_spots.tsv"
    _write_tsv(summary_path, summary_table)
    _write_tsv(spots_path, spots_table)

    manifest = base_manifest(
        mode="boundary",
        output_directory=output,
        seed=int(config.get("seed", 614611)),
    )
    manifest.update(
        {
            "command": command,
            "config_sha256": sha256_file(config_path),
            "parameters": {
                key: config[key]
                for key in (
                    "sample_column",
                    "spot_column",
                    "label_column",
                    "x_column",
                    "y_column",
                    "score_columns",
                    "positive_label",
                    "negative_label",
                )
            }
            | {
                key: config[key]
                for key in (
                    "radius",
                    "n_neighbors",
                    "knn_algorithm",
                    "summary_columns",
                )
                if key in config
            },
            "inputs": {
                "spot_table": {
                    "path": str(config["input"]),
                    "sha256": sha256_file(input_path),
                }
            },
            "outputs": {
                summary_path.name: sha256_file(summary_path),
                spots_path.name: sha256_file(spots_path),
            },
        }
    )
    manifest_path = output / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def _read_signatures(path: Path) -> dict[str, list[str]]:
    table = pd.read_csv(path, sep="\t")
    if not {"module", "gene"}.issubset(table.columns):
        raise ValueError("Signature table must contain module and gene columns")
    return {
        str(module): group["gene"].astype(str).tolist()
        for module, group in table.groupby("module", sort=False)
    }


def _score_command(args: argparse.Namespace) -> int:
    expression_path = Path(args.expression).expanduser().resolve()
    signatures_path = Path(args.signatures).expanduser().resolve()
    expression = pd.read_csv(expression_path, sep="\t", index_col=0)
    signatures = _read_signatures(signatures_path)
    scores, coverage = score_modules(
        expression,
        signatures,
        min_detected_genes=args.min_detected_genes,
    )
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    scores_path = output / "scores.tsv"
    coverage_path = output / "signature_coverage.json"
    scores.to_csv(scores_path, sep="\t", lineterminator="\n")
    coverage_path.write_text(
        json.dumps(coverage, indent=2, sort_keys=True) + "\n"
    )
    manifest = base_manifest(
        mode="score",
        output_directory=output,
        seed=614611,
    )
    manifest.update(
        {
            "command": "sstba score",
            "parameters": {
                "min_detected_genes": int(args.min_detected_genes),
                "n_spots": int(expression.shape[0]),
                "n_expression_genes": int(expression.shape[1]),
            },
            "inputs": {
                "expression": {
                    "path": expression_path.name,
                    "sha256": sha256_file(expression_path),
                },
                "signatures": {
                    "path": signatures_path.name,
                    "sha256": sha256_file(signatures_path),
                },
            },
            "signature_coverage": coverage,
            "outputs": {
                scores_path.name: sha256_file(scores_path),
                coverage_path.name: sha256_file(coverage_path),
            },
        }
    )
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sstba",
        description="Frozen spatial-state scoring and boundary analysis",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate", help="validate a run configuration")
    validate.add_argument("--config", required=True)

    run = subparsers.add_parser("run", help="execute a configured workflow")
    run.add_argument("--config", required=True)
    run.add_argument("--output", required=True)

    boundary = subparsers.add_parser("boundary", help="execute boundary analysis")
    boundary.add_argument("--config", required=True)
    boundary.add_argument("--output", required=True)

    score = subparsers.add_parser("score", help="score frozen gene modules")
    score.add_argument("--expression", required=True)
    score.add_argument("--signatures", required=True)
    score.add_argument("--output", required=True)
    score.add_argument("--min-detected-genes", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "score":
            return _score_command(args)
        config_path = Path(args.config).expanduser().resolve()
        config = _load_config(config_path)
        _validate_boundary_config(config, config_path)
        if args.command == "validate":
            print(f"Valid configuration: {config_path}")
            return 0
        _run_boundary(
            config,
            config_path,
            Path(args.output),
            command=f"sstba {args.command}",
        )
        return 0
    except (OSError, ValueError) as error:
        print(f"sstba: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
