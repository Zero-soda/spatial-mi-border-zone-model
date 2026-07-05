#!/usr/bin/env python3
"""Batch-map D3/D7 GSE214611 MI Visium prototype fibrotic-risk scores."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import tarfile
import zipfile
from dataclasses import dataclass
from pathlib import Path

from map_gse214611_d7_1_spatial_risk import ROOT, run_mapping, sample_slug


RAW_VISIUM_DIR = ROOT / "data" / "raw" / "gse214611" / "visium"
OUT_TABLE_DIR = ROOT / "results" / "tables"


@dataclass(frozen=True)
class SampleConfig:
    sample: str
    gsm: str
    archive_name: str

    @property
    def url(self) -> str:
        return (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6613nnn/"
            f"{self.gsm}/suppl/{self.archive_name}"
        )

    @property
    def archive_path(self) -> Path:
        return RAW_VISIUM_DIR / self.archive_name

    @property
    def unpack_dir(self) -> Path:
        stem = self.archive_name
        for suffix in (".tar.gz", ".zip"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
        return RAW_VISIUM_DIR / stem


SAMPLES: dict[str, SampleConfig] = {
    "D3_1": SampleConfig("D3_1", "GSM6613084", "GSM6613084_V_d1_1.tar.gz"),
    "D3_2": SampleConfig("D3_2", "GSM6613085", "GSM6613085_V_d1_2.tar.gz"),
    "D3_3": SampleConfig("D3_3", "GSM6613086", "GSM6613086_V_d1_3.tar.gz"),
    "D7_1": SampleConfig("D7_1", "GSM6613087", "GSM6613087_V_d7_1.zip"),
    "D7_2": SampleConfig("D7_2", "GSM6613088", "GSM6613088_V_d7_2.tar.gz"),
    "D7_3": SampleConfig("D7_3", "GSM6613089", "GSM6613089_V_d7_3.tar.gz"),
}


def download_archive(config: SampleConfig, force: bool = False) -> None:
    RAW_VISIUM_DIR.mkdir(parents=True, exist_ok=True)
    if config.archive_path.exists() and config.archive_path.stat().st_size > 1024 and not force:
        print(f"[skip-download] {config.archive_name}")
        return

    print(f"[download] {config.sample}: {config.url}")
    subprocess.run(
        [
            "curl",
            "-L",
            "--fail",
            "--retry",
            "3",
            "--retry-delay",
            "5",
            "-o",
            str(config.archive_path),
            config.url,
        ],
        check=True,
    )


def verify_archive(config: SampleConfig) -> None:
    print(f"[verify] {config.archive_name}")
    if config.archive_name.endswith(".zip"):
        with zipfile.ZipFile(config.archive_path) as archive:
            bad_file = archive.testzip()
        if bad_file is not None:
            raise RuntimeError(f"ZIP integrity check failed at {bad_file}")
        return

    if config.archive_name.endswith(".tar.gz"):
        with tarfile.open(config.archive_path, "r:gz") as archive:
            for member in archive:
                if member.name.startswith("/") or ".." in Path(member.name).parts:
                    raise RuntimeError(f"Unsafe archive member path: {member.name}")
        return

    raise ValueError(f"Unsupported archive type: {config.archive_name}")


def expected_visium_dir(path: Path) -> bool:
    standard_spatial = (
        (path / "spatial" / "tissue_positions_list.csv").exists()
        and (path / "spatial" / "tissue_lowres_image.png").exists()
    )
    root_spatial = (
        (path / "tissue_positions_list.csv").exists()
        and (path / "tissue_lowres_image.png").exists()
    )
    return standard_spatial or root_spatial


def unpack_archive(config: SampleConfig, force: bool = False) -> Path:
    if config.unpack_dir.exists() and not force:
        try:
            visium_dir = locate_visium_dir(config.unpack_dir)
            print(f"[skip-unpack] {visium_dir}")
            return visium_dir
        except RuntimeError:
            pass

    if config.unpack_dir.exists():
        shutil.rmtree(config.unpack_dir)
    config.unpack_dir.mkdir(parents=True, exist_ok=True)

    print(f"[unpack] {config.archive_name}")
    if config.archive_name.endswith(".zip"):
        with zipfile.ZipFile(config.archive_path) as archive:
            archive.extractall(config.unpack_dir)
    elif config.archive_name.endswith(".tar.gz"):
        with tarfile.open(config.archive_path, "r:gz") as archive:
            try:
                archive.extractall(config.unpack_dir, filter="data")
            except TypeError:
                archive.extractall(config.unpack_dir)
    else:
        raise ValueError(f"Unsupported archive type: {config.archive_name}")

    return locate_visium_dir(config.unpack_dir)


def locate_visium_dir(unpack_dir: Path) -> Path:
    if expected_visium_dir(unpack_dir):
        return unpack_dir

    candidates: list[Path] = []
    for match in unpack_dir.glob("**/tissue_positions_list.csv"):
        candidate = match.parents[1] if match.parent.name == "spatial" else match.parent
        if expected_visium_dir(candidate):
            candidates.append(candidate)

    unique_candidates = sorted(set(candidates))
    if len(unique_candidates) == 1:
        return unique_candidates[0]
    if not unique_candidates:
        raise RuntimeError(f"No Visium spatial directory found under {unpack_dir}")
    raise RuntimeError(
        f"Multiple Visium spatial directories found under {unpack_dir}: {unique_candidates}"
    )


def write_batch_summary(results: list[dict[str, str | int]]) -> Path:
    OUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_TABLE_DIR / "gse214611_d3_d7_batch_spatial_mapping_summary.tsv"
    fieldnames = [
        "sample",
        "n_spots",
        "domains",
        "visium_dir",
        "spatial_table",
        "domain_summary",
        "top_risk",
        "figure",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows({field: result[field] for field in fieldnames} for result in results)
    return out_path


def sample_stage(sample: str) -> str:
    if sample.startswith("D3"):
        return "day3_mi"
    if sample.startswith("D7"):
        return "day7_mi"
    return "other"


def write_combined_domain_summary(results: list[dict[str, str | int]]) -> Path:
    OUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_TABLE_DIR / "gse214611_d3_d7_batch_domain_summary.tsv"
    combined_rows: list[dict[str, str]] = []

    for result in results:
        sample = str(result["sample"])
        with Path(str(result["domain_summary"])).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                combined_rows.append({"sample": sample, "stage": sample_stage(sample), **row})

    if not combined_rows:
        raise ValueError("No domain summaries to combine")

    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(combined_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(combined_rows)

    return out_path


def write_stage_domain_summary(combined_domain_summary: Path) -> Path:
    out_path = OUT_TABLE_DIR / "gse214611_d3_d7_batch_stage_domain_summary.tsv"
    score_columns = [
        "mean_prototype_fibrotic_risk",
        "mean_prototype_repair_proxy",
        "mean_prototype_immune_proxy",
        "mean_prototype_border_activation",
        "mean_BZ1_genes",
        "mean_BZ2_genes",
        "mean_RZ_genes",
        "mean_X33_Fib_Postn",
        "mean_X36_Fib_Rep",
        "mean_CM_Score",
    ]
    sums: dict[tuple[str, str], dict[str, float]] = {}
    counts: dict[tuple[str, str], int] = {}

    with combined_domain_summary.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["stage"], row["annotated"])
            sums.setdefault(key, {col: 0.0 for col in score_columns})
            counts[key] = counts.get(key, 0) + 1
            for col in score_columns:
                sums[key][col] += float(row[col])

    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["stage", "annotated", "n_replicates", *score_columns])
        for key in sorted(counts):
            values = [f"{sums[key][col] / counts[key]:.6g}" for col in score_columns]
            writer.writerow([*key, counts[key], *values])

    return out_path


def write_contact_sheet(results: list[dict[str, str | int]]) -> Path:
    from PIL import Image, ImageDraw

    out_path = ROOT / "results" / "figures" / "gse214611_d3_d7_spatial_risk_contact_sheet.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    resample = getattr(Image, "Resampling", Image).LANCZOS

    panels: list[tuple[str, Image.Image]] = []
    for result in results:
        sample = str(result["sample"])
        image = Image.open(Path(str(result["figure"]))).convert("RGB")
        target_width = 1100
        target_height = round(image.height * target_width / image.width)
        panels.append((sample, image.resize((target_width, target_height), resample)))

    columns = 2
    rows = (len(panels) + columns - 1) // columns
    margin = 36
    gap = 30
    title_height = 40
    panel_width = max(image.width for _, image in panels)
    panel_height = max(image.height for _, image in panels)
    canvas = Image.new(
        "RGB",
        (
            columns * panel_width + (columns - 1) * gap + 2 * margin,
            rows * (panel_height + title_height) + (rows - 1) * gap + 2 * margin,
        ),
        "white",
    )
    draw = ImageDraw.Draw(canvas)

    for idx, (sample, image) in enumerate(panels):
        row = idx // columns
        col = idx % columns
        x = margin + col * (panel_width + gap)
        y = margin + row * (panel_height + title_height + gap)
        draw.text((x, y), sample, fill=(25, 25, 25))
        canvas.paste(image, (x, y + title_height))

    canvas.save(out_path, dpi=(600, 600))
    canvas.save(out_path.with_suffix(".tiff"), dpi=(600, 600), compression="tiff_lzw")
    return out_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--samples",
        nargs="+",
        default=list(SAMPLES),
        choices=sorted(SAMPLES),
        help="GSE214611 orig.ident values to process.",
    )
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--force-unpack", action="store_true")
    parser.add_argument("--min-spots", type=int, default=2000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results: list[dict[str, str | int]] = []

    for sample in args.samples:
        config = SAMPLES[sample]
        download_archive(config, force=args.force_download)
        verify_archive(config)
        visium_dir = unpack_archive(config, force=args.force_unpack)

        if args.download_only:
            continue

        print(f"[map] {sample}: {visium_dir}")
        result = run_mapping(sample, visium_dir, min_spots=args.min_spots)
        results.append(result)
        print(
            f"[done] {sample}: {result['n_spots']} spots, "
            f"domains {result['domains']}, figure {Path(str(result['figure'])).name}"
        )

    if results:
        summary_path = write_batch_summary(results)
        domain_summary_path = write_combined_domain_summary(results)
        stage_domain_summary_path = write_stage_domain_summary(domain_summary_path)
        contact_sheet_path = write_contact_sheet(results)
        print(f"[summary] {summary_path}")
        print(f"[domain-summary] {domain_summary_path}")
        print(f"[stage-domain-summary] {stage_domain_summary_path}")
        print(f"[contact-sheet] {contact_sheet_path}")


if __name__ == "__main__":
    main()
