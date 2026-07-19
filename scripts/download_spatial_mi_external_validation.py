#!/usr/bin/env python3
"""Download and checksum public external spatial MI validation files."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import time
from contextlib import contextmanager
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from project_paths import project_root
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from external_spatial_validation_utils import read_manifest


ROOT = project_root(__file__)
DEFAULT_MANIFEST = (
    Path(__file__).resolve().parents[1]
    / "config"
    / "external_validation_manifest.tsv"
)
DEFAULT_DATA_ROOT = ROOT / "data" / "raw" / "external_validation"
DEFAULT_LOG = ROOT / "results" / "tables" / "external_validation" / "retrieval_manifest.tsv"


def local_filename(row: dict[str, str]) -> str:
    if row["file_role"] == "h5ad":
        return f"{row['sample']}.h5ad"
    name = Path(urlparse(row["url"]).path).name
    if not name or name == "content":
        raise ValueError(f"cannot derive local file name from {row['url']}")
    return name


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_to_part(
    url: str,
    part: Path,
    expected_size: int | None,
    attempts: int = 3,
    opener=urlopen,
    sleep=time.sleep,
) -> None:
    """Download with range-based retries while preserving verified byte progress."""
    part.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, attempts + 1):
        offset = part.stat().st_size if part.exists() else 0
        if expected_size is not None and offset > expected_size:
            part.unlink()
            offset = 0
        headers = {"User-Agent": "SpatialMIValidation/0.4"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = Request(url, headers=headers)
        try:
            with opener(request, timeout=120) as response:
                status = getattr(response, "status", None)
                append = bool(offset and status == 206)
                mode = "ab" if append else "wb"
                with part.open(mode) as output:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        output.write(chunk)
            if expected_size is not None and part.stat().st_size != expected_size:
                raise IOError(
                    f"size mismatch for {part.name}: {part.stat().st_size} != {expected_size}"
                )
            return
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            if attempt == attempts:
                raise RuntimeError(f"failed to download {url}: {exc}") from exc
            sleep(2**attempt)


@contextmanager
def _exclusive_download_lock(destination: Path):
    lock = destination.with_suffix(destination.suffix + ".download.lock")
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError as exc:
        raise RuntimeError(f"download already active for {destination}") from exc
    try:
        yield
    finally:
        os.close(descriptor)
        lock.unlink(missing_ok=True)


def download_one(row: dict[str, str], data_root: Path, attempts: int = 3) -> dict[str, object]:
    destination = data_root / row["dataset"] / row["sample"] / local_filename(row)
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected_size = int(row["expected_size"]) if row.get("expected_size", "").strip() else None
    status = "downloaded"
    if destination.exists() and (expected_size is None or destination.stat().st_size == expected_size):
        status = "cached"
    else:
        part = destination.with_suffix(destination.suffix + ".part")
        with _exclusive_download_lock(destination):
            legacy_parts = sorted(
                destination.parent.glob(destination.name + ".part.*"),
                key=lambda path: path.stat().st_size,
                reverse=True,
            )
            if not part.exists() and legacy_parts:
                os.replace(legacy_parts[0], part)
            for stale in legacy_parts[1:]:
                stale.unlink(missing_ok=True)
            _download_to_part(row["url"], part, expected_size, attempts=attempts)
            os.replace(part, destination)
    return {
        **row,
        "local_path": str(destination.relative_to(ROOT)),
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": status,
    }


def write_log(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dataset",
        "study_role",
        "species",
        "sample",
        "patient",
        "stage_region",
        "file_role",
        "url",
        "expected_size",
        "source_citation",
        "local_path",
        "bytes",
        "sha256",
        "retrieved_at_utc",
        "status",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (str(row["dataset"]), str(row["sample"]), str(row["file_role"]))))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--dataset", action="append", help="Limit to one or more dataset identifiers")
    parser.add_argument("--sample", action="append", help="Limit to one or more sample identifiers")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_manifest(args.manifest)
    if args.dataset:
        rows = [row for row in rows if row["dataset"] in set(args.dataset)]
    if args.sample:
        rows = [row for row in rows if row["sample"] in set(args.sample)]
    if args.dry_run:
        total = sum(int(row["expected_size"] or 0) for row in rows)
        print(f"files={len(rows)} expected_bytes={total}")
        for row in rows:
            print(row["dataset"], row["sample"], row["file_role"], row["url"])
        return

    completed: list[dict[str, object]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = {executor.submit(download_one, row, args.data_root): row for row in rows}
        for future in as_completed(futures):
            result = future.result()
            completed.append(result)
            print(
                f"[{len(completed)}/{len(rows)}] {result['status']} "
                f"{result['dataset']} {result['sample']} {result['file_role']} {result['bytes']}"
            )
    write_log(args.log, completed)
    print(f"wrote {args.log}")


if __name__ == "__main__":
    main()
