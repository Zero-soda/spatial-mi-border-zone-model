#!/usr/bin/env python3
"""Generate the release-wide checksum and table-row manifest."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


EXCLUDED_PARTS = {".git", ".pytest_cache", "__pycache__"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def table_rows(path: Path) -> str:
    if path.suffix.lower() not in {".tsv", ".csv"}:
        return "NA"
    with path.open("rb") as handle:
        lines = sum(1 for _ in handle)
    return str(max(lines - 1, 0))


def included(path: Path, manifest: Path) -> bool:
    if path == manifest:
        return False
    if any(part in EXCLUDED_PARTS for part in path.parts):
        return False
    return path.suffix.lower() not in EXCLUDED_SUFFIXES


def build_manifest(root: Path, output: Path) -> None:
    rows = ["file\trows_if_table\tsha256"]
    for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        if not included(path, output):
            continue
        relative = path.relative_to(root).as_posix()
        rows.append(f"{relative}\t{table_rows(path)}\t{sha256(path)}")
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", default="FILE_MANIFEST.tsv")
    args = parser.parse_args()
    root = Path(args.root).expanduser().resolve()
    output = Path(args.output).expanduser()
    if not output.is_absolute():
        output = root / output
    build_manifest(root, output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
