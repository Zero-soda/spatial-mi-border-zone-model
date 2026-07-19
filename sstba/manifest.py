"""Run-manifest helpers for reproducible SSTBA executions."""

from __future__ import annotations

import hashlib
import importlib.metadata
import platform
from datetime import datetime, timezone
from pathlib import Path


def sha256_file(path: Path) -> str:
    """Return the SHA-256 checksum of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dependency_versions(names: tuple[str, ...] = ("numpy", "pandas", "scipy")) -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not-installed"
    return versions


def base_manifest(*, mode: str, output_directory: Path, seed: int) -> dict[str, object]:
    """Create the stable manifest envelope for one run."""

    # Outputs are always described relative to the directory containing the
    # manifest so that identical runs remain comparable across machines.
    del output_directory
    return {
        "sstba_version": "0.5.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode": mode,
        "seed": int(seed),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "dependencies": dependency_versions(),
        "output_directory": ".",
        "warnings": [],
    }
