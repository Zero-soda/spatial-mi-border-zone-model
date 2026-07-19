#!/usr/bin/env python3
"""Resolve the analysis root in both the monorepo and standalone release."""

from __future__ import annotations

from pathlib import Path


def project_root(script_file: str | Path) -> Path:
    """Return the root containing analysis inputs or the standalone release.

    The manuscript was developed in a larger workspace, whereas public users
    receive the release directory as the repository root. Prefer a parent that
    already contains ``data/raw`` and ``results``; otherwise fall back to the
    parent containing the release ``requirements.txt`` and ``scripts`` directory.
    """
    script_path = Path(script_file).resolve()
    for parent in script_path.parents:
        if (parent / "data" / "raw").is_dir() and (parent / "results").is_dir():
            return parent
    for parent in script_path.parents:
        if (parent / "requirements.txt").is_file() and (parent / "scripts").is_dir():
            return parent
    raise RuntimeError(f"Could not resolve project root from {script_path}")
