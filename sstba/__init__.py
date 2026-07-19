"""Spatial-state transfer and boundary analysis (SSTBA)."""

from .boundary import BoundaryResult, analyse_boundary
from .scoring import score_modules

__all__ = ["BoundaryResult", "analyse_boundary", "score_modules"]
__version__ = "0.5.0"
