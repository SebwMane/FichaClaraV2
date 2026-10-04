"""Recableado con grados preservados (OMEGA_1_1_DESIGN §1.10).

Envuelve `degree_preserving_null` de Ω-1.0 con n_swaps = swaps_per_edge · E
(Maslov & Sneppen 2002, Science 296:910).
"""

from __future__ import annotations

import numpy as np

from omega.network.weights import validate_weight_matrix
from omega.phases.stability import degree_preserving_null
from omega.types import FloatArray

__all__ = ["degree_preserving_rewire"]


def degree_preserving_rewire(
    w: FloatArray, w_min: float, rng: np.random.Generator, swaps_per_edge: int = 10
) -> FloatArray:
    """Recableado de A=(W>w_min) con grados exactos y multiconjunto de pesos conservados."""
    validate_weight_matrix(w)
    if isinstance(swaps_per_edge, bool) or not isinstance(swaps_per_edge, int) or swaps_per_edge < 0:
        raise ValueError("swaps_per_edge debe ser entero >= 0")
    n_edges = int((np.triu(w, k=1) > w_min).sum())
    return degree_preserving_null(w, w_min, swaps_per_edge * n_edges, rng)
