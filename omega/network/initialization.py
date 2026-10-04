"""Inicializacion aleatoria de W0 (D-0, M§12): U(0,1) en la triangular superior, espejada."""

from __future__ import annotations

import numpy as np

from omega.network.weights import from_upper_triangle
from omega.types import FloatArray

__all__ = ["random_uniform_weights"]


def random_uniform_weights(n: int, rng: np.random.Generator) -> FloatArray:
    """W0 simetrica (n,n), diag 0, con las n(n-1)/2 aristas i.i.d. U(0,1); n>=2."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError(f"n debe ser entero >= 2, recibido {n!r}")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    return from_upper_triangle(rng.random(n * (n - 1) // 2), n)
