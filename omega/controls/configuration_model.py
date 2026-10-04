"""Modelo de configuración "erased" (OMEGA_1_1_DESIGN §1.10).

Referencia: Molloy & Reed (1995); Newman (2010), "Networks", cap. 13. Se usa networkx con
semilla derivada del rng y se borran lazos y multiaristas (grados resultantes <= objetivo).
"""

from __future__ import annotations

import networkx as nx
import numpy as np

from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray, IntArray

__all__ = ["configuration_model", "matched_configuration_model"]


def configuration_model(degrees: IntArray, rng: np.random.Generator) -> FloatArray:
    """Grafo binario simple con grados <= `degrees` (suma par), sin lazos ni multiaristas."""
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    d = np.asarray(degrees)
    if d.ndim != 1 or d.size < 2:
        raise ValueError("degrees debe ser 1D con al menos 2 nodos")
    if not np.issubdtype(d.dtype, np.integer):
        raise TypeError("degrees debe ser entero")
    if (d < 0).any() or (d > d.size - 1).any():
        raise ValueError("grados fuera de [0, n-1]")
    if int(d.sum()) % 2 != 0:
        raise ValueError("la suma de grados debe ser par")
    n = int(d.size)
    out = np.zeros((n, n), dtype=np.float64)
    if int(d.sum()) == 0:
        return out
    seed = int(rng.integers(0, 2**31 - 1))
    multi = nx.configuration_model([int(x) for x in d], seed=seed)
    for a, b in multi.edges():
        if a != b:
            out[a, b] = out[b, a] = 1.0
    _repair(out, d, rng)
    return out


def _repair(out: FloatArray, target: IntArray, rng: np.random.Generator) -> None:
    """Reasigna stubs perdidos por el borrado de lazos/multiaristas uniendo pares de nodos
    con déficit que aún no son vecinos (sigue siendo simple y con grados <= objetivo)."""
    for _ in range(50):
        deficit = target - out.sum(axis=1).astype(np.int64)
        stubs = np.repeat(np.arange(target.size), np.clip(deficit, 0, None))
        if stubs.size < 2:
            return
        stubs = rng.permutation(stubs)
        added = False
        for a, b in zip(stubs[0::2], stubs[1::2]):
            if a != b and out[a, b] == 0.0:
                out[a, b] = out[b, a] = 1.0
                added = True
        if not added:
            return


def matched_configuration_model(w: FloatArray, w_min: float, rng: np.random.Generator) -> FloatArray:
    """Configuration model con la secuencia de grados binaria de A = (W > w_min)."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min debe estar en [0,1]")
    deg = (w > w_min).sum(axis=1).astype(np.int64)
    if int(deg.sum()) % 2 != 0:  # imposible en A simétrica; defensa
        raise ValueError("A no es simétrica")
    return configuration_model(deg, rng)
