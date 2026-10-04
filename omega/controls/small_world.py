"""Watts–Strogatz propio con rng (OMEGA_1_1_DESIGN §1.10).

Referencia: Watts & Strogatz (1998), Nature 393:440.
"""

from __future__ import annotations

import numpy as np

from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray

__all__ = ["ring_lattice", "watts_strogatz", "matched_small_world"]


def ring_lattice(n: int, k: int) -> FloatArray:
    """Retícula en anillo: cada nodo unido a sus k/2 vecinos a cada lado (k par, 2 <= k < n)."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 3:
        raise ValueError(f"n debe ser entero >= 3, recibido {n!r}")
    if isinstance(k, bool) or not isinstance(k, int) or k < 2 or k % 2 != 0 or k >= n:
        raise ValueError(f"k debe ser par con 2 <= k < n, recibido {k!r}")
    a = np.zeros((n, n), dtype=np.float64)
    idx = np.arange(n)
    for s in range(1, k // 2 + 1):
        a[idx, (idx + s) % n] = 1.0
        a[(idx + s) % n, idx] = 1.0
    return a


def watts_strogatz(n: int, k: int, p: float, rng: np.random.Generator) -> FloatArray:
    """Watts–Strogatz: cada arista (i, i+s) se recablea con prob. p a un destino uniforme
    que no sea i ni vecino actual de i. Conserva el número de aristas; p=0 es el anillo."""
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    if isinstance(p, bool) or not isinstance(p, (int, float)) or not 0.0 <= p <= 1.0:
        raise ValueError(f"p debe estar en [0, 1], recibido {p!r}")
    a = ring_lattice(n, k)
    if p == 0.0:
        return a
    for s in range(1, k // 2 + 1):
        for i in range(n):
            if rng.random() >= p:
                continue
            j = (i + s) % n
            if a[i, j] == 0.0:  # ya recableada por otro paso
                continue
            free = np.flatnonzero(a[i] == 0.0)
            free = free[free != i]
            if free.size == 0:
                continue
            t = int(free[rng.integers(0, free.size)])
            a[i, j] = a[j, i] = 0.0
            a[i, t] = a[t, i] = 1.0
    return a


def matched_small_world(w: FloatArray, w_min: float, p: float, rng: np.random.Generator) -> FloatArray:
    """WS con k = par más cercano al grado medio binario de (W > w_min) (mínimo 2, máximo n-2 par)."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min debe estar en [0,1]")
    n = w.shape[0]
    if n < 4:
        raise ValueError("matched_small_world requiere n >= 4")
    kbar = float((w > w_min).sum()) / n
    k = 2 * int(round(kbar / 2.0))
    k = max(2, k)
    kmax = (n - 1) if (n - 1) % 2 == 0 else n - 2
    k = min(k, kmax)
    return watts_strogatz(n, k, p, rng)
