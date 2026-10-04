"""Nulos Erdős–Rényi (OMEGA_1_1_DESIGN §1.10, P§10).

Referencia: Erdős & Rényi (1959), "On random graphs I", Publ. Math. Debrecen 6.
Funciones puras con RNG explícito; devuelven matrices simétricas (n, n) con diagonal 0.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from omega.network.weights import from_upper_triangle, validate_weight_matrix
from omega.types import FloatArray

__all__ = ["erdos_renyi_gnp", "erdos_renyi_gnm", "matched_erdos_renyi"]


def _check(n: int, rng: np.random.Generator) -> None:
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError(f"n debe ser entero >= 2, recibido {n!r}")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")


def erdos_renyi_gnp(n: int, p: float, rng: np.random.Generator) -> FloatArray:
    """G(n, p) binario: cada par se conecta con probabilidad p (peso 1.0)."""
    _check(n, rng)
    if isinstance(p, bool) or not isinstance(p, (int, float)) or not 0.0 <= p <= 1.0:
        raise ValueError(f"p debe estar en [0, 1], recibido {p!r}")
    n_pairs = n * (n - 1) // 2
    v = (rng.random(n_pairs) < p).astype(np.float64)
    return from_upper_triangle(v, n)


def erdos_renyi_gnm(n: int, m: int, rng: np.random.Generator) -> FloatArray:
    """G(n, m) binario: exactamente m aristas elegidas uniformemente sin reemplazo."""
    _check(n, rng)
    n_pairs = n * (n - 1) // 2
    if isinstance(m, bool) or not isinstance(m, int) or not 0 <= m <= n_pairs:
        raise ValueError(f"m debe ser entero en [0, {n_pairs}], recibido {m!r}")
    v = np.zeros(n_pairs, dtype=np.float64)
    v[rng.choice(n_pairs, size=m, replace=False)] = 1.0
    return from_upper_triangle(v, n)


def matched_erdos_renyi(
    w: FloatArray,
    w_min: float,
    rng: np.random.Generator,
    weights: Literal["binary", "shuffled"] = "binary",
) -> FloatArray:
    """G(n, m) con el mismo m = #{W > w_min} que el candidato.

    `binary`: pesos 1.0. `shuffled`: los pesos de las m aristas del candidato se reparten
    al azar entre las m aristas nuevas (conserva el multiconjunto de pesos).
    """
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError("w_min debe estar en [0,1]")
    if weights not in ("binary", "shuffled"):
        raise ValueError(f"weights debe ser 'binary' o 'shuffled', recibido {weights!r}")
    n = w.shape[0]
    iu = np.triu_indices(n, k=1)
    vals = w[iu]
    keep = vals > w_min
    m = int(keep.sum())
    out = erdos_renyi_gnm(n, m, rng)
    if weights == "shuffled" and m > 0:
        picked = rng.permutation(vals[keep])
        v = np.zeros(vals.shape[0], dtype=np.float64)
        mask = out[iu] > 0.0
        v[mask] = picked
        out = from_upper_triangle(v, n)
    return out
