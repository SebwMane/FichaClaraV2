"""Grafos geométricos aleatorios: controles positivos (OMEGA_1_1_DESIGN §1.10).

RGG en el toro T^d y en la esfera S^d; el ángulo de casquete se obtiene por bisección (1e-12).
Referencia: Penrose (2003), "Random Geometric Graphs", Oxford UP.
"""

from __future__ import annotations

import math
from typing import Literal

import networkx as nx
import numpy as np

from omega.types import FloatArray

__all__ = ["rgg_torus", "cap_angle", "rgg_sphere", "balanced_tree"]


def _check(n: int, dim: int, mean_degree: float, rng: np.random.Generator | None = None) -> None:
    if isinstance(n, bool) or not isinstance(n, int) or n < 3:
        raise ValueError(f"n debe ser entero >= 3, recibido {n!r}")
    if isinstance(dim, bool) or not isinstance(dim, int) or dim < 1:
        raise ValueError(f"dim debe ser entero >= 1, recibido {dim!r}")
    if isinstance(mean_degree, bool) or not isinstance(mean_degree, (int, float)):
        raise TypeError("mean_degree debe ser real")
    if not 0.0 < mean_degree <= n - 1:
        raise ValueError(f"mean_degree debe estar en (0, n-1], recibido {mean_degree!r}")
    if rng is not None and not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")


def rgg_torus(
    n: int,
    dim: int,
    mean_degree: float,
    rng: np.random.Generator,
    weights: Literal["binary", "euclidean"] = "binary",
) -> FloatArray:
    """RGG en el toro [0,1)^dim con métrica periódica y radio r tal que V_d r^d (n-1) = k̄.

    `binary`: peso 1 si dist < r. `euclidean`: peso 1 - dist/r en (0, 1].
    """
    _check(n, dim, mean_degree, rng)
    if weights not in ("binary", "euclidean"):
        raise ValueError(f"weights inválido: {weights!r}")
    vol_unit = math.pi ** (dim / 2.0) / math.gamma(dim / 2.0 + 1.0)
    r = (mean_degree / ((n - 1) * vol_unit)) ** (1.0 / dim)
    if r >= 0.5:
        raise ValueError("radio >= 0.5: mean_degree demasiado grande para el toro")
    x = rng.random((n, dim))
    diff = np.abs(x[:, None, :] - x[None, :, :])
    diff = np.minimum(diff, 1.0 - diff)
    dist = np.sqrt((diff**2).sum(axis=2))
    adj = dist < r
    np.fill_diagonal(adj, False)
    if weights == "binary":
        return np.asarray(adj, dtype=np.float64)
    out: FloatArray = np.where(adj, 1.0 - dist / r, 0.0)
    out = np.maximum(out, np.where(adj, np.finfo(np.float64).tiny, 0.0))
    return np.asarray(out, dtype=np.float64)


def _cap_fraction(r: float, dim: int) -> float:
    """Fracción del área de S^dim dentro de un casquete de ángulo r (dim = 1, 2, 3)."""
    if dim == 1:
        return r / math.pi
    if dim == 2:
        return (1.0 - math.cos(r)) / 2.0
    if dim == 3:
        return (r - math.sin(r) * math.cos(r)) / math.pi
    raise ValueError("S^dim soportada solo para dim en {1, 2, 3}")


def cap_angle(n: int, dim: int, mean_degree: float) -> float:
    """Ángulo r de casquete en S^dim con fracción de área = k̄/(n-1), por bisección (1e-12)."""
    _check(n, dim, mean_degree)
    if dim > 3:
        raise ValueError("S^dim soportada solo para dim en {1, 2, 3}")
    target = mean_degree / (n - 1)
    lo, hi = 0.0, math.pi
    while hi - lo > 1e-12:
        mid = 0.5 * (lo + hi)
        if _cap_fraction(mid, dim) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rgg_sphere(n: int, dim: int, mean_degree: float, rng: np.random.Generator) -> FloatArray:
    """RGG binario en S^dim ⊂ R^(dim+1): arista si el ángulo geodésico < cap_angle."""
    _check(n, dim, mean_degree, rng)
    r = cap_angle(n, dim, mean_degree)
    x = rng.standard_normal((n, dim + 1))
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    cosang = np.clip(x @ x.T, -1.0, 1.0)
    adj = np.arccos(cosang) < r
    np.fill_diagonal(adj, False)
    return adj.astype(np.float64)


def balanced_tree(branching: int, height: int) -> FloatArray:
    """Árbol balanceado r-ario (nulo red-team): matriz de adyacencia 0/1."""
    for name, v, lo in (("branching", branching, 2), ("height", height, 1)):
        if isinstance(v, bool) or not isinstance(v, int) or v < lo:
            raise ValueError(f"{name} debe ser entero >= {lo}, recibido {v!r}")
    g = nx.balanced_tree(branching, height)
    return np.asarray(nx.to_numpy_array(g, nodelist=sorted(g.nodes), weight=None), dtype=np.float64)
