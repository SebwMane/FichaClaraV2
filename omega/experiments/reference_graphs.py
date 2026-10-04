"""Redes de referencia para validacion (D-16, ANALYSIS §5.11). SOLO para experimentos y pruebas:
prohibido importarlo desde omega.dynamics y omega.network."""

from __future__ import annotations

from math import gamma, pi

import networkx as nx
import numpy as np

from omega.types import FloatArray

__all__ = ["periodic_lattice", "open_lattice", "random_geometric_torus", "complete_graph", "random_regular"]


def _check_shape(shape: tuple[int, ...], minimum: int) -> None:
    if not isinstance(shape, tuple) or not 1 <= len(shape) <= 3:
        raise ValueError("shape debe ser tupla de 1 a 3 enteros")
    for s in shape:
        if isinstance(s, bool) or not isinstance(s, int) or s < minimum:
            raise ValueError(f"cada lado debe ser entero >= {minimum}, recibido {shape}")


def _lattice(shape: tuple[int, ...], periodic: bool) -> FloatArray:
    n = int(np.prod(shape))
    w = np.zeros((n, n), dtype=np.float64)
    coords = np.stack(np.unravel_index(np.arange(n), shape), axis=1)
    for axis, size in enumerate(shape):
        nxt = coords.copy()
        nxt[:, axis] += 1
        if periodic:
            nxt[:, axis] %= size
            keep = np.ones(n, dtype=np.bool_)
        else:
            keep = nxt[:, axis] < size
        src = np.flatnonzero(keep)
        dst = np.ravel_multi_index(tuple(nxt[src].T), shape)
        w[src, dst] = 1.0
        w[dst, src] = 1.0
    return w


def periodic_lattice(shape: tuple[int, ...]) -> FloatArray:
    """Anillo (n,), toro 2D (n,m) o 3D con W=1 entre vecinos mas proximos; lados >= 3 (D-16)."""
    _check_shape(shape, 3)
    return _lattice(shape, True)


def open_lattice(shape: tuple[int, ...]) -> FloatArray:
    """Reticula con bordes abiertos y W=1 (solo informe, D-16); lados >= 2."""
    _check_shape(shape, 2)
    return _lattice(shape, False)


def random_geometric_torus(
    n: int, dim: int, mean_degree: float, rng: np.random.Generator, *, euclidean: bool = True
) -> FloatArray:
    """Grafo geometrico aleatorio en el toro unidad [0,1)^dim: arista si distancia toroidal <= R, con
    grado medio esperado `mean_degree`. euclidean=True: W = min(1, R/(4 d)) (longitud ~ distancia
    euclidiana, M§8); False: W=1 (saltos)."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("n debe ser entero >= 2")
    if dim not in (1, 2, 3):
        raise ValueError("dim debe ser 1, 2 o 3")
    if not mean_degree > 0.0:
        raise ValueError("mean_degree debe ser > 0")
    vol = pi ** (dim / 2) / gamma(dim / 2 + 1)
    radius = (mean_degree / (n * vol)) ** (1.0 / dim)
    pts = rng.random((n, dim))
    diff = np.abs(pts[:, None, :] - pts[None, :, :])
    diff = np.minimum(diff, 1.0 - diff)
    dist = np.sqrt(np.sum(diff * diff, axis=2))
    edge = dist <= radius
    np.fill_diagonal(edge, False)
    if not euclidean:
        return np.asarray(edge, dtype=np.float64)
    w = np.where(edge, np.minimum(1.0, radius / (4.0 * np.where(edge, dist, 1.0))), 0.0)
    return np.asarray(w, dtype=np.float64)


def complete_graph(n: int) -> FloatArray:
    """Grafo completo J - I (W=1)."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("n debe ser entero >= 2")
    return np.asarray(np.ones((n, n)) - np.eye(n), dtype=np.float64)


def random_regular(n: int, degree: int, rng: np.random.Generator) -> FloatArray:
    """Grafo d-regular aleatorio con W=1 (networkx sembrado desde `rng`); n*degree par."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("n debe ser entero >= 2")
    if isinstance(degree, bool) or not isinstance(degree, int) or not 1 <= degree < n:
        raise ValueError("degree debe estar en [1, n-1]")
    if (n * degree) % 2:
        raise ValueError("n*degree debe ser par")
    g = nx.random_regular_graph(degree, n, seed=int(rng.integers(0, 2**31 - 1)))
    return np.asarray(nx.to_numpy_array(g, nodelist=range(n), dtype=np.float64), dtype=np.float64)
