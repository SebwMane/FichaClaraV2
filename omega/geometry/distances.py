"""Distancias sobre la red de pesos (M§8, D-1, D-2, D-11): longitud de arista 1/(W+eps), Dijkstra."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components, dijkstra

from omega.network.weights import validate_weight_matrix
from omega.types import BoolArray, FloatArray, IntArray

__all__ = [
    "threshold_adjacency",
    "edge_length_matrix",
    "distance_matrix",
    "hop_distance_matrix",
    "validate_distance_matrix",
    "giant_nodes_from_distances",
    "ball",
    "ball_counts",
    "path_length",
    "diameter",
]


def threshold_adjacency(w: FloatArray, w_min: float) -> BoolArray:
    """A_ij = (W_ij > w_min), diagonal False (M§8: umbral estricto, D-2)."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError(f"w_min debe estar en [0,1], recibido {w_min}")
    a = np.asarray(w > w_min, dtype=np.bool_)
    np.fill_diagonal(a, False)
    return a


def edge_length_matrix(w: FloatArray, a: BoolArray, epsilon: float) -> csr_array:
    """Matriz dispersa de longitudes l_ij = 1/(W_ij+eps) sobre las aristas de `a` (M§8, D-1)."""
    validate_weight_matrix(w)
    if not isinstance(a, np.ndarray) or a.dtype != np.bool_ or a.shape != w.shape:
        raise ValueError("a debe ser ndarray bool con la misma forma que w")
    if not np.array_equal(a, a.T) or np.any(np.diag(a)):
        raise ValueError("a debe ser simetrica con diagonal False")
    if not epsilon > 0.0:
        raise ValueError(f"epsilon debe ser > 0, recibido {epsilon}")
    lengths = np.where(a, 1.0 / (w + epsilon), 0.0)
    return csr_array(lengths)


def distance_matrix(w: FloatArray, w_min: float, epsilon: float) -> FloatArray:
    """Distancias de camino minimo (Dijkstra) sobre aristas con W>w_min; inf entre componentes (M§8)."""
    a = threshold_adjacency(w, w_min)
    g = edge_length_matrix(w, a, epsilon)
    d = np.asarray(dijkstra(g, directed=False), dtype=np.float64)
    return np.asarray(0.5 * (d + d.T), dtype=np.float64)  # exacta simetria (orden de suma en Dijkstra)


def hop_distance_matrix(a: BoolArray) -> FloatArray:
    """Distancia en saltos (BFS) sobre la adyacencia binaria; inf entre componentes (D-22)."""
    if not isinstance(a, np.ndarray) or a.dtype != np.bool_:
        raise TypeError("a debe ser ndarray bool")
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"a debe ser cuadrada, forma {a.shape}")
    if not np.array_equal(a, a.T) or np.any(np.diag(a)):
        raise ValueError("a debe ser simetrica con diagonal False")
    d = dijkstra(csr_array(a.astype(np.float64)), directed=False, unweighted=True)
    return np.asarray(d, dtype=np.float64)


def validate_distance_matrix(d: FloatArray) -> None:
    """Exige (N,N) float64, N>=1, sin NaN, >=0 (inf permitido), simetrica (rtol 1e-9) y diagonal 0."""
    if not isinstance(d, np.ndarray) or d.dtype != np.float64:
        raise TypeError("d debe ser ndarray float64")
    if d.ndim != 2 or d.shape[0] != d.shape[1] or d.shape[0] < 1:
        raise ValueError(f"d debe ser cuadrada no vacia, forma {d.shape}")
    if np.any(np.isnan(d)) or np.any(d < 0.0):
        raise ValueError("d contiene NaN o valores negativos")
    if np.any(np.diag(d) != 0.0):
        raise ValueError("la diagonal de d debe ser 0")
    if not np.allclose(d, d.T, rtol=1e-9, atol=1e-12):
        raise ValueError("d no es simetrica")


def giant_nodes_from_distances(d: FloatArray) -> IntArray:
    """Nodos (ordenados) de la mayor clase de pares con distancia finita; desempate: menor etiqueta (D-11)."""
    validate_distance_matrix(d)
    fin = np.isfinite(d)
    np.fill_diagonal(fin, False)
    _, labels = connected_components(csr_array(fin), directed=False)
    giant = int(np.argmax(np.bincount(labels)))
    return np.asarray(np.flatnonzero(labels == giant), dtype=np.int64)


def ball(d: FloatArray, i: int, r: float) -> IntArray:
    """Indices j con d_ij <= r (incluye i); bola B_i(r) de M§9."""
    validate_distance_matrix(d)
    if isinstance(i, bool) or not isinstance(i, (int, np.integer)) or not 0 <= i < d.shape[0]:
        raise ValueError(f"i fuera de rango: {i!r}")
    return np.asarray(np.flatnonzero(d[i] <= r), dtype=np.int64)


def ball_counts(d: FloatArray, radii: FloatArray) -> IntArray:
    """Matriz (N,K) de N_i(r_k)=|B_i(r_k)|; monotona no decreciente en r si radii es ascendente."""
    validate_distance_matrix(d)
    r = np.asarray(radii, dtype=np.float64)
    if r.ndim != 1 or r.size == 0 or not np.all(np.isfinite(r)):
        raise ValueError("radii debe ser vector 1D no vacio y finito")
    n = d.shape[0]
    out = np.empty((n, r.size), dtype=np.int64)
    srt = np.sort(d, axis=1)
    for i in range(n):
        out[i] = np.searchsorted(srt[i], r, side="right")
    return out


def _finite_offdiag(d: FloatArray) -> FloatArray:
    """Distancias finitas entre pares distintos de la componente gigante."""
    nodes = giant_nodes_from_distances(d)
    sub = d[np.ix_(nodes, nodes)]
    n = sub.shape[0]
    return np.asarray(sub[~np.eye(n, dtype=np.bool_)], dtype=np.float64)


def path_length(d: FloatArray) -> float:
    """Media de d sobre pares distintos de la componente gigante (D-11); nan si hay <2 nodos."""
    v = _finite_offdiag(d)
    return float(v.mean()) if v.size else float("nan")


def diameter(d: FloatArray) -> float:
    """Maximo de d sobre pares de la componente gigante (D-11); nan si hay <2 nodos."""
    v = _finite_offdiag(d)
    return float(v.max()) if v.size else float("nan")
