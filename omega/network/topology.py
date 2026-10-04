"""Observables topologicos de W (D-1, D-2, D-10, D-11; R2, R7)."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

from omega.network.weights import upper_triangle, validate_weight_matrix
from omega.types import BoolArray, FloatArray, IntArray, TopologyObservables

__all__ = [
    "adjacency",
    "strength",
    "binary_degree",
    "component_labels",
    "component_sizes",
    "giant_component_nodes",
    "submatrix",
    "weighted_clustering",
    "binary_clustering",
    "topology_observables",
]


def _check_adjacency(a: BoolArray) -> None:
    if not isinstance(a, np.ndarray) or a.dtype != np.bool_:
        raise TypeError("a debe ser ndarray bool")
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"a debe ser cuadrada, forma {a.shape}")
    if np.any(np.diag(a)) or not np.array_equal(a, a.T):
        raise ValueError("a debe ser simetrica con diagonal False")


def adjacency(w: FloatArray, w_min: float) -> BoolArray:
    """A_ij = (W_ij >= w_min) con diagonal False; w_min en [0,1] (D-2)."""
    validate_weight_matrix(w)
    if not 0.0 <= w_min <= 1.0:
        raise ValueError(f"w_min debe estar en [0,1], recibido {w_min}")
    a = np.asarray(w >= w_min, dtype=np.bool_)
    np.fill_diagonal(a, False)
    return a


def strength(w: FloatArray) -> FloatArray:
    """Fuerza k_i = sum_j W_ij."""
    validate_weight_matrix(w)
    return np.asarray(w.sum(axis=1), dtype=np.float64)


def binary_degree(a: BoolArray) -> IntArray:
    """Grado binario por nodo."""
    _check_adjacency(a)
    return np.asarray(a.sum(axis=1), dtype=np.int64)


def component_labels(a: BoolArray) -> tuple[int, IntArray]:
    """(numero de componentes, etiqueta por nodo) de la adyacencia binaria."""
    _check_adjacency(a)
    n_comp, labels = connected_components(csr_array(a), directed=False)
    return int(n_comp), np.asarray(labels, dtype=np.int64)


def component_sizes(labels: IntArray) -> IntArray:
    """Tamano de cada componente indexado por etiqueta (0..n_comp-1)."""
    lab = np.asarray(labels)
    if lab.ndim != 1 or lab.size == 0 or lab.min() < 0:
        raise ValueError("labels debe ser vector no vacio de enteros >= 0")
    return np.asarray(np.bincount(lab), dtype=np.int64)


def giant_component_nodes(labels: IntArray) -> IntArray:
    """Nodos (ordenados) de la mayor componente; desempate por menor etiqueta."""
    sizes = component_sizes(labels)
    g = int(np.argmax(sizes))
    return np.asarray(np.flatnonzero(np.asarray(labels) == g), dtype=np.int64)


def submatrix(w: FloatArray, nodes: IntArray) -> FloatArray:
    """Submatriz W[nodes, nodes] (copia)."""
    n = np.asarray(nodes)
    if n.ndim != 1 or n.size == 0 or not np.issubdtype(n.dtype, np.integer):
        raise ValueError("nodes debe ser vector entero no vacio")
    if n.min() < 0 or n.max() >= w.shape[0] or np.unique(n).size != n.size:
        raise ValueError("nodes fuera de rango o con repetidos")
    return np.asarray(w[np.ix_(n, n)], dtype=np.float64).copy()


def weighted_clustering(w: FloatArray) -> float:
    """C_W = tr(W^3) / sum_i (k_i^2 - sum_j W_ij^2) en [0,1]; 0 si no hay tripletes (D-10)."""
    validate_weight_matrix(w)
    k = w.sum(axis=1)
    open_triplets = float(np.sum(k * k) - np.sum(w * w))
    if open_triplets <= 0.0:
        return 0.0
    closed = float(np.sum((w @ w) * w))  # tr(W^3)
    return float(min(1.0, max(0.0, closed / open_triplets)))


def binary_clustering(a: BoolArray) -> float:
    """Media de coeficientes locales sobre los N nodos (0 si grado<2) (D-10)."""
    _check_adjacency(a)
    af = a.astype(np.float64)
    k = af.sum(axis=1)
    tri = np.sum((af @ af) * af, axis=1)  # 2 * triangulos por nodo
    possible = k * (k - 1.0)
    local = np.divide(tri, possible, out=np.zeros_like(tri), where=possible > 0.0)
    return float(local.mean())


def topology_observables(w: FloatArray, w_min: float, large_component_frac: float) -> TopologyObservables:
    """Observables topologicos; componente 'grande' si su tamano >= frac*N (D-17, §4.1)."""
    validate_weight_matrix(w)
    if not 0.0 <= large_component_frac <= 1.0:
        raise ValueError("large_component_frac debe estar en [0,1]")
    n = w.shape[0]
    a = adjacency(w, w_min)
    k = strength(w)
    n_comp, labels = component_labels(a)
    sizes = component_sizes(labels)
    giant = int(sizes.max())
    mean_k = float(k.mean())
    std_k = float(k.std())
    edges = upper_triangle(w)
    return TopologyObservables(
        n=n,
        mean_strength=mean_k,
        std_strength=std_k,
        cv_strength=std_k / mean_k if mean_k > 0.0 else 0.0,
        mean_binary_degree=float(binary_degree(a).mean()),
        binary_density=float(a.sum()) / (n * (n - 1)),
        mean_weight=float(edges.mean()),
        giant_size=giant,
        giant_fraction=giant / n,
        n_components=n_comp,
        n_large_components=int(np.sum(sizes >= large_component_frac * n)),
        clustering_weighted=weighted_clustering(w),
        clustering_binary=binary_clustering(a),
        frac_at_zero=float(np.mean(edges == 0.0)),
        frac_at_one=float(np.mean(edges == 1.0)),
    )
