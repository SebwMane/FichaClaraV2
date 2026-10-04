"""Filtracion por umbral theta de la matriz de pesos W (OMEGA_1_1_DESIGN §1.2, §3.3).

G_theta = (V, {W > theta}) con umbral estricto, como en las distancias de Omega-1.0
(`omega.geometry.distances.threshold_adjacency`). Funciones puras, sin E/S.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

from omega.contracts import FiltrationLevel
from omega.network.weights import validate_weight_matrix
from omega.types import BoolArray, FloatArray, IntArray

__all__ = ["threshold_graph", "descending_edges", "component_labels", "filtration_levels"]


def threshold_graph(w: FloatArray, theta: float) -> BoolArray:
    """Adyacencia binaria A_ij = (W_ij > theta), simetrica y con diagonal False.

    theta debe estar en [0, 1]. Umbral estricto (M§8, D-2).
    """
    validate_weight_matrix(w)
    if isinstance(theta, bool) or not 0.0 <= float(theta) <= 1.0:
        raise ValueError(f"theta debe estar en [0,1], recibido {theta}")
    a = np.asarray(w > theta, dtype=np.bool_)
    np.fill_diagonal(a, False)
    return a


def descending_edges(w: FloatArray) -> tuple[IntArray, IntArray, FloatArray]:
    """Aristas (i<j) con W>0 ordenadas por (-W, i, j): orden de la filtracion descendente.

    Devuelve (i, j, peso) como arrays int64/int64/float64. El desempate por (i, j) hace el
    orden determinista (usado por la persistencia H0 por Kruskal).
    """
    validate_weight_matrix(w)
    iu, ju = np.triu_indices(w.shape[0], k=1)
    vals = w[iu, ju]
    keep = vals > 0.0
    i = iu[keep].astype(np.int64)
    j = ju[keep].astype(np.int64)
    v = np.asarray(vals[keep], dtype=np.float64)
    order = np.lexsort((j, i, -v))
    return i[order], j[order], v[order]


def component_labels(a: BoolArray) -> tuple[int, IntArray]:
    """Numero de componentes conexas (incluidos nodos aislados) y etiqueta por nodo."""
    n_comp, labels = connected_components(csr_array(a.astype(np.int8)), directed=False)
    return int(n_comp), np.asarray(labels, dtype=np.int64)


def filtration_levels(w: FloatArray, thetas: Sequence[float]) -> tuple[FiltrationLevel, ...]:
    """Resumen de G_theta para cada theta: aristas, beta0, fraccion gigante y grado medio."""
    validate_weight_matrix(w)
    n = w.shape[0]
    out: list[FiltrationLevel] = []
    for theta in thetas:
        a = threshold_graph(w, float(theta))
        n_comp, labels = component_labels(a)
        giant = int(np.bincount(labels).max())
        n_edges = int(a.sum()) // 2
        out.append(
            FiltrationLevel(
                theta=float(theta),
                n_edges=n_edges,
                n_components=n_comp,
                giant_fraction=giant / n,
                mean_degree=2.0 * n_edges / n,
            )
        )
    return tuple(out)
