"""Curvatura de Ollivier perezosa sobre una arista, version dispersa (auditoria K, OMEGA_K_AUDIT_PRERREGISTRO §1).

Equivalente exacta a `omega.curvature.ollivier.ollivier_edge` con W = A binaria y coste en saltos, sin matriz densa.
Soporte de m_x: {x} U N(x); idem y. Para una arista (x, y), toda pareja (a, b) con a en el soporte de x y b en el de y
cumple d(a, b) <= 3 (a-x-y-b), luego la distancia de saltos acotada a profundidad 3 se obtiene exactamente de:
0 si a == b; 1 si a ~ b; 2 si comparten un vecino; 3 en otro caso (equivale al BFS limitado a profundidad 3).
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from omega.curvature.ollivier import wasserstein1
from omega.types import FloatArray, IntArray

__all__ = ["ollivier_edge_sparse"]


def _support(adj: sparse.csr_array, x: int, idleness: float) -> tuple[IntArray, FloatArray]:
    nb = np.asarray(adj.indices[adj.indptr[x] : adj.indptr[x + 1]], dtype=np.int64)
    nb = nb[nb != x]
    if nb.size == 0:
        return np.asarray([x], dtype=np.int64), np.asarray([1.0], dtype=np.float64)
    support = np.concatenate([[x], nb]).astype(np.int64)
    mass = np.concatenate([[idleness], np.full(nb.size, (1.0 - idleness) / nb.size)]).astype(np.float64)
    return support, mass


def _hop_cost(adj: sparse.csr_array, sx: IntArray, sy: IntArray) -> FloatArray:
    """Distancias de saltos (<= 3 garantizado para soportes de una arista) entre sx y sy."""
    ax = adj[sx, :]
    ay = adj[sy, :]
    one = np.asarray((ax[:, sy] != 0).toarray(), dtype=bool)
    two = np.asarray(((ax @ ay.T) != 0).toarray(), dtype=bool)
    cost = np.full((sx.size, sy.size), 3.0, dtype=np.float64)
    cost[two] = 2.0
    cost[one] = 1.0
    cost[sx[:, None] == sy[None, :]] = 0.0
    return cost


def ollivier_edge_sparse(adj: sparse.csr_array, x: int, y: int, idleness: float) -> float:
    """kappa(x, y) = 1 - W1(m_x, m_y) / d(x, y) con d(x, y) = 1 (arista), adyacencia binaria dispersa simetrica."""
    if not 0.0 <= idleness < 1.0:
        raise ValueError("idleness debe estar en [0,1)")
    if x == y or adj[x, y] == 0:
        raise ValueError("(x, y) debe ser una arista")
    sx, mx = _support(adj, x, idleness)
    sy, my = _support(adj, y, idleness)
    return 1.0 - wasserstein1(mx, my, _hop_cost(adj, sx, sy))
