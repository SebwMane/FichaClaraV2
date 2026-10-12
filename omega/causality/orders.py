"""Ordenes parciales finitos: 2-ordenes aleatorios y controles no geometricos (bloque C0)."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["random_2order", "kr_order", "transitive_percolation", "transitive_closure", "is_partial_order"]

BoolMatrix = npt.NDArray[np.bool_]


def transitive_closure(adj: BoolMatrix) -> BoolMatrix:
    """Cierre transitivo estricto de un DAG por cuadrados booleanos repetidos."""
    if adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("adj debe ser cuadrada")
    c = adj.astype(np.float32)
    while True:
        nxt = ((c @ c) > 0) | (c > 0)
        if np.array_equal(nxt, c > 0):
            out: BoolMatrix = nxt
            return out
        c = nxt.astype(np.float32)


def is_partial_order(c: BoolMatrix) -> bool:
    """Irreflexivo, antisimetrico y transitivo."""
    if np.any(np.diag(c)):
        return False
    if np.any(c & c.T):
        return False
    cf = c.astype(np.float32)
    return not bool(np.any(((cf @ cf) > 0) & ~c))


def random_2order(n: int, rng: np.random.Generator) -> BoolMatrix:
    """Interseccion de dos ordenes totales uniformes: i < j si u_i < u_j y v_i < v_j."""
    u = rng.permutation(n)
    v = rng.permutation(n)
    rel: BoolMatrix = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
    return rel


def kr_order(n: int, rng: np.random.Generator, p: float = 0.5) -> BoolMatrix:
    """Orden de Kleitman-Rothschild: capas n/4, n/2, n/4; relaciones entre capas consecutivas con prob. p."""
    n1 = n // 4
    n3 = n // 4
    n2 = n - n1 - n3
    layer = np.concatenate([np.zeros(n1, int), np.ones(n2, int), np.full(n3, 2)])
    adj = (layer[:, None] + 1 == layer[None, :]) & (rng.random((n, n)) < p)
    perm = rng.permutation(n)  # los indices no deben codificar la capa
    adj = adj[np.ix_(perm, perm)]
    return transitive_closure(adj)


def transitive_percolation(n: int, p: float, rng: np.random.Generator) -> BoolMatrix:
    """Aristas i<j con prob. p y cierre transitivo; los indices se permutan al final."""
    upper = np.triu(rng.random((n, n)) < p, k=1)
    clo = transitive_closure(upper)
    perm = rng.permutation(n)
    out: BoolMatrix = clo[np.ix_(perm, perm)]
    return out
