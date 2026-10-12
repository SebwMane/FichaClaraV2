"""Ordenes parciales finitos: 2-ordenes aleatorios y controles no geometricos (bloque C0)."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = [
    "random_2order",
    "kr_order",
    "transitive_percolation",
    "thick_layer_order",
    "layered_2order",
    "random_bipartite",
    "transitive_closure",
    "is_partial_order",
]

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


def thick_layer_order(n: int, layers: int, p: float, rng: np.random.Generator) -> BoolMatrix:
    """Capas de tamano n // layers (el resto en la ultima); relaciones entre capas consecutivas con prob. p."""
    if layers < 2 or n < layers:
        raise ValueError("se requieren layers >= 2 y n >= layers")
    size = n // layers
    layer = np.minimum(np.arange(n) // size, layers - 1)
    adj = (layer[:, None] + 1 == layer[None, :]) & (rng.random((n, n)) < p)
    perm = rng.permutation(n)
    return transitive_closure(adj[np.ix_(perm, perm)])


def layered_2order(n: int, k: int, rng: np.random.Generator) -> BoolMatrix:
    """k antichains apiladas y totalmente relacionadas, construidas como 2-orden; etiquetas permutadas."""
    if k < 1 or n < k:
        raise ValueError("se requieren k >= 1 y n >= k")
    size = n // k
    layer = np.minimum(np.arange(n) // size, k - 1)
    u = np.arange(n)
    v = np.empty(n, dtype=np.int64)
    for g in range(k):
        idx = np.flatnonzero(layer == g)
        v[idx] = idx[::-1]  # dentro del bloque: u creciente, v decreciente
    perm = rng.permutation(n)
    u, v = u[perm], v[perm]
    rel: BoolMatrix = (u[:, None] < u[None, :]) & (v[:, None] < v[None, :])
    return rel


def random_bipartite(n: int, p: float, rng: np.random.Generator) -> BoolMatrix:
    """n//2 elementos abajo y el resto arriba; cada par abajo-arriba relacionado con prob. p."""
    lo = n // 2
    rel = np.zeros((n, n), dtype=bool)
    rel[:lo, lo:] = rng.random((lo, n - lo)) < p
    perm = rng.permutation(n)
    out: BoolMatrix = rel[np.ix_(perm, perm)]
    return out
