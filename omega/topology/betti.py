"""Numeros de Betti sobre Z/2 (OMEGA_1_1_DESIGN §1.2, §3.3).

* `clique_complex_betti`: complejo de cliques truncado; beta_k exactos hasta `max_dim`.
* `short_cycle_betti1`: beta1^(L) = E - rank(d1) - rank(caras), con caras = triangulos y 4-ciclos
  rellenados. En retículas cubicas el complejo de cliques no tiene triangulos (T3 9^3 da beta1=1459);
  el relleno de ciclos cortos da beta1(T3)=3, beta1(T2)=2 y anillo 1 (hallazgo 7).

El rango sobre GF(2) usa mascaras de bits en enteros de Python (convertir a `int` antes de
desplazar: `np.int64` desborda). Referencias: Hatcher, *Algebraic Topology*; Edelsbrunner-Harer,
*Computational Topology*.
"""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations

import numpy as np

from omega.contracts import BettiResult, ShortCycleBetti
from omega.types import BoolArray, IntArray

__all__ = [
    "gf2_rank",
    "count_triangles",
    "count_four_cycles",
    "clique_complex_betti",
    "short_cycle_betti1",
]


def _check_adjacency(a: BoolArray) -> None:
    if not isinstance(a, np.ndarray) or a.dtype != np.bool_:
        raise TypeError("a debe ser ndarray bool")
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"a debe ser cuadrada, forma {a.shape}")
    if not np.array_equal(a, a.T) or np.any(np.diag(a)):
        raise ValueError("a debe ser simetrica con diagonal False")


def gf2_rank(columns: Sequence[int]) -> int:
    """Rango sobre GF(2) de columnas dadas como enteros (mascaras de bits).

    Reduccion gaussiana con base indexada por el bit mas significativo.
    """
    basis: dict[int, int] = {}
    rank = 0
    for col in columns:
        v = int(col)
        while v:
            h = v.bit_length()
            b = basis.get(h)
            if b is None:
                basis[h] = v
                rank += 1
                break
            v ^= b
    return rank


def count_triangles(a: BoolArray) -> int:
    """Numero de triangulos T = tr(A^3)/6."""
    _check_adjacency(a)
    af = a.astype(np.float64)
    a2 = af @ af
    return int(round(float(np.sum(a2 * af)) / 6.0))


def count_four_cycles(a: BoolArray) -> int:
    """Numero de 4-ciclos (subgrafos) C4 = (tr(A^4) - 2 sum k_i^2 + 2m) / 8 (K4 da 3)."""
    _check_adjacency(a)
    af = a.astype(np.float64)
    a2 = af @ af
    tr4 = float(np.sum(a2 * a2))
    k = af.sum(axis=1)
    m = float(af.sum()) / 2.0
    return int(round((tr4 - 2.0 * float(np.sum(k * k)) + 2.0 * m) / 8.0))


def _components(a: BoolArray) -> int:
    from omega.topology.filtration import component_labels

    return component_labels(a)[0]


def clique_complex_betti(a: BoolArray, max_dim: int = 2, max_simplices: int = 300_000) -> BettiResult:
    """Betti (beta_0..beta_max_dim) del complejo de cliques sobre Z/2.

    Se enumeran simplices hasta dimension max_dim+1 (necesaria para beta_max_dim exacto).
    `counts` = (V, E, T, ...) hasta dimension max_dim+1; chi_beta = sum (-1)^k beta_k y
    chi_cuentas = sum (-1)^k counts_k (solo informe). Si #simplices de dimension >= 2
    supera `max_simplices` devuelve status "over_budget" (betti ceros, counts parciales).
    """
    _check_adjacency(a)
    if isinstance(max_dim, bool) or not isinstance(max_dim, int) or max_dim < 0:
        raise ValueError("max_dim debe ser entero >= 0")
    if isinstance(max_simplices, bool) or not isinstance(max_simplices, int) or max_simplices < 1:
        raise ValueError("max_simplices debe ser entero >= 1")
    n = a.shape[0]
    top = max_dim + 1
    nbrs_up = [np.flatnonzero(a[i, i + 1 :]) + i + 1 for i in range(n)]
    simplices: list[list[tuple[int, ...]]] = [[(i,) for i in range(n)]]
    big = 0
    for dim in range(1, top + 1):
        nxt: list[tuple[int, ...]] = []
        for s in simplices[-1]:
            cand = nbrs_up[s[-1]]
            for v in s[:-1]:
                if cand.size == 0:
                    break
                cand = cand[a[v, cand]]
            for c in cand:
                nxt.append((*s, int(c)))
                if dim >= 2:
                    big += 1
                    if big > max_simplices:
                        counts = tuple(len(x) for x in simplices) + (0,) * (top + 1 - len(simplices))
                        return BettiResult(
                            betti=(0,) * (max_dim + 1),
                            counts=counts,
                            euler_betti=0,
                            euler_counts=0,
                            status="over_budget",
                        )
        simplices.append(nxt)
        if not nxt:
            simplices.extend([] for _ in range(top - dim))
            break
    counts_t = tuple(len(x) for x in simplices)
    # rangos de d_k: simplices k -> (k-1); d_0 = 0
    ranks = [0] * (top + 2)
    for k in range(1, top + 1):
        if not simplices[k]:
            continue
        index = {s: i for i, s in enumerate(simplices[k - 1])}
        cols = []
        for s in simplices[k]:
            m = 0
            for r in range(len(s)):
                m |= 1 << index[s[:r] + s[r + 1 :]]
            cols.append(m)
        ranks[k] = gf2_rank(cols)
    betti = tuple(counts_t[k] - ranks[k] - ranks[k + 1] for k in range(max_dim + 1))
    euler_b = sum((-1) ** k * b for k, b in enumerate(betti))
    euler_c = sum((-1) ** k * c for k, c in enumerate(counts_t))
    return BettiResult(betti=betti, counts=counts_t, euler_betti=int(euler_b), euler_counts=int(euler_c), status="ok")


def _edge_index(a: BoolArray) -> tuple[IntArray, int]:
    iu, ju = np.nonzero(np.triu(a, k=1))
    n_e = int(iu.size)
    idx = np.full(a.shape, -1, dtype=np.int64)
    idx[iu, ju] = np.arange(n_e)
    idx[ju, iu] = np.arange(n_e)
    return idx, n_e


def short_cycle_betti1(a: BoolArray, max_length: int = 4, max_faces: int = 500_000) -> ShortCycleBetti:
    """beta1^(L) sobre Z/2 rellenando triangulos (L>=3) y 4-ciclos (L=4).

    beta1 = E - rank(d1) - rank(caras), con rank(d1) = N - beta0. `max_length` en {3, 4}.
    Si T + C4 > max_faces (conteos de forma cerrada) devuelve "over_budget" con b1=0 y n_faces
    = numero de caras contadas. b1_density = b1/E (0 si E=0).
    """
    _check_adjacency(a)
    if max_length not in (3, 4):
        raise ValueError("max_length debe ser 3 o 4")
    if isinstance(max_faces, bool) or not isinstance(max_faces, int) or max_faces < 1:
        raise ValueError("max_faces debe ser entero >= 1")
    n = a.shape[0]
    b0 = _components(a)
    idx, n_e = _edge_index(a)
    t = count_triangles(a)
    c4 = count_four_cycles(a) if max_length >= 4 else 0
    n_faces = t + c4
    if n_faces > max_faces:
        return ShortCycleBetti(b0=b0, b1=0, n_edges=n_e, n_faces=n_faces, b1_density=0.0,
                               max_length=max_length, status="over_budget")
    cols: list[int] = []
    nbrs = [np.flatnonzero(a[i]) for i in range(n)]
    for i in range(n):  # triangulos i<j<k
        up = nbrs[i][nbrs[i] > i]
        for j, k in combinations(up.tolist(), 2):
            if a[j, k]:
                cols.append((1 << int(idx[i, j])) | (1 << int(idx[i, k])) | (1 << int(idx[j, k])))
    if max_length >= 4:
        a2 = a.astype(np.int64) @ a.astype(np.int64)
        for i in range(n):  # i = nodo minimo del ciclo, j opuesto
            for j in range(i + 1, n):
                if a2[i, j] < 2:
                    continue
                cn = np.flatnonzero(a[i] & a[j])
                cn = cn[cn > i].tolist()
                for k, l in combinations(cn, 2):
                    cols.append(
                        (1 << int(idx[i, k])) | (1 << int(idx[k, j])) | (1 << int(idx[j, l])) | (1 << int(idx[l, i]))
                    )
    b1 = n_e - (n - b0) - gf2_rank(cols)
    return ShortCycleBetti(b0=b0, b1=b1, n_edges=n_e, n_faces=len(cols), b1_density=(b1 / n_e if n_e else 0.0),
                           max_length=max_length, status="ok")
