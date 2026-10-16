"""Referencias de Omega-C0 (prerregistro §3), todas N = 216 por defecto; adyacencias 0/1 simetricas float64.

La amplitud t en (0,1] se aplica despues (`optimal_amplitude`). Aleatoriedad solo con PCG64 explicito.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import networkx as nx
import numpy as np
from scipy.optimize import minimize_scalar

from omega.c0.functional import C0Params
from omega.landscape.references import clique_union, erdos_renyi_m, rgg3_torus_binary, torus_lattice_3d
from omega.types import FloatArray

__all__ = [
    "rng_from_key",
    "torus_lattice_3d",
    "rgg3_torus_binary",
    "clique_union",
    "erdos_renyi_m",
    "random_regular",
    "complete_bipartite_half",
    "random_tree",
    "ring_lattice",
    "watts_strogatz",
    "triangular_torus",
    "uniform_complete",
    "AmpProfile",
    "amplitude_profile",
    "optimal_amplitude_from_profile",
    "optimal_amplitude",
]

AMP_GRID = 2000


def rng_from_key(key: Sequence[int]) -> np.random.Generator:
    """PCG64 determinista a partir de una clave entera (SeedSequence)."""
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([int(k) for k in key])))


def _nx_seed(rng: np.random.Generator) -> int:
    return int(rng.integers(2**31))


def _to_adj(g: nx.Graph, n: int) -> FloatArray:
    return np.asarray(nx.to_numpy_array(g, nodelist=range(n), dtype=np.float64), dtype=np.float64)


def random_regular(n: int, k: int, rng: np.random.Generator) -> FloatArray:
    """Grafo k-regular aleatorio (n*k par)."""
    return _to_adj(nx.random_regular_graph(k, n, seed=_nx_seed(rng)), n)


def complete_bipartite_half(n: int = 216) -> FloatArray:
    """K_{n/2,n/2}: primera mitad contra segunda mitad."""
    if n % 2:
        raise ValueError("n debe ser par")
    h = n // 2
    a = np.zeros((n, n), dtype=np.float64)
    a[:h, h:] = 1.0
    a[h:, :h] = 1.0
    return a


def random_tree(n: int, rng: np.random.Generator) -> FloatArray:
    """Arbol etiquetado uniforme (secuencia de Prufer con el rng)."""
    if n < 3:
        raise ValueError("n >= 3")
    seq = [int(x) for x in rng.integers(0, n, size=n - 2)]
    return _to_adj(nx.from_prufer_sequence(seq), n)


def ring_lattice(n: int = 216, k: int = 6) -> FloatArray:
    """Anillo: cada nodo unido a k/2 vecinos por lado (k par)."""
    if k % 2 or k < 2 or k >= n:
        raise ValueError("k debe ser par y < n")
    a = np.zeros((n, n), dtype=np.float64)
    idx = np.arange(n)
    for d in range(1, k // 2 + 1):
        a[idx, (idx + d) % n] = 1.0
        a[(idx + d) % n, idx] = 1.0
    return a


def watts_strogatz(n: int, k: int, beta: float, rng: np.random.Generator) -> FloatArray:
    """Watts-Strogatz (anillo k, recableado beta)."""
    return _to_adj(nx.watts_strogatz_graph(n, k, beta, seed=_nx_seed(rng)), n)


def triangular_torus(rows: int = 12, cols: int = 18) -> FloatArray:
    """Reticula triangular 6-regular en toro rows x cols (vecinos: +-(0,1), +-(1,0), +-(1,1))."""
    n = rows * cols
    a = np.zeros((n, n), dtype=np.float64)
    for r in range(rows):
        for c in range(cols):
            i = r * cols + c
            for dr, dc in ((0, 1), (1, 0), (1, 1)):
                j = ((r + dr) % rows) * cols + (c + dc) % cols
                a[i, j] = 1.0
                a[j, i] = 1.0
    return a


def uniform_complete(n: int = 216) -> FloatArray:
    """Todo-uno fuera de la diagonal."""
    a = np.ones((n, n), dtype=np.float64)
    np.fill_diagonal(a, 0.0)
    return a


# ----------------------------------------------------------- amplitud optima


@dataclass(frozen=True)
class AmpProfile:
    """Resumen de t*A: valores unicos de C0 = A^2 sobre aristas (con multiplicidad), A y sum k^2."""

    a_vals: FloatArray
    c_vals: FloatArray
    counts: FloatArray
    sum_k2: float


def amplitude_profile(adj: FloatArray) -> AmpProfile:
    """Compacta una adyacencia 0/1 (o con pesos fijos) para evaluar S(t A) rapido."""
    iu = np.triu_indices(adj.shape[0], k=1)
    c0 = (adj @ adj)[iu]
    av = adj[iu]
    m = av > 0.0
    key = np.stack([np.round(av[m], 9), np.round(c0[m], 9)], axis=1)
    uniq, cnt = np.unique(key, axis=0, return_counts=True)
    k = adj.sum(axis=1)
    return AmpProfile(
        a_vals=np.asarray(uniq[:, 0], dtype=np.float64),
        c_vals=np.asarray(uniq[:, 1], dtype=np.float64),
        counts=np.asarray(cnt, dtype=np.float64),
        sum_k2=float(np.dot(k, k)),
    )


def _s_of_t(t: FloatArray, prof: AmpProfile, p: C0Params) -> FloatArray:
    tt = t[:, None]
    c = tt**2 * prof.c_vals[None, :]
    psi_v = p.a + p.b * np.log1p(c) - p.lam * c
    first = -np.sum(prof.counts[None, :] * tt * prof.a_vals[None, :] * psi_v, axis=1)
    return np.asarray(first + p.kappa * t**2 * prof.sum_k2, dtype=np.float64)


def optimal_amplitude_from_profile(prof: AmpProfile, p: C0Params) -> tuple[float, float]:
    """min_t S(tA), t en (0,1]: rejilla de 2000 valores y refinamiento acotado."""
    grid = np.arange(1, AMP_GRID + 1, dtype=np.float64) / AMP_GRID
    s = _s_of_t(grid, prof, p)
    i = int(np.argmin(s))
    best_t, best_s = float(grid[i]), float(s[i])
    lo = float(grid[max(i - 1, 0)]) if i > 0 else 1e-9
    hi = float(grid[min(i + 1, AMP_GRID - 1)])
    res = minimize_scalar(
        lambda x: float(_s_of_t(np.array([x]), prof, p)[0]),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": 1e-12},
    )
    if float(res.fun) < best_s:
        best_t, best_s = float(res.x), float(res.fun)
    return best_t, best_s


def optimal_amplitude(adj: FloatArray, p: C0Params) -> tuple[float, float]:
    """(t, S) minimizando action_c0(t*adj) con t en (0,1]."""
    return optimal_amplitude_from_profile(amplitude_profile(adj), p)
