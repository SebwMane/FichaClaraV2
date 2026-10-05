"""L-1 (Consejo Omega Rev. 2, R2.5): identidades y cotas exactas del funcional S0.

Alcance. Se prueban SOLO cinco afirmaciones: (1) S_smooth frente a su definicion directa;
(2) dS = -(alpha + 12 eta) dT bajo double-edge swaps que preservan grados; (3) cota de vecindario
T <= N C(k,2)/3 en grafos k-regulares (igualdad sii union disjunta de K_{k+1}), exhaustiva en el atlas
de <= 7 nodos; (4) cota de Lovasz T <= C(x,3) con C(x,2) = m; (5) tr W^3 <= (sum W^2)^{3/2} <= (sum W)^{3/2}
(sumas sobre la matriz completa, ambos triangulos). NO prueba nada sobre grados libres, pesos
optimos ni minimos locales.
"""

from __future__ import annotations

import math

import networkx as nx
import numpy as np
import pytest

from omega.config.settings import FunctionalParams
from omega.dynamics.functional import action, smoothness, triangles
from omega.network.initialization import random_uniform_weights
from omega.types import FloatArray

ATLAS = [g for g in nx.graph_atlas_g() if g.number_of_nodes() >= 1]


def _adj(g: nx.Graph) -> FloatArray:
    return np.asarray(nx.to_numpy_array(g, nodelist=sorted(g.nodes)), dtype=np.float64)


def _smooth_direct(w: FloatArray) -> float:
    """sum_{ij} W_ij sum_k (W_ik - W_jk)^2, definicion directa O(N^3)."""
    n = w.shape[0]
    total = 0.0
    for i in range(n):
        for j in range(n):
            if i != j and w[i, j] != 0.0:
                total += w[i, j] * float(np.sum((w[i, :] - w[j, :]) ** 2))
    return total


def _rel_close(a: float, b: float, rel: float = 1e-10) -> bool:
    return abs(a - b) <= rel * max(1.0, abs(a), abs(b))


def _random_weighted(n: int, seed: int) -> FloatArray:
    return random_uniform_weights(n, np.random.Generator(np.random.PCG64(seed)))


# ------------------------------------------------------------------ punto 1


@pytest.mark.parametrize("n", [3, 5, 9, 17])
def test_smooth_matches_direct_definition_weighted(n: int) -> None:
    for seed in range(5):
        w = _random_weighted(n, 100 * n + seed)
        k = w.sum(axis=1)
        s = np.diag(w @ w)
        assert _rel_close(smoothness(w), _smooth_direct(w))
        assert _rel_close(smoothness(w), 2.0 * float(np.sum(k * s)) - 2.0 * float(np.trace(w @ w @ w)))


def test_smooth_binary_equals_2sumk2_minus_12T_on_atlas() -> None:
    for g in ATLAS:
        if g.number_of_nodes() < 2:
            continue
        a = _adj(g)
        k = a.sum(axis=1)
        t = sum(nx.triangles(g).values()) / 3.0
        expected = 2.0 * float(np.sum(k * k)) - 12.0 * t
        assert _rel_close(smoothness(a), expected)
        assert _rel_close(smoothness(a), _smooth_direct(a))
        assert _rel_close(triangles(a), t)


# ------------------------------------------------------------------ punto 2


def _degree_preserving_swaps(w: FloatArray, n_swaps: int, seed: int) -> list[FloatArray]:
    rng = np.random.Generator(np.random.PCG64(seed))
    cur = w.copy()
    n = cur.shape[0]
    out: list[FloatArray] = []
    attempts = 0
    while len(out) < n_swaps and attempts < 20_000:
        attempts += 1
        edges = np.argwhere(np.triu(cur, 1) > 0)
        (a, b), (c, d) = edges[rng.integers(len(edges), size=2)]
        if len({a, b, c, d}) < 4 or cur[a, d] > 0 or cur[c, b] > 0:
            continue
        nxt = cur.copy()
        for x, y, v in ((a, b, 0.0), (c, d, 0.0), (a, d, 1.0), (c, b, 1.0)):
            nxt[x, y] = nxt[y, x] = v
        assert np.array_equal(nxt.sum(axis=1), cur.sum(axis=1))
        out.append(nxt)
        cur = nxt
    assert len(out) == n_swaps and n == cur.shape[0]
    return out


@pytest.mark.parametrize("eta", [0.05, 0.3, 1.7])
def test_degree_preserving_swaps_delta_action(eta: float) -> None:
    g = nx.gnp_random_graph(24, 0.3, seed=11)
    w0 = _adj(g)
    p = FunctionalParams(alpha=0.37, beta=1.3, gamma=0.21, eta=eta)
    prev = w0
    for nxt in _degree_preserving_swaps(w0, 60, seed=5):
        d_s = action(nxt, p) - action(prev, p)
        d_t = triangles(nxt) - triangles(prev)
        assert abs(d_s + (p.alpha + 12.0 * p.eta) * d_t) <= 1e-9 * max(1.0, abs(d_s))
        prev = nxt


# ------------------------------------------------------------------ puntos 3 y 4


def _is_union_of_complete(g: nx.Graph, size: int) -> bool:
    return all(len(c) == size and g.subgraph(c).number_of_edges() == size * (size - 1) // 2
               for c in nx.connected_components(g))


def test_neighbourhood_bound_k_regular_atlas() -> None:
    n_regular = n_equal = 0
    for g in ATLAS:
        n = g.number_of_nodes()
        degs = {d for _, d in g.degree}
        if len(degs) != 1:
            continue
        k = degs.pop()
        n_regular += 1
        t = sum(nx.triangles(g).values()) // 3
        bound = n * k * (k - 1) / 2 / 3
        assert t <= bound + 1e-12, (n, k, t)
        equal = abs(t - bound) < 1e-12
        n_equal += equal
        assert equal == _is_union_of_complete(g, k + 1), (n, k, t)
    assert n_regular > 20 and n_equal > 5


def test_lovasz_bound_atlas() -> None:
    for g in ATLAS:
        m = g.number_of_edges()
        t = sum(nx.triangles(g).values()) // 3
        x = (1.0 + math.sqrt(1.0 + 8.0 * m)) / 2.0  # C(x,2) = m
        assert t <= x * (x - 1.0) * (x - 2.0) / 6.0 + 1e-9, (g.number_of_nodes(), m, t)


# ------------------------------------------------------------------ punto 5


def _chain(w: FloatArray) -> tuple[float, float, float]:
    tr3 = float(np.trace(w @ w @ w))
    s2 = float(np.sum(w * w))
    s1 = float(np.sum(w))
    return tr3, s2 ** 1.5, s1 ** 1.5


def test_trace_cube_chain_random_and_extreme() -> None:
    rng = np.random.Generator(np.random.PCG64(20261005))
    mats: list[FloatArray] = []
    for _ in range(200):
        n = int(rng.integers(3, 41))
        mats.append(random_uniform_weights(n, rng))
    n = 12
    star = np.zeros((n, n))
    star[0, 1:] = star[1:, 0] = 1.0
    single = np.zeros((n, n))
    single[0, 1] = single[1, 0] = 1.0
    mats += [np.zeros((n, n)), np.ones((n, n)) - np.eye(n), star, single, _adj(nx.cycle_graph(n)),
             _adj(nx.complete_bipartite_graph(5, 7)), np.full((n, n), 0.5) - 0.5 * np.eye(n)]
    for w in mats:
        a, b, c = _chain(w)
        assert a <= b * (1.0 + 1e-12) + 1e-12
        assert b <= c * (1.0 + 1e-12) + 1e-12
    # el extremo K_n satura la segunda desigualdad solo en el caso 0/1 (W^2 = W): sum W^2 = sum W
    _, b, c = _chain(np.ones((n, n)) - np.eye(n))
    assert _rel_close(b, c)
