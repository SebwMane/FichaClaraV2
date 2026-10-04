"""Pruebas de la regla de coarse-graining heavy_edge_matching/max/v1 (WP-E)."""

from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from omega.config.seeds import SeedKey as _SK, make_rng as _mk

from omega.coarse_graining.rule import aggregate, coarse_grain, coarse_grain_hierarchy, heavy_edge_matching
from omega.config.settings11 import CoarseGrainConfig


def _rng(seed: int) -> np.random.Generator:
    return _mk(_SK(12345, (seed,)))


def _random_w(n: int, seed: int, p: float = 0.2, binary: bool = False) -> np.ndarray:
    rng = _rng(seed)
    iu = np.triu_indices(n, 1)
    mask = rng.random(iu[0].size) < p
    v = np.where(mask, 1.0 if binary else rng.random(iu[0].size), 0.0)
    w = np.zeros((n, n))
    w[iu] = v
    return w + w.T


def _ring(n: int) -> np.ndarray:
    w = np.zeros((n, n))
    i = np.arange(n)
    w[i, (i + 1) % n] = 1.0
    w[(i + 1) % n, i] = 1.0
    return w


def test_partition_and_maximal() -> None:
    w = _random_w(50, 0, 0.1, binary=True)
    blocks = heavy_edge_matching(w, 0.0, _rng(1))
    flat = sorted(i for b in blocks for i in b)
    assert flat == list(range(50))
    assert all(1 <= len(b) <= 2 and list(b) == sorted(b) for b in blocks)
    assert [b[0] for b in blocks] == sorted(b[0] for b in blocks)
    singles = [b[0] for b in blocks if len(b) == 1]
    assert not (w[np.ix_(singles, singles)] > 0).any()


def test_max_binary_equals_quotient_graph() -> None:
    w = _random_w(40, 2, 0.15, binary=True)
    blocks = heavy_edge_matching(w, 0.0, _rng(3))
    wp = aggregate(w, blocks, "max")
    g = nx.from_numpy_array(w)
    q = nx.quotient_graph(g, [set(b) for b in blocks], relabel=False)
    ref = np.zeros_like(wp)
    index = {frozenset(b): k for k, b in enumerate(blocks)}
    for a, b in q.edges():
        ref[index[a], index[b]] = ref[index[b], index[a]] = 1.0
    assert np.array_equal(wp, ref)


def test_mean_and_range() -> None:
    w = _random_w(30, 4, 0.3)
    blocks = heavy_edge_matching(w, 0.0, _rng(0))
    for how in ("max", "mean"):
        wp = aggregate(w, blocks, how)
        assert wp.min() >= 0 and wp.max() <= 1 and np.array_equal(wp, wp.T)
        assert np.all(np.diag(wp) == 0)
    with pytest.raises(ValueError):
        aggregate(w, [(0, 1)], "max")


def test_deterministic_by_seed() -> None:
    w = _random_w(60, 5, 0.1, binary=True)
    a = heavy_edge_matching(w, 0.0, _rng(9))
    b = heavy_edge_matching(w, 0.0, _rng(9))
    assert a == b


def test_ring_stays_ring_through_levels() -> None:
    cfg = CoarseGrainConfig()
    levels = coarse_grain_hierarchy(_ring(800), 0.0, cfg, _rng(0))
    assert 1 <= len(levels) <= 4
    ns = [800] + [lv.n for lv in levels]
    assert all(b < a and b >= 50 for a, b in zip(ns, ns[1:]))
    for lv in levels:
        a = lv.w > 0
        assert np.all(a.sum(axis=1) == 2)
        assert nx.is_connected(nx.from_numpy_array(lv.w))


def test_hierarchy_stops_below_n_min() -> None:
    levels = coarse_grain_hierarchy(_ring(80), 0.0, CoarseGrainConfig(), _rng(0))
    assert levels == ()  # 80 -> n' < 50 se descarta
    lv = coarse_grain(_ring(10), 0.0, CoarseGrainConfig(), _rng(0), level=3)
    assert 5 <= lv.n < 10 and lv.level == 3


def test_permutation_equivariance_distinct_weights() -> None:
    n = 40
    w = _random_w(n, 6, 0.2)  # pesos continuos distintos: sin empates
    rng = _rng(1)
    perm = rng.permutation(n)  # nodo nuevo p -> viejo perm[p]
    wp = w[np.ix_(perm, perm)]
    b0 = heavy_edge_matching(w, 0.0, _rng(10))
    b1 = heavy_edge_matching(wp, 0.0, _rng(99))  # otra semilla: u no influye
    inv = np.empty(n, dtype=int)
    inv[perm] = np.arange(n)
    mapped = {frozenset(int(inv[i]) for i in b) for b in b0}
    assert mapped == {frozenset(b) for b in b1}
    w0 = aggregate(w, b0, "max")
    w1 = aggregate(wp, b1, "max")
    # reordenar bloques de b0 según su imagen en b1
    pos = {frozenset(b): k for k, b in enumerate(b1)}
    order = np.array([pos[frozenset(int(inv[i]) for i in b)] for b in b0])
    assert np.array_equal(w0, w1[np.ix_(order, order)])
