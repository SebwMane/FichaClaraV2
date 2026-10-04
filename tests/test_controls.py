"""Pruebas de controles nulos (WP-E): invariantes, determinismo y casos límite."""

from __future__ import annotations

import networkx as nx
import numpy as np
import pytest

from omega.config.seeds import SeedKey as _SK, make_rng as _mk

from omega.config.seeds import SeedKey
from omega.config.settings11 import NullModel
from omega.controls.battery import make_null, null_battery, null_seed_key
from omega.controls.configuration_model import configuration_model, matched_configuration_model
from omega.controls.degree_preserving_rewire import degree_preserving_rewire
from omega.controls.erdos_renyi import erdos_renyi_gnm, erdos_renyi_gnp, matched_erdos_renyi
from omega.controls.random_geometric import balanced_tree, cap_angle, rgg_sphere, rgg_torus
from omega.controls.small_world import matched_small_world, ring_lattice, watts_strogatz


def _rng(seed: int) -> np.random.Generator:
    return _mk(_SK(12345, (seed,)))


def _weighted(n: int, seed: int, p: float = 0.15) -> np.ndarray:
    rng = _rng(seed)
    iu = np.triu_indices(n, 1)
    v = rng.random(iu[0].size) * (rng.random(iu[0].size) < p)
    w = np.zeros((n, n))
    w[iu] = v
    out = np.asarray(w + w.T, dtype=np.float64)
    return out


def _deg(w: np.ndarray, w_min: float = 0.0) -> np.ndarray:
    return np.asarray((w > w_min).sum(axis=1))


def test_gnm_exact_edges_and_symmetry() -> None:
    a = erdos_renyi_gnm(40, 123, _rng(0))
    assert a.sum() == 2 * 123
    assert np.array_equal(a, a.T) and np.all(np.diag(a) == 0)


def test_gnp_density() -> None:
    n, p = 200, 0.1
    a = erdos_renyi_gnp(n, p, _rng(1))
    dens = a.sum() / (n * (n - 1))
    assert abs(dens - p) < 0.01


def test_matched_er_binary_and_shuffled() -> None:
    w = _weighted(40, 2)
    m = int((np.triu(w, 1) > 0.05).sum())
    b = matched_erdos_renyi(w, 0.05, _rng(3))
    assert b.sum() == 2 * m
    s = matched_erdos_renyi(w, 0.05, _rng(3), "shuffled")
    assert int((np.triu(s, 1) > 0).sum()) == m
    assert np.allclose(np.sort(np.triu(s, 1)[np.triu(s, 1) > 0]), np.sort(np.triu(w, 1)[np.triu(w, 1) > 0.05]))
    with pytest.raises(ValueError):
        erdos_renyi_gnm(5, 11, _rng(0))


def test_rewire_preserves_degrees_and_weights() -> None:
    w = _weighted(60, 4, 0.2)
    r = degree_preserving_rewire(w, 0.0, _rng(5))
    assert np.array_equal(_deg(r), _deg(w))
    assert np.array_equal(r, r.T)
    assert np.allclose(np.sort(np.triu(r, 1).ravel()), np.sort(np.triu(w, 1).ravel()))
    assert not np.array_equal(r, w)


def test_configuration_model_degrees() -> None:
    w = _weighted(120, 6, 0.1)
    deg = _deg(w)
    c = matched_configuration_model(w, 0.0, _rng(7))
    dc = _deg(c)
    assert np.all(dc <= deg)
    assert dc.sum() >= 0.95 * deg.sum()
    assert np.all(np.diag(c) == 0) and np.array_equal(c, c.T)
    with pytest.raises(ValueError):
        configuration_model(np.array([1, 1, 1]), _rng(0))


def test_ws_p0_is_ring_lattice() -> None:
    n, k = 30, 4
    ws = watts_strogatz(n, k, 0.0, _rng(0))
    ring = ring_lattice(n, k)
    assert np.array_equal(ws, ring)
    assert np.all(ring.sum(axis=1) == k)
    g = nx.from_numpy_array(ring)
    assert nx.is_isomorphic(g, nx.circulant_graph(n, [1, 2]))


def test_ws_preserves_edge_count_and_matched() -> None:
    ws = watts_strogatz(60, 6, 0.3, _rng(1))
    assert ws.sum() == 60 * 6
    assert np.array_equal(ws, ws.T) and np.all(np.diag(ws) == 0)
    w = _weighted(60, 8, 0.15)
    kbar = (w > 0).sum() / 60
    m = matched_small_world(w, 0.0, 0.1, _rng(2))
    assert abs(m.sum() / 60 - kbar) <= 1.0 + 1e-9


def test_rgg_torus_mean_degree() -> None:
    n, k = 400, 8.0
    a = rgg_torus(n, 2, k, _rng(0))
    assert abs(a.sum() / n - k) < 0.15 * k
    e = rgg_torus(n, 2, k, _rng(0), "euclidean")
    assert np.array_equal(e > 0, a > 0) and e.max() <= 1.0
    with pytest.raises(ValueError):
        rgg_torus(10, 2, 9.0, _rng(0))


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_rgg_sphere_mean_degree(dim: int) -> None:
    n, k = 500, 10.0
    a = rgg_sphere(n, dim, k, _rng(dim))
    assert abs(a.sum() / n - k) < 0.1 * k


def test_cap_angle_analytic() -> None:
    n, k = 1001, 10.0
    r2 = cap_angle(n, 2, k)
    assert abs((1 - np.cos(r2)) / 2 - k / (n - 1)) < 1e-10
    r3 = cap_angle(n, 3, k)
    assert abs((r3 - np.sin(r3) * np.cos(r3)) / np.pi - k / (n - 1)) < 1e-10


def test_balanced_tree() -> None:
    a = balanced_tree(2, 3)
    assert a.shape == (15, 15) and a.sum() == 2 * 14


def test_determinism_by_seed() -> None:
    w = _weighted(40, 9)
    for kind in NullModel:
        a = make_null(w, 0.0, kind, _rng(11))
        b = make_null(w, 0.0, kind, _rng(11))
        assert np.array_equal(a, b), kind


def test_seed_keys_distinct_and_battery() -> None:
    base = SeedKey(7, (1, 2, 3))
    keys = {null_seed_key(base, k) for k in NullModel}
    assert len(keys) == len(NullModel)
    assert null_seed_key(base, NullModel.ERDOS_RENYI).spawn_key == (1, 2, 3, 1000)
    w = _weighted(40, 10)
    kinds = [NullModel.DEGREE_PRESERVING_REWIRE, NullModel.CONFIGURATION_MODEL, NullModel.ERDOS_RENYI,
             NullModel.SHUFFLED_WEIGHTS]
    b1 = null_battery(w, 0.0, kinds, base)
    b2 = null_battery(w, 0.0, kinds, base)
    assert list(b1) == kinds
    for k in kinds:
        assert np.array_equal(b1[k], b2[k])
