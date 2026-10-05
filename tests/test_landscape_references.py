"""Pruebas de las referencias del paisaje L-2 (omega.landscape.references)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.experiments.reference_graphs import periodic_lattice
from omega.landscape import references as R
from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray

N = 216
PAIRS = N * (N - 1) // 2


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([20261005, seed])))


def _mass(w: FloatArray) -> float:
    return float(w.sum() / 2.0)


def test_torus_lattice_matches_reference_and_degree() -> None:
    w = R.torus_lattice_3d(6)
    assert w.shape == (N, N)
    assert np.array_equal(w, periodic_lattice((6, 6, 6)))
    validate_weight_matrix(w)
    assert np.all(w.sum(axis=1) == 6.0)
    assert _mass(w) == 648.0 == 6 / 215 * PAIRS


def test_decorated_lattice() -> None:
    w = R.decorated_lattice_3d(3, 8)
    validate_weight_matrix(w)
    assert w.shape == (216, 216)
    assert np.all(w.sum(axis=1) == 13.0)
    assert set(np.unique(w)) == {0.0, 1.0}
    # bloque K8 completo en el sitio 0 y union vertice a vertice
    assert np.all(w[:8, :8] == 1.0 - np.eye(8))
    sites = R.torus_lattice_3d(3)
    assert np.array_equal(np.kron(sites, np.eye(8)), w - np.kron(np.eye(27), 1.0 - np.eye(8)))
    # 6 vecinos distintos por sitio en lado 3
    assert np.all(sites.sum(axis=1) == 6.0)
    with pytest.raises(ValueError):
        R.decorated_lattice_3d(2, 8)


def test_rgg_binary_vs_points_and_determinism() -> None:
    k = 12
    w1 = R.rgg3_torus_binary(N, k, _rng(3))
    w2 = R.rgg3_torus_binary(N, k, _rng(3))
    assert np.array_equal(w1, w2)
    assert not np.array_equal(w1, R.rgg3_torus_binary(N, k, _rng(4)))
    validate_weight_matrix(w1)
    pts = R.rgg3_torus_points(N, _rng(3))
    step = R.weighted_rgg_from_points(pts, R.rgg3_radius(N, k), "step", 1.0)
    assert np.array_equal(step, w1)


@pytest.mark.parametrize("profile", R.PROFILES)
def test_weighted_rgg_mass_exact(profile: R.Profile) -> None:
    pts = R.rgg3_torus_points(N, _rng(0))
    r12 = R.rgg3_radius(N, 12)
    for total in (648.0, 1296.0):
        r = 2.0 * r12
        a = R.fit_amplitude_for_mass(pts, r, profile, total)
        w = R.weighted_rgg_from_points(pts, r, profile, a)
        validate_weight_matrix(w)
        assert abs(_mass(w) - total) < 1e-9
    # inalcanzable: muy pocos pares dentro de 0.5 r12 para masa 2322
    with pytest.raises(ValueError):
        R.fit_amplitude_for_mass(pts, 0.5 * r12, profile, 2322.0)


def test_weighted_rgg_support_and_profiles() -> None:
    pts = R.rgg3_torus_points(60, _rng(1))
    d = R.torus_distances(pts)
    r = 0.3
    for profile in R.PROFILES:
        w = R.weighted_rgg_from_points(pts, r, profile, 0.7)
        assert np.all(w[d >= r] == 0.0)
        assert w.max() <= 0.7 + 1e-15
    lin = R.weighted_rgg_from_points(pts, r, "linear", 50.0)
    assert lin.max() == 1.0


def test_colex_clique_masses_and_structure() -> None:
    for total in (648.0, 1296.0, 2322.0, 100.5, 0.0):
        w = R.colex_clique_with_mass(N, total)
        validate_weight_matrix(w)
        assert abs(_mass(w) - total) < 1e-9
    w = R.colex_clique_with_mass(N, 6.0)  # K4 sobre los nodos 0..3
    assert np.array_equal(w[:4, :4], 1.0 - np.eye(4))
    assert w[4:].sum() == 0.0
    w = R.colex_clique_with_mass(N, 7.5)
    assert w[0, 4] == 1.0 and w[1, 4] == 0.5 and w[2, 4] == 0.0


def test_clique_union_and_caveman() -> None:
    w = R.clique_union(N, 7)
    validate_weight_matrix(w)
    assert np.all(w[:7].sum(axis=1)[:7] == 6.0)
    assert w[N - 1].sum() == 12.0  # ultima absorbe 216 - 29*7 = 13 nodos
    w27 = R.clique_union(N, 8)
    assert np.all(w27.sum(axis=1) == 7.0) and _mass(w27) == 27 * 28
    cav = R.connected_caveman(N, 8)
    validate_weight_matrix(cav)
    import networkx as nx

    g = nx.from_numpy_array(cav)
    assert nx.is_connected(g)
    ref = nx.to_numpy_array(nx.connected_caveman_graph(27, 8), nodelist=range(N))
    assert np.array_equal(cav, ref)
    assert _mass(cav) == 27 * 28


def test_er_and_uniform() -> None:
    w = R.erdos_renyi_m(N, 1296, _rng(0))
    validate_weight_matrix(w)
    assert _mass(w) == 1296.0
    assert np.array_equal(w, R.erdos_renyi_m(N, 1296, _rng(0)))
    u = R.uniform_with_mass(N, 0.1 * PAIRS)
    validate_weight_matrix(u)
    assert abs(_mass(u) - 0.1 * PAIRS) < 1e-9
    assert np.allclose(u[np.triu_indices(N, 1)], 0.1)
    with pytest.raises(ValueError):
        R.uniform_with_mass(N, PAIRS + 5.0)


@pytest.mark.parametrize("total", [648.0, 1296.0, 2322.0, 1000.25])
def test_match_mass_binary_both_directions(total: float) -> None:
    for base in (R.clique_union(N, 7), R.decorated_lattice_3d(), R.torus_lattice_3d(6)):
        out = R.match_mass(base, total)
        assert out is not None
        validate_weight_matrix(out)
        assert abs(_mass(out) - total) < 1e-9
        frac = out[np.triu_indices(N, 1)]
        assert np.count_nonzero((frac > 0) & (frac < 1)) <= 1  # una sola arista parcial
        # no muta la entrada
    base = R.torus_lattice_3d(6)
    before = base.copy()
    R.match_mass(base, 100.0)
    assert np.array_equal(base, before)


def test_match_mass_not_matchable_and_noop() -> None:
    w = R.torus_lattice_3d(6)
    assert R.match_mass(w, PAIRS + 10.0) is None
    assert R.match_mass(w, -1.0) is None
    same = R.match_mass(w, 648.0)
    assert same is not None and np.array_equal(same, w)


def test_colex_order() -> None:
    i, j = R.colex_pairs(5)
    assert list(zip(i[:6].tolist(), j[:6].tolist())) == [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]
