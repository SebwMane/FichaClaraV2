"""Pruebas de omega.geometry.distance_suite (DESIGN §1.1, §3.2 WP-B)."""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
import pytest

from omega.config.settings import DimensionConfig, GraphConfig, SpectralConfig
from omega.config.settings11 import GEODESIC_MODES, DistanceMode, DistanceSuiteConfig
from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.geometry.distance_suite import (
    distance_suite,
    edge_length_matrix_mode,
    metric_spread,
    mode_dimension,
    mode_distance_matrix,
    resistance_exponent,
    resistance_matrix,
)
from omega.geometry.distances import distance_matrix, hop_distance_matrix, threshold_adjacency
from omega.geometry.observables import geometry_observables
from omega.types import DimensionEstimate, FloatArray

GRAPH = GraphConfig()
DIM = DimensionConfig()
SUITE = DistanceSuiteConfig()


def _weighted(shape: tuple[int, ...], seed: int, lo: float = 0.3) -> FloatArray:
    rng = np.random.Generator(np.random.PCG64(seed))
    u = np.triu(rng.uniform(lo, 1.0, size=(int(np.prod(shape)),) * 2), 1)
    return np.asarray(periodic_lattice(shape) * (u + u.T), dtype=np.float64)


@lru_cache(maxsize=None)
def _suite(name: str):  # type: ignore[no-untyped-def]
    shapes = {"ring": (300,), "t17": (17, 17), "t20": (20, 20), "t7": (7, 7, 7)}
    return distance_suite(periodic_lattice(shapes[name]), GRAPH, DIM, SUITE)


def _est(value: float, status: str = "ok") -> DimensionEstimate:
    z = np.zeros(0)
    return DimensionEstimate(value, 0.0, status, None, z, z, z, True, "t")  # type: ignore[arg-type]


def test_edge_lengths_modes() -> None:
    w = np.array([[0.0, 1.0, 0.5], [1.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    a = threshold_adjacency(w, 0.1)
    hop = edge_length_matrix_mode(w, a, DistanceMode.HOP, 1e-9, 1e-12).toarray()
    inv = edge_length_matrix_mode(w, a, DistanceMode.WEIGHTED_INVERSE, 1e-9, 1e-12).toarray()
    log = edge_length_matrix_mode(w, a, DistanceMode.WEIGHTED_LOG, 1e-9, 1e-12).toarray()
    assert hop[0, 1] == 1.0 and hop[0, 2] == 1.0 and hop[1, 2] == 0.0
    assert inv[0, 2] == pytest.approx(1.0 / (0.5 + 1e-9))
    assert log[0, 1] == 1.0  # W=1 da longitud 1 (sin colapso)
    assert log[0, 2] == pytest.approx(1.0 - math.log(0.5))
    with pytest.raises(ValueError):
        edge_length_matrix_mode(w, a, DistanceMode.RESISTANCE, 1e-9, 1e-12)


def test_strict_threshold() -> None:
    w = np.array([[0.0, 0.1], [0.1, 0.0]])  # W == w_min: sin arista (estricto)
    d = mode_distance_matrix(w, 0.1, DistanceMode.WEIGHTED_LOG, GRAPH, SUITE)
    assert math.isinf(d[0, 1])


def test_lattices_binary_modes_identical() -> None:
    for shape in [(40,), (9, 9), (5, 5, 5)]:
        w = periodic_lattice(shape)
        dh = mode_distance_matrix(w, GRAPH.w_min, DistanceMode.HOP, GRAPH, SUITE)
        dl = mode_distance_matrix(w, GRAPH.w_min, DistanceMode.WEIGHTED_LOG, GRAPH, SUITE)
        assert np.array_equal(dh, dl)
        di = mode_distance_matrix(w, GRAPH.w_min, DistanceMode.WEIGHTED_INVERSE, GRAPH, SUITE)
        assert np.allclose(dh, di, rtol=1e-8)
    res = distance_suite(periodic_lattice((15, 15)), GRAPH, DIM, SUITE)
    vals = [res.estimates[m] for m in GEODESIC_MODES]
    assert all(e.status == "ok" for e in vals)
    assert vals[0].value == vals[1].value == vals[2].value
    assert res.metric_spread == 0.0 and not res.distance_sensitive


def test_ring_resistance_and_zeta() -> None:
    r = _suite("ring")
    e = r.estimates[DistanceMode.RESISTANCE]
    assert e.status == "ok" and e.value == pytest.approx(1.0, abs=0.05)
    assert r.resistance_exponent > 0.8


def test_resistance_no_window_in_2d_3d_and_zeta() -> None:
    for name in ("t17", "t7"):
        r = _suite(name)
        assert r.estimates[DistanceMode.RESISTANCE].status == "no_window"
        assert r.fallback_used[DistanceMode.RESISTANCE] is True
    assert 0.35 <= _suite("t20").resistance_exponent <= 0.65
    assert _suite("t7").resistance_exponent < 0.35


def test_heavy_9cube_not_distance_sensitive_and_matches_omega10() -> None:
    w = _weighted((9, 9, 9), seed=3)
    suite = SUITE.__class__(modes=GEODESIC_MODES)
    r = distance_suite(w, GRAPH, DIM, suite)
    assert r.metric_spread <= 0.5 and not r.distance_sensitive
    obs = geometry_observables(w, GRAPH, DIM, SpectralConfig())
    inv = r.estimates[DistanceMode.WEIGHTED_INVERSE]
    assert obs.d_eff.status == "ok" == inv.status
    assert inv.value == obs.d_eff.value and inv.window == obs.d_eff.window
    hop = r.estimates[DistanceMode.HOP]
    assert hop.value == obs.d_eff_hops.value and hop.window == obs.d_eff_hops.window
    assert np.array_equal(
        mode_distance_matrix(w, GRAPH.w_min, DistanceMode.WEIGHTED_INVERSE, GRAPH, SUITE),
        distance_matrix(w, GRAPH.w_min, GRAPH.epsilon),
    )
    a = threshold_adjacency(w, GRAPH.w_min)
    assert np.array_equal(
        mode_distance_matrix(w, GRAPH.w_min, DistanceMode.HOP, GRAPH, SUITE), hop_distance_matrix(a)
    )


def test_distance_sensitive_when_modes_missing_window() -> None:
    est = {DistanceMode.HOP: _est(2.0), DistanceMode.WEIGHTED_INVERSE: _est(2.0), DistanceMode.WEIGHTED_LOG: _est(math.nan, "no_window")}
    assert metric_spread(est) == math.inf
    # solo INV con ventana
    only_inv = {DistanceMode.HOP: _est(math.nan, "no_window"), DistanceMode.WEIGHTED_INVERSE: _est(2.0), DistanceMode.WEIGHTED_LOG: _est(math.nan, "no_window")}
    assert metric_spread(only_inv) == math.inf
    ok3 = {m: _est(v) for m, v in zip(GEODESIC_MODES, (2.0, 2.3, 2.5), strict=True)}
    assert metric_spread(ok3) == pytest.approx(0.5)


def test_sparse_hop_window_is_distance_sensitive() -> None:
    rng = np.random.Generator(np.random.PCG64(0))
    pts = rng.random((120, 2))
    d = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    w = np.exp(-d / 0.05)
    np.fill_diagonal(w, 0.0)
    r = distance_suite(w, GraphConfig(w_min=0.0), DIM, SuiteGeo)
    assert r.estimates[DistanceMode.HOP].status == "no_window"
    assert r.distance_sensitive and r.metric_spread == math.inf


SuiteGeo = DistanceSuiteConfig(modes=GEODESIC_MODES)


def test_complete_graph_no_window() -> None:
    r = distance_suite(complete_graph(30), GRAPH, DIM, SUITE)
    assert all(e.status == "no_window" for e in r.estimates.values())
    assert r.distance_sensitive and math.isnan(r.resistance_exponent)


def test_trivial_giant() -> None:
    r = distance_suite(np.zeros((5, 5)), GRAPH, DIM, SUITE)
    assert all(e.status == "insufficient_component" for e in r.estimates.values())
    assert r.distance_sensitive


def test_resistance_triangle_inequality_and_known_values() -> None:
    rng = np.random.Generator(np.random.PCG64(5))
    n = 25
    u = np.triu(rng.uniform(0.2, 1.0, (n, n)) * (rng.random((n, n)) < 0.3), 1)
    w = u + u.T
    a = threshold_adjacency(w, 0.1)
    r = resistance_matrix(w, a)
    fin = np.isfinite(r)
    assert np.all(np.diag(r) == 0.0) and np.allclose(r, r.T)
    for i in range(n):
        for j in range(n):
            if not fin[i, j]:
                continue
            bound = r[i, :] + r[:, j]
            ok = np.isfinite(bound)
            assert np.all(r[i, j] <= bound[ok] + 1e-9)
    # cadena de 3 nodos con W=1: R_02 = 2; dos aristas en paralelo (triangulo): R=2/3
    c = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.float64)
    assert resistance_matrix(c, c > 0)[0, 2] == pytest.approx(2.0)
    t = complete_graph(3)
    assert resistance_matrix(t, t > 0)[0, 1] == pytest.approx(2.0 / 3.0)


def test_resistance_components_and_limits() -> None:
    w = np.zeros((5, 5))
    w[0, 1] = w[1, 0] = 1.0
    w[2, 3] = w[3, 2] = 0.5
    r = resistance_matrix(w, w > 0)
    assert r[0, 1] == pytest.approx(1.0) and r[2, 3] == pytest.approx(2.0)
    assert math.isinf(r[0, 2]) and math.isinf(r[4, 0]) and r[4, 4] == 0.0
    with pytest.raises(ValueError):
        resistance_matrix(complete_graph(10), complete_graph(10) > 0, max_nodes=5)


def test_resistance_exponent_nan_and_ring_value() -> None:
    w = periodic_lattice((60,))
    hop = hop_distance_matrix(w > 0)
    r = resistance_matrix(w, w > 0)
    assert math.isnan(resistance_exponent(r, hop, 1.0))
    assert 0.8 < resistance_exponent(r, hop, 10.0) < 1.05


def test_mode_dimension_fallback_flag() -> None:
    k = complete_graph(30)
    d = hop_distance_matrix(k > 0)
    est, fb = mode_dimension(d, DIM, SUITE, integer_metric=True)
    assert est.status == "no_window" and fb is False
    w = _weighted((7, 7, 7), seed=1)
    r = mode_distance_matrix(w, 0.1, DistanceMode.RESISTANCE, GRAPH, SUITE)
    est, fb = mode_dimension(r, DIM, SUITE, integer_metric=False)
    assert est.status == "no_window" and fb is True


def test_resistance_consistent() -> None:
    from omega.certificate.taxonomy import resistance_consistent

    s = SUITE
    assert resistance_consistent(0.9, 1, s) and not resistance_consistent(0.5, 1, s)
    assert resistance_consistent(0.35, 2, s) and resistance_consistent(0.65, 2, s)
    assert not resistance_consistent(0.7, 2, s)
    assert resistance_consistent(0.3, 3, s) and not resistance_consistent(0.35, 3, s)
    assert not resistance_consistent(0.5, None, s) and not resistance_consistent(math.nan, 2, s)
