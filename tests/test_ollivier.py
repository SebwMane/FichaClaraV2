"""Tests de omega.curvature.ollivier (WP-C)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

from omega.config.settings import GraphConfig
from omega.config.settings11 import CurvatureConfig, DistanceMode, DistanceSuiteConfig
from omega.curvature.ollivier import neighbor_measure, ollivier_curvature, wasserstein1
from omega.experiments.reference_graphs import complete_graph, periodic_lattice, random_geometric_torus

G = GraphConfig()
S = DistanceSuiteConfig()


def _kappa(w: np.ndarray, cfg: CurvatureConfig | None = None, seed: int = 0):  # type: ignore[no-untyped-def]
    return ollivier_curvature(w, 0.1, G, S, cfg or CurvatureConfig(), np.random.Generator(np.random.PCG64(seed)))


def test_wasserstein_matches_assignment() -> None:
    rng = np.random.Generator(np.random.PCG64(0))
    for n in (3, 6):
        cost = rng.random((n, n))
        mu = np.full(n, 1.0 / n)
        r, c = linear_sum_assignment(cost)
        assert wasserstein1(mu, mu, cost) == pytest.approx(cost[r, c].sum() / n, abs=1e-9)


def test_wasserstein_validation() -> None:
    with pytest.raises(ValueError):
        wasserstein1(np.array([1.0]), np.array([0.5]), np.zeros((1, 1)))
    with pytest.raises(ValueError):
        wasserstein1(np.array([1.0]), np.array([1.0]), np.zeros((2, 2)))


def test_neighbor_measure() -> None:
    w = np.zeros((3, 3))
    w[0, 1] = w[1, 0] = 0.3
    w[0, 2] = w[2, 0] = 0.1
    s, m = neighbor_measure(w, w > 0, 0, 0.5)
    assert s.tolist() == [0, 1, 2] and np.allclose(m, [0.5, 0.375, 0.125])
    s, m = neighbor_measure(w, np.zeros((3, 3), dtype=np.bool_), 1, 0.5)
    assert s.tolist() == [1] and m.tolist() == [1.0]


@pytest.mark.parametrize("shape", [(6, 6), (6, 6, 6)])
def test_hypercubic_tori_flat(shape: tuple[int, ...]) -> None:
    s = _kappa(periodic_lattice(shape), CurvatureConfig(max_edges=60))
    assert np.max(np.abs(s.edge_values)) <= 1e-9 and s.ok and s.sampled


def test_small_torus_measured_deviation() -> None:
    # Medido: T3 4^3 da kappa = 1/6 (lazos de envoltura de longitud 4), no 0; T2 4x4 da 1/4.
    s = _kappa(periodic_lattice((4, 4, 4)), CurvatureConfig(max_edges=40))
    assert s.mean == pytest.approx(1 / 6, abs=1e-9)


def test_cycle_flat() -> None:
    s = _kappa(periodic_lattice((8,)))
    assert np.max(np.abs(s.edge_values)) <= 1e-9 and s.n_edges == 8 and not s.sampled


@pytest.mark.parametrize("n", [5, 10])
def test_complete_graph_closed_form(n: int) -> None:
    s = _kappa(complete_graph(n))
    expected = 1.0 - abs(0.5 - 0.5 / (n - 1))
    assert np.max(np.abs(s.edge_values - expected)) <= 1e-9 and s.mean > 0


def test_inverse_mode_and_empty() -> None:
    cfg = CurvatureConfig(distance_mode=DistanceMode.WEIGHTED_INVERSE)
    s = _kappa(periodic_lattice((8,)), cfg)
    assert abs(s.mean) <= 1e-9
    s0 = _kappa(np.zeros((5, 5)))
    assert s0.n_edges == 0 and not s0.ok


def test_sampling_reproducible() -> None:
    w = random_geometric_torus(60, 2, 8, np.random.Generator(np.random.PCG64(0)), euclidean=False)
    cfg = CurvatureConfig(max_edges=15)
    a, b = _kappa(w, cfg, 3), _kappa(w, cfg, 3)
    assert a.sampled and a.n_edges == 15 and np.array_equal(a.edge_values, b.edge_values)


def _sphere_rgg(n: int, k: float, rng: np.random.Generator) -> np.ndarray:
    p = rng.normal(size=(n, 3))
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    cos_t = 1.0 - 2.0 * k / n  # casquete con fraccion de area k/n
    a = (p @ p.T) >= cos_t
    np.fill_diagonal(a, False)
    return a.astype(np.float64)


@pytest.mark.slow
def test_sphere_more_curved_than_torus() -> None:
    cfg = CurvatureConfig(max_edges=300)
    sph, tor = [], []
    for seed in range(5):
        rng = np.random.Generator(np.random.PCG64(100 + seed))
        sph.append(_kappa(_sphere_rgg(300, 20, rng), cfg, seed).mean)
        tor.append(_kappa(random_geometric_torus(300, 2, 20, rng, euclidean=False), cfg, seed).mean)
    se = np.sqrt(np.var(sph, ddof=1) / 5 + np.var(tor, ddof=1) / 5)
    assert np.mean(sph) - np.mean(tor) > 2 * se
