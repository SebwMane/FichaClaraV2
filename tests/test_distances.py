"""Pruebas de omega.geometry.distances (ANALYSIS §6.2: distancias)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.geometry.distances import (
    ball,
    ball_counts,
    diameter,
    distance_matrix,
    edge_length_matrix,
    giant_nodes_from_distances,
    hop_distance_matrix,
    path_length,
    threshold_adjacency,
)
from omega.network.initialization import random_uniform_weights
from omega.types import FloatArray

EPS = 1e-9


def _chain(weights: list[float]) -> FloatArray:
    n = len(weights) + 1
    w = np.zeros((n, n))
    for i, x in enumerate(weights):
        w[i, i + 1] = w[i + 1, i] = x
    return w


def test_chain_analytic_distances() -> None:
    ws = [1.0, 0.5, 0.25, 1.0]
    d = distance_matrix(_chain(ws), 0.0, EPS)
    lengths = [1.0 / (x + EPS) for x in ws]
    for i in range(5):
        for j in range(5):
            expected = abs(sum(lengths[min(i, j) : max(i, j)]))
            assert d[i, j] == pytest.approx(expected, rel=1e-12, abs=1e-12)
    assert np.array_equal(d, d.T)


def test_edge_lengths_and_threshold_is_strict() -> None:
    w = _chain([0.5, 0.2])
    a = threshold_adjacency(w, 0.2)  # W>0.2 excluye la arista de peso 0.2
    assert a[0, 1] and not a[1, 2]
    g = edge_length_matrix(w, a, EPS)
    assert g[0, 1] == pytest.approx(1.0 / (0.5 + EPS))
    assert g[1, 2] == 0.0


def test_triangle_inequality_random() -> None:
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(11)))
    w = random_uniform_weights(60, rng)
    d = distance_matrix(w, 0.6, EPS)
    assert np.all(np.isfinite(d))
    viol = d[:, None, :] - (d[:, :, None] + d[None, :, :])  # d_ik - d_ij - d_jk
    assert viol.max() <= 1e-9
    assert np.all(np.diag(d) == 0.0)


def test_inf_between_components() -> None:
    w = np.zeros((6, 6))
    for i, j in ((0, 1), (1, 2), (3, 4)):
        w[i, j] = w[j, i] = 1.0
    d = distance_matrix(w, 0.5, EPS)
    assert np.isinf(d[0, 3]) and np.isinf(d[2, 4]) and np.isinf(d[5, 0])
    assert np.isfinite(d[0, 2])
    h = hop_distance_matrix(threshold_adjacency(w, 0.5))
    assert h[0, 2] == 2.0 and np.isinf(h[0, 3])
    nodes = giant_nodes_from_distances(d)
    assert nodes.tolist() == [0, 1, 2]
    assert path_length(d) == pytest.approx((2 * (1 + 1 + 2) / 6) * (1 - EPS), rel=1e-6)
    assert diameter(d) == pytest.approx(2.0, rel=1e-6)


def test_ball_counts_monotone_and_consistent() -> None:
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(12)))
    w = random_uniform_weights(40, rng)
    d = distance_matrix(w, 0.5, EPS)
    radii = np.linspace(0.5, 6.0, 25)
    c = ball_counts(d, radii)
    assert c.shape == (40, 25)
    assert np.all(np.diff(c, axis=1) >= 0)
    assert np.all(c >= 1)
    for k in (0, 10, 24):
        assert c[7, k] == ball(d, 7, float(radii[k])).size
    assert 7 in ball(d, 7, 0.0).tolist()


def test_inputs_not_mutated_and_validation() -> None:
    w = _chain([1.0, 1.0])
    w0 = w.copy()
    distance_matrix(w, 0.5, EPS)
    assert np.array_equal(w, w0)
    with pytest.raises(ValueError):
        distance_matrix(w, 1.5, EPS)
    with pytest.raises(ValueError):
        distance_matrix(w, 0.5, 0.0)
    with pytest.raises(ValueError):
        ball_counts(np.zeros((3, 3)), np.array([]))
