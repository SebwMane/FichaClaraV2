"""Pruebas del funcional S0 (ANALYSIS §6.2)."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from omega.config.settings import FunctionalParams
from omega.dynamics.functional import (
    action,
    degree_irregularity,
    density,
    mass,
    smoothness,
    triangles,
)
from omega.network.initialization import random_uniform_weights
from omega.network.weights import permute
from omega.types import FloatArray


def _w(n: int = 9, seed: int = 1) -> FloatArray:
    return random_uniform_weights(n, np.random.Generator(np.random.PCG64(seed)))


def test_triangles_triple_loop() -> None:
    w = _w()
    n = w.shape[0]
    loop = sum(w[i, j] * w[j, k] * w[k, i] for i, j, k in itertools.combinations(range(n), 3))
    assert triangles(w) == pytest.approx(loop, abs=1e-12)
    assert triangles(w) == pytest.approx(np.trace(w @ w @ w) / 6.0, abs=1e-12)


def test_density_degree_by_loops() -> None:
    w = _w()
    n = w.shape[0]
    d = sum(w[i, j] ** 2 for i, j in itertools.combinations(range(n), 2))
    assert density(w) == pytest.approx(d, abs=1e-12)
    k = [sum(w[i, j] for j in range(n) if j != i) for i in range(n)]
    kbar = sum(k) / n
    assert degree_irregularity(w) == pytest.approx(sum((x - kbar) ** 2 for x in k), abs=1e-12)
    assert mass(w) == pytest.approx(sum(w[i, j] for i, j in itertools.combinations(range(n), 2)))


def test_smoothness_by_loops() -> None:
    w = _w(7)
    n = w.shape[0]
    w2 = w @ w
    total = 0.0
    for i in range(n):
        for j in range(n):
            c = w2[i, i] + w2[j, j] - 2.0 * w2[i, j]
            total += w[i, j] * c
    assert smoothness(w) == pytest.approx(total, abs=1e-12)


@pytest.mark.parametrize("n", [3, 6, 10])
def test_complete_graph(n: int) -> None:
    w = np.ones((n, n)) - np.eye(n)
    assert triangles(w) == pytest.approx(n * (n - 1) * (n - 2) / 6)
    assert degree_irregularity(w) == pytest.approx(0.0, abs=1e-12)
    assert density(w) == pytest.approx(n * (n - 1) / 2)


def test_action_composition() -> None:
    w = _w()
    p = FunctionalParams(alpha=0.7, beta=1.3, gamma=0.3, eta=0.2, mu=0.4)
    expected = (
        -0.7 * triangles(w) + 1.3 * density(w) + 0.3 * degree_irregularity(w)
        + 0.2 * smoothness(w) - 0.4 * mass(w)
    )
    assert action(w, p) == pytest.approx(expected, abs=1e-12)


def test_invariance_and_no_mutation() -> None:
    w = _w()
    p = FunctionalParams(alpha=0.7, beta=1.0, gamma=0.3, eta=0.2, mu=0.1)
    copy = w.copy()
    perm = np.random.Generator(np.random.PCG64(5)).permutation(w.shape[0])
    assert action(permute(w, perm), p) == pytest.approx(action(w, p), abs=1e-10)
    assert np.array_equal(w, copy)


def test_validation() -> None:
    bad = _w()
    bad[0, 1] = 2.0
    with pytest.raises(ValueError):
        triangles(bad)
    with pytest.raises(TypeError):
        action(_w(), 1.0)  # type: ignore[arg-type]
