"""Pruebas del gradiente por arista (ANALYSIS §1.2, §1.5, §1.6, D-25)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.settings import FunctionalParams
from omega.dynamics.functional import action
from omega.dynamics.gradient import (
    grad_action,
    grad_degree_irregularity,
    grad_density,
    grad_mass,
    grad_smoothness,
    grad_triangles,
    lipschitz_bound,
    theta_gradient,
)
from omega.network.initialization import random_uniform_weights
from omega.network.weights import from_upper_triangle, permute, upper_triangle
from omega.types import FloatArray


def _w(n: int = 12, seed: int = 3) -> FloatArray:
    return random_uniform_weights(n, np.random.Generator(np.random.PCG64(seed)))


def _fd(w: FloatArray, p: FunctionalParams, h: float = 1e-6) -> FloatArray:
    n = w.shape[0]
    v = upper_triangle(w)
    g = np.zeros_like(v)
    for e in range(v.size):
        vp, vm = v.copy(), v.copy()
        vp[e] += h
        vm[e] -= h
        # fuera de [0,1] la validacion falla: se usa un W interior (0.1..0.9)
        g[e] = (action(from_upper_triangle(vp, n), p) - action(from_upper_triangle(vm, n), p)) / (2 * h)
    return from_upper_triangle(g, n)


def _interior(n: int = 12) -> FloatArray:
    w = _w(n)
    return from_upper_triangle(0.1 + 0.8 * upper_triangle(w), n)


@pytest.mark.parametrize(
    "p",
    [
        FunctionalParams(alpha=0.7, beta=1.0, gamma=0.3, eta=0.2, mu=0.15),
        FunctionalParams(alpha=-0.4, beta=0.5, gamma=0.0),
        FunctionalParams(alpha=1.1, beta=2.0, gamma=1.5, mu=-0.3),
    ],
)
def test_gradient_vs_finite_differences(p: FunctionalParams) -> None:
    w = _interior()
    g = grad_action(w, p)
    fd = _fd(w, p)
    err = np.max(np.abs(g - fd)) / np.max(np.abs(fd))
    assert err < 1e-6


def test_term_gradients_symmetric_zero_diag() -> None:
    w = _w()
    for g in (
        grad_triangles(w), grad_density(w), grad_degree_irregularity(w),
        grad_smoothness(w), grad_mass(w), grad_action(w, FunctionalParams(alpha=1.0, eta=0.1, mu=0.1)),
    ):
        assert np.array_equal(g, g.T)
        assert np.all(np.diag(g) == 0.0)


def test_equivariance() -> None:
    w = _w(10)
    p = FunctionalParams(alpha=0.7, beta=1.0, gamma=0.3, eta=0.2, mu=0.1)
    perm = np.random.Generator(np.random.PCG64(9)).permutation(10)
    g = grad_action(w, p)
    assert np.allclose(grad_action(permute(w, perm), p), g[np.ix_(perm, perm)], atol=1e-12)
    assert action(permute(w, perm), p) == pytest.approx(action(w, p), abs=1e-10)


def test_degree_hessian_max_eigenvalue() -> None:
    n = 8
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    b = np.zeros((len(pairs), n))
    for e, (i, j) in enumerate(pairs):
        b[e, i] = b[e, j] = 1.0
    proj = np.eye(n) - np.ones((n, n)) / n
    gamma = 0.7
    h = 2.0 * gamma * b @ proj @ b.T
    assert np.linalg.eigvalsh(h)[-1] == pytest.approx(2.0 * gamma * (n - 2), rel=1e-10)


def test_lipschitz() -> None:
    p = FunctionalParams(alpha=-2.0, beta=1.5, gamma=0.5)
    assert lipschitz_bound(10, p) == pytest.approx(3.0 + 2 * 2.5 * 8)
    with pytest.raises(ValueError):
        lipschitz_bound(10, FunctionalParams(alpha=1.0, eta=0.1))
    with pytest.raises(ValueError):
        lipschitz_bound(1, p)


def test_theta_gradient() -> None:
    w = _w(6)
    g = grad_action(w, FunctionalParams(alpha=1.0, gamma=0.2))
    t = theta_gradient(g, w)
    assert np.allclose(t, g * w * (1 - w))
    assert np.all(np.diag(t) == 0.0)
    with pytest.raises(ValueError):
        theta_gradient(g[:3, :3], w)


def test_no_mutation() -> None:
    w = _w(6)
    c = w.copy()
    grad_action(w, FunctionalParams(alpha=1.0, eta=0.1, mu=0.1, gamma=0.1))
    assert np.array_equal(w, c)
