"""Tests analiticos del estado vacio W=0 (diseno §2: E1-E6; P§2)."""

from __future__ import annotations

import itertools
import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.seeds import seed_key
from omega.config.settings import FunctionalParams, SeedConfig
from omega.dynamics.functional import action
from omega.dynamics.gradient import grad_action
from omega.network.weights import from_upper_triangle
from omega.phases.scan import default_config, simulate_from, with_params
from omega.types import FloatArray, PhaseLabel, RunStatus


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


def _uniform(n: int, rng: np.random.Generator) -> FloatArray:
    return from_upper_triangle(rng.random(n * (n - 1) // 2), n)


def test_e1_gradient_at_zero_is_exactly_zero() -> None:
    for n, ah, gh, eta in itertools.product((5, 12, 40), (0.0, 1.0, 3.0), (0.0, 1.0, 3.0), (0.0, 0.2)):
        p = reduced_to_raw(ah, gh, n)
        p = FunctionalParams(alpha=p.alpha, beta=p.beta, gamma=p.gamma, eta=eta, mu=0.0)
        g = grad_action(np.zeros((n, n)), p)
        assert np.all(g == 0.0), (n, ah, gh, eta)


def test_e2_action_at_zero_is_exactly_zero() -> None:
    for n, ah, gh in itertools.product((5, 12, 40), (0.0, 1.0, 3.0), (0.0, 1.0, 3.0)):
        assert action(np.zeros((n, n)), reduced_to_raw(ah, gh, n)) == 0.0


def test_e3_hessian_at_zero() -> None:
    n, h = 8, 1e-6
    m = n * (n - 1) // 2
    iu = np.triu_indices(n, k=1)
    # Matriz de incidencia B (N x M): k = B x.
    b = np.zeros((n, m))
    b[iu[0], np.arange(m)] = 1.0
    b[iu[1], np.arange(m)] = 1.0
    proj = np.eye(n) - np.ones((n, n)) / n
    for ah, gh in ((0.0, 0.0), (1.5, 1.0), (3.0, 2.0)):
        p = reduced_to_raw(ah, gh, n)
        hess = np.zeros((m, m))
        for e in range(m):
            x = np.zeros(m)
            x[e] = h
            # En W=0 la rama -h no es un estado valido (W>=0): se usa diferencia hacia adelante;
            # G(0)=0 y (W^2)_{ab}=0 fuera de la diagonal para una sola arista, asi que es exacta
            # (T no contribuye en 0).
            hess[:, e] = grad_action(from_upper_triangle(x, n), p)[iu] / h
        expected = 2.0 * p.beta * np.eye(m) + 2.0 * p.gamma * b.T @ proj @ b
        assert np.max(np.abs(hess - expected)) <= 1e-6 * np.max(np.abs(expected))
        lam_min = float(np.linalg.eigvalsh((hess + hess.T) / 2.0)[0])
        assert lam_min == pytest.approx(2.0 * p.beta, rel=1e-6)


def test_e4_action_nonnegative_without_triangles() -> None:
    n = 12
    rng = _rng(4)
    for gamma in (0.0, 1.0):
        p = FunctionalParams(alpha=0.0, beta=1.0, gamma=gamma)
        for i in range(200):
            w = _uniform(n, rng)
            if i % 2 == 1:  # versiones dispersas
                mask = from_upper_triangle((rng.random(n * (n - 1) // 2) < 0.2).astype(float), n)
                w = w * mask
            s = action(w, p)
            assert s >= 0.0
            if np.any(w != 0.0):
                assert s > 0.0


def test_e5_zero_is_locally_stable() -> None:
    n = 40
    cfg0 = default_config(n, master_entropy=20261004, experiment_id=950, replicates=3)
    for ah, gh in itertools.product((0.5, 1.0, 1.5, 3.0), (0.0, 1.0)):
        p = reduced_to_raw(ah, gh, n)
        cfg = with_params(cfg0, p)
        w0 = 0.05 * _uniform(n, _rng(5))
        res = simulate_from(cfg, p, seed_key(SeedConfig(20261004, 950, 3), 0, 0), w0)
        assert res.trajectory.status is RunStatus.CONVERGED
        assert float(res.trajectory.w_final.max()) <= 1e-6
        assert res.assessment.label is PhaseLabel.A


def test_e6_above_unstable_point_flows_to_zero() -> None:
    n = 30
    cfg0 = default_config(n, master_entropy=20261004, experiment_id=950, replicates=3)
    ones = np.ones((n, n)) - np.eye(n)
    for ah in (1.2, 2.0, 3.0):
        p = reduced_to_raw(ah, 0.0, n)
        res = simulate_from(
            with_params(cfg0, p), p, seed_key(SeedConfig(20261004, 950, 3), 0, 0), (0.9 / ah) * ones
        )
        assert res.trajectory.status is RunStatus.CONVERGED
        assert float(res.trajectory.w_final.max()) <= 1e-6
