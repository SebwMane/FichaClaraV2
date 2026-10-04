"""Tests analiticos del estado completo W=J-I (diseno §2: C1-C5; P§2)."""

from __future__ import annotations

from math import comb

import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.seeds import seed_key
from omega.config.settings import SeedConfig
from omega.dynamics.evolution import projected_step
from omega.dynamics.functional import action, degree_irregularity, density, triangles
from omega.dynamics.gradient import grad_action
from omega.network.weights import from_upper_triangle
from omega.phases.scan import default_config, simulate_from, with_params
from omega.types import FloatArray, PhaseLabel, RunStatus


def _jmi(n: int) -> FloatArray:
    return np.ones((n, n)) - np.eye(n)


def _tol(a: float, b: float) -> float:
    return 1e-12 * max(1.0, abs(a), abs(b))


def test_c1_counts() -> None:
    for n in (5, 40, 200):
        w = _jmi(n)
        assert abs(triangles(w) - comb(n, 3)) <= _tol(triangles(w), comb(n, 3))
        assert abs(density(w) - comb(n, 2)) <= _tol(density(w), comb(n, 2))
        assert abs(degree_irregularity(w)) <= 1e-12


def test_c2_gradient_independent_of_gamma() -> None:
    n = 30
    w = _jmi(n)
    for ah in (0.0, 0.7, 1.0, 2.5):
        for gh in (0.0, 1.0, 100.0):
            p = reduced_to_raw(ah, gh, n)
            g = grad_action(w, p)
            expected = 2.0 * p.beta * (1.0 - ah) * w
            assert np.max(np.abs(g - expected)) <= 1e-12 * max(1.0, float(np.max(np.abs(expected))))


def test_c3_kkt_threshold_at_alpha_hat_one() -> None:
    n, dt = 40, 0.01
    w = _jmi(n)
    for ah in (1.0, 1.1, 3.0):
        for gh in (0.0, 1.0, 100.0):
            p = reduced_to_raw(ah, gh, n)
            out = projected_step(w, grad_action(w, p), dt)
            assert np.array_equal(out, w), (ah, gh)
    p = reduced_to_raw(0.9, 1.0, n)
    out = projected_step(w, grad_action(w, p), dt)
    off = ~np.eye(n, dtype=bool)
    assert np.max(np.abs((w - out)[off] - dt * 2.0 * p.beta * 0.1)) <= 1e-12


def _run(n: int, ah: float, w0: FloatArray) -> tuple[RunStatus, FloatArray, PhaseLabel]:
    p = reduced_to_raw(ah, 0.0, n)
    cfg = with_params(default_config(n, master_entropy=20261004, experiment_id=950, replicates=3), p)
    res = simulate_from(cfg, p, seed_key(SeedConfig(20261004, 950, 3), 0, 0), w0)
    return res.trajectory.status, res.trajectory.w_final, res.assessment.label


def test_c4_local_stability() -> None:
    n = 40
    st, wf, lab = _run(n, 1.05, 0.99 * _jmi(n))
    assert st is RunStatus.CONVERGED and np.array_equal(wf, _jmi(n)) and lab is PhaseLabel.E
    st, wf, _ = _run(n, 0.95, 0.99 * _jmi(n))
    assert st is RunStatus.CONVERGED and float(wf.max()) <= 1e-6
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(44)))
    pert = from_upper_triangle(0.01 * rng.random(n * (n - 1) // 2), n)
    st, wf, _ = _run(n, 1.2, _jmi(n) - pert)
    assert st is RunStatus.CONVERGED and np.array_equal(wf, _jmi(n))


@pytest.mark.parametrize("ah,sign", [(1.4, 1), (1.5, 0), (1.6, -1)])
def test_c5_threshold_three_halves(ah: float, sign: int) -> None:
    for n in (10, 40, 200):
        p = reduced_to_raw(ah, 0.0, n)
        diff = action(_jmi(n), p) - action(np.zeros((n, n)), p)
        expected = p.beta * comb(n, 2) * (1.0 - 2.0 * ah / 3.0)
        assert abs(diff - expected) <= 1e-9 * comb(n, 2)
        if sign == 0:
            assert abs(diff) <= 1e-9 * comb(n, 2)
        else:
            assert np.sign(diff) == sign
