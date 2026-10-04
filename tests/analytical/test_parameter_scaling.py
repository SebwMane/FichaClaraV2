"""Tests analiticos del escalado de parametros (diseno §2: S1-S5; P§2, D-15)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.convert import raw_to_reduced, reduced_to_raw
from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.dynamics.evolution import evolve
from omega.dynamics.functional import action
from omega.dynamics.gradient import grad_action, lipschitz_bound
from omega.network.initialization import random_uniform_weights
from omega.types import FloatArray


def _jmi(n: int) -> FloatArray:
    return np.ones((n, n)) - np.eye(n)


def test_s1_lipschitz_independent_of_n() -> None:
    for ah, gh in ((0.0, 0.0), (1.5, 1.0), (-0.5, 3.0), (3.0, 100.0)):
        for n in (10, 64, 200, 800):
            p = reduced_to_raw(ah, gh, n)
            expected = 2.0 * p.beta * (1.0 + gh + 2.0 * abs(ah))
            assert lipschitz_bound(n, p) == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_s2_homogeneous_trajectories_invariant_in_n() -> None:
    for gh in (0.0, 3.0):
        for u0 in (0.35, 0.45):
            runs = []
            for n in (20, 50, 120):
                p = reduced_to_raw(2.5, gh, n)
                runs.append(evolve(u0 * _jmi(n), p, DynamicsConfig()))
            ref = runs[0]
            for tr in runs[1:]:
                assert tr.status is ref.status and tr.steps == ref.steps
                assert np.array_equal(tr.snapshot_steps, ref.snapshot_steps)
                assert float(np.max(np.abs(tr.scalars["mean_weight"] - ref.scalars["mean_weight"]))) <= 1e-12
                assert float(np.max(np.abs(tr.snapshots[:, :1] - ref.snapshots[:, :1]))) <= 1e-12
                assert abs(float(tr.w_final[0, 1]) - float(ref.w_final[0, 1])) <= 1e-12


@pytest.mark.parametrize("lam", [0.5, 2.0, 10.0])
def test_s3_scaling_all_coefficients(lam: float) -> None:
    n = 30
    p = reduced_to_raw(1.3, 0.7, n)
    ps = FunctionalParams(alpha=lam * p.alpha, beta=lam * p.beta, gamma=lam * p.gamma)
    w0 = random_uniform_weights(n, np.random.Generator(np.random.PCG64(np.random.SeedSequence(33))))
    a0, a1 = action(w0, p), action(w0, ps)
    assert abs(a1 - lam * a0) <= 1e-12 * max(1.0, abs(a1))
    assert float(np.max(np.abs(grad_action(w0, ps) - lam * grad_action(w0, p)))) <= 1e-12 * max(1.0, lam)
    cfg = DynamicsConfig(max_steps=300)
    t0, t1 = evolve(w0, p, cfg), evolve(w0, ps, cfg)
    assert t0.status is t1.status and t0.steps == t1.steps
    assert float(np.max(np.abs(t0.w_final - t1.w_final))) <= 1e-12
    assert float(np.max(np.abs(t0.snapshots - t1.snapshots))) <= 1e-12
    assert t1.dt == pytest.approx(t0.dt / lam, rel=1e-12)


def test_s4_critical_alpha() -> None:
    for beta in (1.0, 2.5):
        for ah_c in (1.0, 1.5, 2.0):
            for n in (64, 100, 200, 800):
                p = reduced_to_raw(ah_c, 0.0, n, beta)
                assert p.alpha == pytest.approx(2.0 * beta * ah_c / (n - 2), rel=1e-12)


def test_s5_raw_grid_is_always_above_threshold() -> None:
    """Malla cruda de M§25 (alpha/beta >= 0.1): alpha_hat = alpha (N-2)/(2 beta) > 2 desde N >= 42."""
    for n in (100, 200):
        for ratio in np.arange(1, 31) * 0.1:
            ah, _ = raw_to_reduced(FunctionalParams(alpha=float(ratio), beta=1.0), n)
            assert ah > 2.0
            assert ah >= 0.05 * (n - 2) - 1e-12
