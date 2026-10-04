"""Tests analiticos de la variedad homogenea W=w(J-I) (diseno §2: U1-U7; P§2)."""

from __future__ import annotations

from math import comb

import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.seeds import seed_key
from omega.config.settings import DynamicsConfig, FunctionalParams, SeedConfig
from omega.dynamics.evolution import evolve, step_size
from omega.dynamics.functional import action
from omega.dynamics.gradient import grad_action, grad_degree_irregularity, lipschitz_bound
from omega.phases.scan import default_config, simulate, simulate_from, with_params
from omega.types import FloatArray, PhaseLabel, RunStatus

SEEDS = SeedConfig(20261004, 950, 3)


def _jmi(n: int) -> FloatArray:
    return np.ones((n, n)) - np.eye(n)


def _hval(w: float, n: int, p: FunctionalParams) -> float:
    """H(w) = S0(w (J-I))."""
    return action(w * _jmi(n), p)


def test_u1_homogeneous_action() -> None:
    n = 12
    for gh in (0.0, 1.0, 100.0):
        for w in (0.1, 0.37, 0.8):
            p = reduced_to_raw(1.3, gh, n)
            expected = -p.alpha * comb(n, 3) * w**3 + p.beta * comb(n, 2) * w**2
            got = action(w * _jmi(n), p)
            assert abs(got - expected) <= 1e-12 * max(1.0, abs(got), abs(expected))


def test_u2_homogeneous_gradient() -> None:
    n = 12
    for gh in (0.0, 1.0, 100.0):
        for w in (0.1, 0.37, 0.8):
            ah = 1.7
            p = reduced_to_raw(ah, gh, n)
            x = w * _jmi(n)
            expected = 2.0 * p.beta * w * (1.0 - ah * w) * _jmi(n)
            assert np.max(np.abs(grad_action(x, p) - expected)) <= 1e-12
            assert float(np.max(np.abs(grad_degree_irregularity(x)))) <= 1e-12


def test_u3_stationary_point() -> None:
    n = 30
    for ah in (1.25, 2.0, 3.0):
        p = reduced_to_raw(ah, 1.0, n)
        g = grad_action((1.0 / ah) * _jmi(n), p)
        assert float(np.max(np.abs(g))) <= 1e-12 * 2.0 * p.beta
    # Para alpha_hat <= 1 el gradiente homogeneo es estrictamente positivo en (0,1).
    for ah in (0.5, 1.0):
        p = reduced_to_raw(ah, 0.0, n)
        for w in np.linspace(0.01, 0.99, 25):
            g = grad_action(float(w) * _jmi(n), p)
            assert float(g[0, 1]) > 0.0


def test_u4_stationary_point_is_a_maximum() -> None:
    n, h = 20, 0.01
    for ah in (1.25, 2.0, 3.0):
        p = reduced_to_raw(ah, 1.0, n)
        ws = 1.0 / ah

        assert _hval(ws + h, n, p) < _hval(ws, n, p)
        assert _hval(ws - h, n, p) < _hval(ws, n, p)
        h2 = (_hval(ws + h, n, p) - 2.0 * _hval(ws, n, p) + _hval(ws - h, n, p)) / h**2
        assert h2 == pytest.approx(-2.0 * p.beta * comb(n, 2), rel=1e-6)


def test_u5_three_distinct_thresholds() -> None:
    got = [round(reduced_to_raw(a, 0.0, 200).alpha, 4) for a in (1.0, 1.5, 2.0)]
    assert got == [0.0101, 0.0152, 0.0202]
    assert round(reduced_to_raw(1.0, 0.0, 100).alpha, 4) == 0.0204
    # alpha_hat = 3/2: H(1) = H(0).
    n = 200
    p = reduced_to_raw(1.5, 0.0, n)
    assert abs(action(_jmi(n), p) - action(np.zeros((n, n)), p)) <= 1e-9 * comb(n, 2)
    # alpha_hat = 2: w0 = 1/2 estacionario; 1.98 -> 0; 2.02 -> J-I.
    w0 = 0.5 * _jmi(n)
    assert float(np.max(np.abs(grad_action(w0, reduced_to_raw(2.0, 0.0, n))))) <= 1e-12
    for ah, target in ((1.98, 0.0), (2.02, 1.0)):
        p = reduced_to_raw(ah, 0.0, n)
        cfg = with_params(default_config(n, master_entropy=20261004, experiment_id=950, replicates=3), p)
        res = simulate_from(cfg, p, seed_key(SEEDS, 0, 0), w0)
        assert res.trajectory.status is RunStatus.CONVERGED
        assert float(np.max(np.abs(res.trajectory.w_final - target * _jmi(n)))) <= 1e-6


def _labels_from_uniform(n: int, ah: float) -> list[PhaseLabel]:
    p = reduced_to_raw(ah, 0.0, n)
    cfg = with_params(default_config(n, master_entropy=20261004, experiment_id=950, replicates=3), p)
    return [simulate(cfg, p, seed_key(SEEDS, 0, r)).assessment.label for r in range(3)]


def _check_u5_random(n: int) -> None:
    for ah in (1.2, 1.6, 1.9):
        assert _labels_from_uniform(n, ah) == [PhaseLabel.A] * 3, ah
    assert _labels_from_uniform(n, 2.1) == [PhaseLabel.E] * 3
    # J-I es localmente estable (y preferido globalmente si alpha_hat > 3/2) pero inaccesible
    # desde U(0,1) para alpha_hat < 2.
    for ah in (1.6, 1.9):
        p = reduced_to_raw(ah, 0.0, n)
        cfg = with_params(default_config(n, master_entropy=20261004, experiment_id=950, replicates=3), p)
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(55)))
        pert = np.triu(0.01 * rng.random((n, n)), 1)
        w0 = _jmi(n) - (pert + pert.T)
        res = simulate_from(cfg, p, seed_key(SEEDS, 0, 0), w0)
        assert np.array_equal(res.trajectory.w_final, _jmi(n))


def test_u5_from_uniform_n100() -> None:
    _check_u5_random(100)


@pytest.mark.slow
def test_u5_from_uniform_n200() -> None:
    _check_u5_random(200)


def test_u6_gamma_does_not_change_homogeneous_trajectory() -> None:
    n = 50
    cases = [(0.5, 0.9), (1.2, 0.8), (1.2, 0.9), (2.0, 0.3), (2.0, 0.49), (2.0, 0.51)]
    for ah, u0 in cases:
        dt = 0.5 / lipschitz_bound(n, reduced_to_raw(ah, 100.0, n))
        cfg = DynamicsConfig(dt_mode="fixed", dt_fixed=dt, max_steps=20000)
        runs = [
            evolve(u0 * _jmi(n), reduced_to_raw(ah, gh, n), cfg) for gh in (0.0, 1.0, 100.0)
        ]
        ref = runs[0]
        for tr in runs[1:]:
            assert tr.status is ref.status and tr.steps == ref.steps
            assert np.array_equal(tr.snapshot_steps, ref.snapshot_steps)
            assert float(np.max(np.abs(tr.snapshots - ref.snapshots))) <= 1e-12
        # Recursion escalar u_{t+1} = clip(u_t - dt 2 beta u_t (1 - alpha_hat u_t), 0, 1).
        u = u0
        last = 0
        for step, snap in zip(ref.snapshot_steps, ref.snapshots):
            while last < int(step):
                u = min(1.0, max(0.0, u - dt * 2.0 * u * (1.0 - ah * u)))
                last += 1
            assert float(np.max(np.abs(snap - u))) <= 1e-12


def test_u7_lipschitz_and_auto_dt_depend_on_gamma() -> None:
    n = 40
    dts = []
    for ah in (0.5, 2.0):
        for gh in (0.0, 1.0, 100.0):
            p = reduced_to_raw(ah, gh, n)
            lip = lipschitz_bound(n, p)
            assert lip == pytest.approx(2.0 * p.beta * (1.0 + gh + 2.0 * ah), rel=1e-12)
            dt = step_size(n, p, DynamicsConfig())
            assert dt == pytest.approx(0.5 / lip, rel=1e-12)
            dts.append(dt)
    assert len(set(dts[:3])) == 3 and len(set(dts[3:])) == 3
