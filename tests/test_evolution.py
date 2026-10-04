"""Pruebas de la evolucion (ANALYSIS §6.2, D-3, D-4, D-13)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.dynamics.evolution import (
    evolve,
    is_converged,
    logit,
    projected_step,
    sigmoid,
    sigmoid_step,
    snapshot_steps,
    step_scalars,
    step_size,
)
from omega.dynamics.functional import action
from omega.dynamics.gradient import grad_action, lipschitz_bound
from omega.network.initialization import random_uniform_weights
from omega.types import FloatArray, RunStatus

N = 12


def _w0(n: int = N, seed: int = 11) -> FloatArray:
    return random_uniform_weights(n, np.random.Generator(np.random.PCG64(seed)))


def test_step_size() -> None:
    p = FunctionalParams(alpha=0.3, beta=1.0, gamma=0.05)
    cfg = DynamicsConfig(dt_safety=0.4)
    assert step_size(N, p, cfg) == pytest.approx(0.4 / lipschitz_bound(N, p))
    assert step_size(N, p, DynamicsConfig(dt_mode="fixed", dt_fixed=0.01)) == 0.01
    with pytest.raises(ValueError):
        step_size(N, FunctionalParams(alpha=1.0, eta=0.1), cfg)
    assert step_size(N, FunctionalParams(alpha=1.0, eta=0.1), DynamicsConfig(dt_mode="fixed", dt_fixed=1e-3)) == 1e-3


def test_helpers() -> None:
    assert list(snapshot_steps(10, "log2")) == [0, 1, 2, 4, 8]
    assert snapshot_steps(10, "none").size == 0
    assert is_converged([1.0, 1e-12, 1e-12], 1e-10, 2)
    assert not is_converged([1e-12], 1e-10, 2)
    assert not is_converged([1e-12, 1e-3], 1e-10, 2)
    x = np.array([0.0, 0.3, 1.0])
    assert np.allclose(sigmoid(logit(x, 1e-6)), [1e-6, 0.3, 1 - 1e-6])
    w = _w0()
    g = grad_action(w, FunctionalParams(alpha=1.0))
    c = w.copy()
    out = projected_step(w, g, 0.1)
    assert np.array_equal(w, c) and out.min() >= 0 and out.max() <= 1 and np.all(np.diag(out) == 0)


def test_action_nonincreasing_auto_dt() -> None:
    p = FunctionalParams(alpha=0.2, beta=1.0, gamma=0.1, mu=0.05)
    traj = evolve(_w0(), p, DynamicsConfig(max_steps=200, tol_step=1e-14))
    s = traj.scalars["action"]
    assert s.shape == (traj.steps + 1,)
    assert np.all(np.diff(s) <= 1e-12)


def test_alpha_zero_goes_to_zero_clip() -> None:
    traj = evolve(_w0(), FunctionalParams(alpha=0.0, beta=1.0), DynamicsConfig(max_steps=20000))
    assert traj.status is RunStatus.CONVERGED
    assert np.max(traj.w_final) < 1e-9


def test_alpha_zero_sigmoid_decays() -> None:
    w0 = _w0()
    traj = evolve(w0, FunctionalParams(alpha=0.0, beta=1.0, gamma=0.1),
                  DynamicsConfig(integrator="sigmoid", max_steps=500))
    assert np.max(traj.w_final) < np.max(w0)
    assert traj.scalars["mean_weight"][-1] < traj.scalars["mean_weight"][0]


def test_alpha_hat_three_goes_to_complete() -> None:
    p = reduced_to_raw(3.0, 0.0, N)
    traj = evolve(_w0(), p, DynamicsConfig(max_steps=20000))
    assert traj.status is RunStatus.CONVERGED
    assert np.array_equal(traj.w_final, np.ones((N, N)) - np.eye(N))
    assert traj.scalars["frac_one"][-1] == 1.0


def test_uniform_point_is_stationary() -> None:
    alpha_hat = 2.5
    p = reduced_to_raw(alpha_hat, 0.0, N)
    w = (1.0 / alpha_hat) * (np.ones((N, N)) - np.eye(N))
    assert np.max(np.abs(grad_action(w, p))) < 1e-12
    traj = evolve(w, p, DynamicsConfig(max_steps=50))
    assert np.max(np.abs(traj.w_final - w)) < 1e-12
    assert traj.status is RunStatus.CONVERGED


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_nonfinite() -> None:
    w0 = _w0()
    for integ in ("clip", "sigmoid"):
        p = FunctionalParams(alpha=1e308, beta=1.0)
        cfg = DynamicsConfig(integrator=integ, dt_mode="fixed", dt_fixed=1e10)  # type: ignore[arg-type]
        traj = evolve(w0, p, cfg)
        assert traj.status is RunStatus.NONFINITE
        assert np.all(np.isfinite(traj.w_final))
        assert traj.scalars["action"].shape == (traj.steps + 1,)
    with pytest.raises(ValueError):
        bad = w0.copy()
        bad[0, 1] = bad[1, 0] = np.nan
        evolve(bad, p, DynamicsConfig())


def test_record_lengths_and_snapshots() -> None:
    p = reduced_to_raw(3.0, 1.0, N)
    w0 = _w0()
    c = w0.copy()
    traj = evolve(w0, p, DynamicsConfig(max_steps=5, tol_step=1e-300))
    assert traj.status is RunStatus.MAX_STEPS and traj.steps == 5
    for v in traj.scalars.values():
        assert v.shape == (6,)
    assert list(traj.snapshot_steps) == [0, 1, 2, 4, 5]
    assert traj.snapshots.shape == (5, N * (N - 1) // 2)
    assert traj.tau == pytest.approx(5 * traj.dt)
    assert np.array_equal(w0, c)
    assert np.allclose(traj.snapshots[-1], traj.w_final[np.triu_indices(N, 1)])
    none = evolve(w0, p, DynamicsConfig(max_steps=5, snapshot_schedule="none", tol_step=1e-300))
    assert list(none.snapshot_steps) == [5]


def test_step_scalars_and_eta_requires_fixed() -> None:
    p = FunctionalParams(alpha=0.3, gamma=0.1)
    w = _w0()
    sc = step_scalars(w, w, p)
    assert sc["max_dw"] == 0.0 and sc["action"] == pytest.approx(action(w, p))
    pe = FunctionalParams(alpha=0.3, eta=0.01)
    with pytest.raises(ValueError):
        evolve(w, pe, DynamicsConfig())
    traj = evolve(w, pe, DynamicsConfig(dt_mode="fixed", dt_fixed=1e-3, max_steps=20))
    assert np.all(np.isfinite(traj.w_final))


def test_sigmoid_step_decreases_action() -> None:
    p = FunctionalParams(alpha=0.2, beta=1.0, gamma=0.1)
    w = _w0()
    theta = logit(w, 1e-6)
    dt = 0.5 / lipschitz_bound(N, p)
    th2 = sigmoid_step(theta, p, dt)
    w1 = sigmoid(theta)
    w2 = sigmoid(th2)
    np.fill_diagonal(w1, 0.0)
    np.fill_diagonal(w2, 0.0)
    assert action(w2, p) <= action(w1, p) + 1e-12
