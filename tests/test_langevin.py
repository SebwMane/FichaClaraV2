"""WP-D: Langevin (Theta es temperatura estadistica, no tiempo fisico)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, LangevinConfig
from omega.contracts import ChainResult
from omega.dynamics.gradient import lipschitz_bound
from omega.network.weights import from_upper_triangle, validate_weight_matrix
from omega.statistics.ensemble import effective_sample_size
from omega.statistics.langevin import (
    chain_observables,
    langevin_step,
    reflect_unit,
    run_langevin,
    temperature,

)


def _rng(seed: int) -> np.random.Generator:
    """Generador PCG64 explícito (prohibido np.random global)."""
    return np.random.Generator(np.random.PCG64(seed))


def _w0(n: int, seed: int) -> np.ndarray:
    v = _rng(seed).uniform(size=n * (n - 1) // 2)
    return from_upper_triangle(v, n)


def test_reflect_unit_values_and_purity() -> None:
    x = np.array([-0.3, 0.0, 0.4, 1.0, 1.2, 2.5, 3.1, -2.2])
    x0 = x.copy()
    y = reflect_unit(x)
    np.testing.assert_allclose(y, [0.3, 0.0, 0.4, 1.0, 0.8, 0.5, 0.9, 0.2], atol=1e-12)
    np.testing.assert_array_equal(x, x0)
    z = _rng(0).uniform(-50, 50, 1000)
    r = reflect_unit(z)
    assert r.min() >= 0.0 and r.max() <= 1.0


def test_temperature_scales_with_beta() -> None:
    assert temperature(0.3, 2.0) == pytest.approx(0.6)
    with pytest.raises(ValueError):
        temperature(0.0, 1.0)


def test_eta_nonzero_rejected() -> None:
    p = FunctionalParams(alpha=0.1, beta=1.0, eta=0.1)
    w0 = _w0(5, 0)
    with pytest.raises(ValueError):
        run_langevin(w0, p, LangevinConfig(theta_hat=0.1, n_steps=10), 0.1, _rng(0))
    with pytest.raises(ValueError):
        langevin_step(w0[np.triu_indices(5, 1)], 5, p, 0.01, 0.1, _rng(0), "reflect")


def test_step_pure_bounded_and_deterministic() -> None:
    p = FunctionalParams(alpha=0.5, beta=1.0, gamma=0.2)
    u = _w0(6, 1)[np.triu_indices(6, 1)]
    u0 = u.copy()
    a = langevin_step(u, 6, p, 0.01, 5.0, _rng(7), "reflect")
    b = langevin_step(u, 6, p, 0.01, 5.0, _rng(7), "reflect")
    np.testing.assert_array_equal(u, u0)
    np.testing.assert_array_equal(a, b)
    assert a.min() >= 0.0 and a.max() <= 1.0 and a.shape == u.shape


def test_run_langevin_contract_and_invariants() -> None:
    p = FunctionalParams(alpha=0.5, beta=1.0, gamma=0.1)
    w0 = _w0(6, 2)
    w0c = w0.copy()
    cfg = LangevinConfig(theta_hat=0.3, n_steps=400, thin=10)
    r = run_langevin(w0, p, cfg, 0.1, _rng(3), n_states=3)
    np.testing.assert_array_equal(w0, w0c)
    assert isinstance(r, ChainResult) and r.engine is Engine.LANGEVIN
    assert r.theta == pytest.approx(0.3) and r.acceptance == 1.0 and r.n_steps == 400
    assert r.step_size == pytest.approx(0.1 / lipschitz_bound(6, p))
    assert r.sample_steps.shape == (40,) and r.sample_steps[0] == 10 and r.sample_steps[-1] == 400
    for key in ("action", "mean_weight", "binary_density", "triangle_density", "strength_cv"):
        assert r.samples[key].shape == (40,)
    validate_weight_matrix(r.w_final)
    assert len(r.states) == 3
    for s in r.states:
        validate_weight_matrix(s)
    r2 = run_langevin(w0, p, cfg, 0.1, _rng(3), n_states=3)
    np.testing.assert_array_equal(r.w_final, r2.w_final)


def test_chain_observables_known_values() -> None:
    n = 5
    w = np.ones((n, n)) - np.eye(n)
    o = chain_observables(w, FunctionalParams(alpha=0.0, beta=1.0), 0.5)
    assert o["mean_weight"] == pytest.approx(1.0)
    assert o["binary_density"] == pytest.approx(1.0)
    assert o["triangle_density"] == pytest.approx(1.0)
    assert o["strength_cv"] == pytest.approx(0.0)
    assert o["action"] == pytest.approx(10.0)


def test_gaussian_variance_far_from_borders() -> None:
    """alpha=gamma=0, H = w^2 - w/... con mu=1: minimo en 1/2, varianza Theta/(2 beta) = Theta_hat/2."""
    n, theta_hat = 6, 0.01
    p = FunctionalParams(alpha=0.0, beta=1.0, gamma=0.0, mu=1.0)
    w = np.full((n, n), 0.5) - 0.5 * np.eye(n)
    u = w[np.triu_indices(n, 1)]
    rng = _rng(11)
    dt = 0.02 / lipschitz_bound(n, p)
    theta = temperature(theta_hat, p.beta)
    xs = []
    for t in range(40_000):
        u = langevin_step(u, n, p, dt, theta, rng, "reflect")
        if t >= 2000 and t % 4 == 0:
            xs.append(u.copy())
    x = np.array(xs)
    assert abs(x.mean() - 0.5) < 0.005
    var = float(x.var())
    # Euler-Maruyama: sesgo relativo O(k*dt/2) ~ 1% (k=2 beta); SE estadistico ~1-2%
    assert var == pytest.approx(theta_hat / 2.0, rel=0.06)


def _quad_mean(theta_hat: float) -> float:
    num = quad(lambda w: w * np.exp(-w * w / theta_hat), 0, 1)[0]
    den = quad(lambda w: np.exp(-w * w / theta_hat), 0, 1)[0]
    return float(num / den)


def _lan_mean(w0: np.ndarray, p: FunctionalParams, theta_hat: float, dts: float, steps: int,
              seed: int) -> tuple[float, float]:
    cfg = LangevinConfig(theta_hat=theta_hat, n_steps=steps, dt_safety=dts, thin=5)
    x = run_langevin(w0, p, cfg, 0.1, _rng(seed)).samples["mean_weight"]
    x = x[len(x) // 10 :]
    return float(x.mean()), float(x.std() / np.sqrt(effective_sample_size(x)))


@pytest.mark.parametrize("theta_hat", [0.1, 1.0])
def test_gibbs_independent_edges_quadrature(theta_hat: float) -> None:
    """N=8 (M=28), alpha=gamma=0: <w> = cuadratura; Euler-Maruyama (orden 1 en dt) se
    extrapola a dt->0 por Richardson (2 m(h) - m(2h)) y se compara con 3 SE combinado."""
    p = FunctionalParams(alpha=0.0, beta=1.0, gamma=0.0)
    w0 = _w0(8, 5)
    m1, s1 = _lan_mean(w0, p, theta_hat, 0.02, 25_000, 21)
    m2, s2 = _lan_mean(w0, p, theta_hat, 0.04, 25_000, 22)
    extrap = 2.0 * m1 - m2
    se = float(np.sqrt(4.0 * s1**2 + s2**2))
    exact = _quad_mean(theta_hat)
    assert abs(extrap - exact) < 3.0 * se
    # el sesgo de discretizacion existe y crece con dt: m(2h) esta mas lejos de la exacta que m(h)
    assert abs(m2 - exact) > abs(m1 - exact) - 3.0 * s1


def test_project_creates_atoms_reflect_does_not() -> None:
    n = 6
    p = FunctionalParams(alpha=0.0, beta=1.0, gamma=0.0)
    dt = 0.1 / lipschitz_bound(n, p)
    theta = temperature(1.0, p.beta)
    frac = {}
    for bnd in ("reflect", "project"):
        rng = _rng(4)
        u = np.full(n * (n - 1) // 2, 0.3)
        zeros = 0
        total = 0
        for t in range(3000):
            u = langevin_step(u, n, p, dt, theta, rng, bnd)
            if t >= 500:
                zeros += int(np.sum(u == 0.0))
                total += u.size
        frac[bnd] = zeros / total
    assert frac["project"] > 0.01
    assert frac["reflect"] == 0.0
