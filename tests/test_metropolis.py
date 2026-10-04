"""WP-D: Metropolis por aristas (balance detallado, adaptacion solo en burn-in)."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, MetropolisConfig
from omega.dynamics.functional import action
from omega.network.weights import from_upper_triangle, validate_weight_matrix
from omega.statistics.ensemble import effective_sample_size
from omega.statistics.metropolis import delta_action_edge, metropolis_sweep, run_metropolis


def _w0(n: int, seed: int) -> np.ndarray:
    return from_upper_triangle(np.random.default_rng(seed).uniform(size=n * (n - 1) // 2), n)


def test_delta_action_matches_action_difference_200() -> None:
    rng = np.random.default_rng(0)
    for _ in range(200):
        n = int(rng.integers(3, 10))
        p = FunctionalParams(
            alpha=float(rng.normal()), beta=float(rng.uniform(0.2, 2)), gamma=float(rng.uniform(0, 1.5)),
            mu=float(rng.normal()),
        )
        w = _w0(n, int(rng.integers(1 << 30)))
        a, b = (int(i) for i in rng.choice(n, 2, replace=False))
        new = float(rng.uniform())
        w2 = w.copy()
        w2[a, b] = w2[b, a] = new
        exact = action(w2, p) - action(w, p)
        got = delta_action_edge(w, a, b, new, p, w.sum(axis=1))
        assert got == pytest.approx(exact, rel=1e-10, abs=1e-12)


def test_eta_rejected() -> None:
    p = FunctionalParams(alpha=0.1, beta=1.0, eta=0.2)
    w = _w0(4, 1)
    with pytest.raises(ValueError):
        delta_action_edge(w, 0, 1, 0.5, p, w.sum(axis=1))
    with pytest.raises(ValueError):
        metropolis_sweep(w, p, 0.1, 0.2, np.random.default_rng(0))
    with pytest.raises(ValueError):
        run_metropolis(w, p, MetropolisConfig(theta_hat=0.1, n_sweeps=5), 0.1, np.random.default_rng(0))


def test_sweep_pure_symmetric_bounded() -> None:
    p = FunctionalParams(alpha=0.5, beta=1.0, gamma=0.3)
    w = _w0(7, 2)
    w0 = w.copy()
    out, acc = metropolis_sweep(w, p, 0.3, 0.5, np.random.default_rng(1))
    np.testing.assert_array_equal(w, w0)
    validate_weight_matrix(out)
    assert 0 <= acc <= 21
    out2, acc2 = metropolis_sweep(w, p, 0.3, 0.5, np.random.default_rng(1))
    np.testing.assert_array_equal(out, out2)
    assert acc == acc2


def _q(x: float, y: float, s: float) -> float:
    """Densidad de la propuesta B(x + s U(-1,1)) en y (imagenes por reflexion)."""
    tot = 0.0
    for m in range(-3, 4):
        for img in (y + 2 * m, -y + 2 * m):
            if abs(img - x) < s:
                tot += 1.0 / (2.0 * s)
    return tot


def test_detailed_balance_single_edge() -> None:
    """pi(x) q(x,y) a(x->y) == pi(y) q(y,x) a(y->x) con pi ~ exp(-H/Theta), propuesta reflejada simetrica."""
    p = FunctionalParams(alpha=0.0, beta=1.0, gamma=0.0, mu=0.7)
    theta = 0.4
    rng = np.random.default_rng(3)
    for _ in range(100):
        x, y = rng.uniform(size=2)
        s = float(rng.uniform(0.1, 1.0))
        wx = np.array([[0.0, x], [x, 0.0]])
        wy = np.array([[0.0, y], [y, 0.0]])
        dxy = delta_action_edge(wx, 0, 1, float(y), p, wx.sum(axis=1))
        dyx = delta_action_edge(wy, 0, 1, float(x), p, wy.sum(axis=1))
        assert dxy == pytest.approx(-dyx, abs=1e-12)
        pi_x, pi_y = np.exp(-action(wx, p) / theta), np.exp(-action(wy, p) / theta)
        lhs = pi_x * _q(x, y, s) * min(1.0, np.exp(-dxy / theta))
        rhs = pi_y * _q(y, x, s) * min(1.0, np.exp(-dyx / theta))
        assert lhs == pytest.approx(rhs, rel=1e-10, abs=1e-14)


def test_reflected_proposal_empirically_symmetric() -> None:
    from omega.statistics.langevin import reflect_unit

    rng = np.random.default_rng(5)
    s, x, y, h = 0.5, 0.15, 0.4, 0.02
    n = 400_000
    px = np.mean(np.abs(reflect_unit(x + s * rng.uniform(-1, 1, n)) - y) < h)
    py = np.mean(np.abs(reflect_unit(y + s * rng.uniform(-1, 1, n)) - x) < h)
    assert px == pytest.approx(py, rel=0.05)
    assert px == pytest.approx(2 * h * _q(x, y, s), rel=0.05)


def test_step_fixed_after_burn_in_and_bounds() -> None:
    p = FunctionalParams(alpha=0.0, beta=1.0)
    w0 = _w0(6, 1)
    r0 = run_metropolis(w0, p, MetropolisConfig(theta_hat=0.1, n_sweeps=50, step_init=0.3,
                                                burn_in_fraction=0.0), 0.1, np.random.default_rng(2))
    assert r0.step_size == pytest.approx(0.3)
    r1 = run_metropolis(w0, p, MetropolisConfig(theta_hat=0.01, n_sweeps=400, step_init=1.0,
                                                burn_in_fraction=0.5), 0.1, np.random.default_rng(2))
    assert 1e-4 <= r1.step_size <= 1.0 and r1.step_size < 1.0
    assert r1.engine is Engine.METROPOLIS and r1.theta == pytest.approx(0.01)
    assert 0.0 < r1.acceptance <= 1.0
    # adaptacion: paso inicial demasiado grande a Theta baja -> aceptacion baja -> reduce
    validate_weight_matrix(r1.w_final)


def test_run_contract_states_determinism() -> None:
    p = FunctionalParams(alpha=0.5, beta=1.0, gamma=0.1)
    w0 = _w0(6, 3)
    cfg = MetropolisConfig(theta_hat=0.3, n_sweeps=100, thin=2)
    a = run_metropolis(w0, p, cfg, 0.1, np.random.default_rng(9), n_states=4)
    b = run_metropolis(w0, p, cfg, 0.1, np.random.default_rng(9), n_states=4)
    assert len(a.states) == 4 and a.sample_steps.shape == (50,)
    np.testing.assert_array_equal(a.w_final, b.w_final)
    np.testing.assert_array_equal(a.samples["action"], b.samples["action"])


@pytest.mark.parametrize("theta_hat", [0.1, 1.0])
def test_gibbs_independent_edges_quadrature(theta_hat: float) -> None:
    p = FunctionalParams(alpha=0.0, beta=1.0, gamma=0.0)
    cfg = MetropolisConfig(theta_hat=theta_hat, n_sweeps=4000)
    r = run_metropolis(_w0(8, 5), p, cfg, 0.1, np.random.default_rng(3))
    x = r.samples["mean_weight"][2000:]
    se = float(x.std() / np.sqrt(effective_sample_size(x)))
    num = quad(lambda w: w * np.exp(-w * w / theta_hat), 0, 1)[0]
    den = quad(lambda w: np.exp(-w * w / theta_hat), 0, 1)[0]
    assert abs(x.mean() - num / den) < 3.0 * se
