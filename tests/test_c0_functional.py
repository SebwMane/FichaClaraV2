"""Omega-C0: funcional, gradiente, cota inferior y estados que la alcanzan (C0-T1/T2)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.c0.functional import C0Params, action_and_grad, action_c0, grad_c0, kkt_box, params_from_targets
from omega.c0.references import (
    complete_bipartite_half,
    clique_union,
    optimal_amplitude,
    torus_lattice_3d,
)

N = 216


def _rand_w(n: int, density: float, rng: np.random.Generator) -> np.ndarray:
    u = np.triu(rng.random((n, n)) * (rng.random((n, n)) < density), k=1)
    return np.asarray(u + u.T, dtype=np.float64)


@pytest.mark.parametrize("c_star,k_star,a", [(0, 6, 1.0), (2, 8, 0.25), (4, 5, 4.0), (8, 12, 1.0)])
def test_gradient_matches_finite_differences(c_star: int, k_star: float, a: float) -> None:
    rng = np.random.Generator(np.random.PCG64(3))
    p = params_from_targets(c_star, k_star, a)
    w = _rand_w(12, 0.6, rng)
    g = grad_c0(w, p)
    s0, g2 = action_and_grad(w, p)
    assert np.allclose(g, g2)
    assert s0 == pytest.approx(action_c0(w, p))
    h = 1e-6
    for _ in range(15):
        i, j = rng.choice(12, size=2, replace=False)
        wp = w.copy()
        wm = w.copy()
        wp[i, j] = wp[j, i] = w[i, j] + h
        wm[i, j] = wm[j, i] = w[i, j] - h
        fd = (action_c0(wp, p) - action_c0(wm, p)) / (2 * h)
        assert fd == pytest.approx(g[i, j], rel=1e-5, abs=1e-6)
    assert np.all(np.diag(g) == 0.0) and np.allclose(g, g.T)


def test_params_from_targets() -> None:
    p = params_from_targets(0, 6, 1.0)
    assert p.lam == 2.0 and p.c_star == 0.0 and p.psi_star == pytest.approx(1.0) and p.k_star == pytest.approx(6.0)
    q = params_from_targets(4, 5, 1.0)
    assert q.lam == pytest.approx(0.2) and q.c_star == pytest.approx(4.0) and q.k_star == pytest.approx(5.0)
    assert q.lower_bound(N) == pytest.approx(-N * q.psi_star**2 / (16 * q.kappa))
    with pytest.raises(ValueError):
        params_from_targets(2, 6, -5.0)  # psi* <= 0


@pytest.mark.parametrize("c_star,k_star,a", [(0, 6, 1.0), (2, 8, 0.25), (8, 4, 4.0)])
def test_lower_bound_never_violated(c_star: int, k_star: float, a: float) -> None:
    rng = np.random.Generator(np.random.PCG64(11))
    p = params_from_targets(c_star, k_star, a)
    lb = p.lower_bound(40)
    for dens in (0.02, 0.1, 0.3, 0.6, 1.0):
        for _ in range(40):
            w = _rand_w(40, dens, rng)
            assert action_c0(w, p) >= lb - 1e-9


def test_bipartite_and_t3_attain_bound_at_favorable_cell() -> None:
    p = params_from_targets(0, 6, 1.0)
    lb = p.lower_bound(N)
    t, s = optimal_amplitude(complete_bipartite_half(N), p)
    assert abs(s - lb) <= 1e-6 * abs(lb)
    assert t == pytest.approx(2 * 6 / N, rel=1e-3)
    t3 = torus_lattice_3d(6)
    assert action_c0(t3, p) == pytest.approx(lb, rel=1e-9)
    assert abs(action_c0(t3, p)) == pytest.approx(324.0, rel=1e-9)


def test_clique_union_attains_bound_when_kstar_is_cstar_plus_one() -> None:
    p = params_from_targets(4, 5, 1.0)
    w = clique_union(N, 6)
    assert action_c0(w, p) == pytest.approx(p.lower_bound(N), rel=1e-9)


def test_kkt_box_t3_favorable_cell() -> None:
    p = params_from_targets(0, 6, 1.0)
    r = kkt_box(torus_lattice_3d(6), p)
    assert r["is_kkt"]  # G = 0 en las aristas (degenerado): KKT no estricto
    assert r["n_upper"] == 648
    assert not kkt_box(np.full((N, N), 0.3) - 0.3 * np.eye(N), p)["is_kkt"]
