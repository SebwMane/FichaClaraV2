"""Pruebas de la rama Omega-B de densidad fija (diseno §1.9, §3.1; P§9)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.config.settings11 import FixedDensityConfig
from omega.dynamics.evolution import step_size
from omega.dynamics.fixed_density import (
    density_to_total,
    evolve_fixed_density,
    fixed_density_step,
    project_capped_simplex,
    project_fixed_density,
    uniform_state_threshold,
    uniform_state_unstable,
)
from omega.dynamics.gradient import grad_action
from omega.network.weights import from_upper_triangle, upper_triangle
from omega.types import FloatArray, RunStatus


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


def _jmi(n: int) -> FloatArray:
    return np.ones((n, n)) - np.eye(n)


def _perturbed(n: int, rho: float, rng: np.random.Generator, amp: float = 1e-3) -> FloatArray:
    x = rng.standard_normal(n * (n - 1) // 2)
    x -= x.mean()
    x *= amp / float(np.abs(x).max())
    return from_upper_triangle(rho + x, n)


def test_density_to_total() -> None:
    assert density_to_total(0.1, 40) == pytest.approx(78.0)
    assert density_to_total(0.0, 5) == 0.0 and density_to_total(1.0, 5) == 10.0
    with pytest.raises(ValueError):
        density_to_total(1.5, 10)
    with pytest.raises(ValueError):
        density_to_total(0.5, 1)
    with pytest.raises(TypeError):
        density_to_total(0.5, 10.0)  # type: ignore[arg-type]


def test_projection_feasibility_and_optimality() -> None:
    rng = _rng(1)
    m = 60
    for scale in (0.3, 1.0, 3.0):
        v = scale * rng.standard_normal(m) + 0.4
        for total in (0.0, 1.0, 12.5, 30.0, 59.0, float(m)):
            x = project_capped_simplex(v, total)
            assert abs(float(x.sum()) - total) <= 1e-9 * max(1.0, total)
            assert float(x.min()) >= 0.0 and float(x.max()) <= 1.0
            d0 = float(np.linalg.norm(x - v))
            for _ in range(200):
                y = rng.random(m)
                dirn = y - y.mean()
                c = np.full(m, total / m)
                with np.errstate(divide="ignore"):
                    pos = np.where(dirn > 0, (1.0 - c) / np.where(dirn > 0, dirn, 1.0), np.inf)
                    neg = np.where(dirn < 0, (0.0 - c) / np.where(dirn < 0, dirn, 1.0), np.inf)
                lam_max = float(min(pos.min(), neg.min()))
                z = c + rng.random() * lam_max * dirn if np.isfinite(lam_max) else c
                assert d0 <= float(np.linalg.norm(z - v)) + 1e-12


def test_projection_is_clip_with_common_shift_and_idempotent() -> None:
    rng = _rng(2)
    v = rng.standard_normal(100) * 0.5 + 0.3
    x = project_capped_simplex(v, 25.0)
    free = (x > 0.0) & (x < 1.0)
    shifts = (v - x)[free]
    assert float(shifts.max() - shifts.min()) <= 1e-12
    assert np.allclose(project_capped_simplex(x, 25.0), x, atol=1e-12)
    inside = np.full(10, 0.3)
    assert np.allclose(project_capped_simplex(inside, 3.0), inside, atol=1e-12)


def test_projection_validation() -> None:
    with pytest.raises(ValueError):
        project_capped_simplex(np.zeros(5), -0.1)
    with pytest.raises(ValueError):
        project_capped_simplex(np.zeros(5), 5.1)
    with pytest.raises(ValueError):
        project_capped_simplex(np.array([0.0, np.nan]), 1.0)
    with pytest.raises(ValueError):
        project_capped_simplex(np.zeros((2, 2)), 1.0)
    with pytest.raises(ValueError):
        project_fixed_density(_jmi(5), 11.0)


def test_project_fixed_density_matrix() -> None:
    n = 12
    w = _perturbed(n, 0.3, _rng(3), 0.2)
    out = project_fixed_density(w, 15.0)
    assert np.array_equal(out, out.T) and np.all(np.diag(out) == 0.0)
    assert abs(float(upper_triangle(out).sum()) - 15.0) <= 1e-9 * 15.0
    assert float(out.min()) >= 0.0 and float(out.max()) <= 1.0


def test_uniform_state_is_fixed_point() -> None:
    n, rho = 40, 0.1
    w = rho * _jmi(n)
    total = density_to_total(rho, n)
    for ah, gh in ((0.0, 0.0), (5.0, 0.0), (15.0, 1.0)):
        p = reduced_to_raw(ah, gh, n)
        out = fixed_density_step(w, grad_action(w, p), 0.01, total)
        assert float(np.max(np.abs(out - w))) <= 1e-12


def test_threshold_formula() -> None:
    assert uniform_state_threshold(0.0, 0.1, 40) == pytest.approx(38.0 / (0.1 * 36.0))
    assert uniform_state_threshold(0.0, 0.1, 40) == pytest.approx(10.556, abs=1e-3)
    assert uniform_state_threshold(1.0, 0.1, 40) == pytest.approx(2.0 * uniform_state_threshold(0.0, 0.1, 40))
    ac = uniform_state_threshold(0.0, 0.1, 40)
    assert uniform_state_unstable(1.01 * ac, 0.0, 0.1, 40)
    assert not uniform_state_unstable(0.99 * ac, 0.0, 0.1, 40)
    assert not uniform_state_unstable(ac, 0.0, 0.1, 40)
    with pytest.raises(ValueError):
        uniform_state_threshold(0.0, 0.1, 4)
    with pytest.raises(ValueError):
        uniform_state_threshold(-1.0, 0.1, 40)
    with pytest.raises(ValueError):
        uniform_state_threshold(0.0, 0.0, 40)


@pytest.mark.parametrize("gh", [0.0, 1.0])
def test_threshold_dynamics(gh: float) -> None:
    n, rho = 40, 0.1
    ac = uniform_state_threshold(gh, rho, n)
    if gh == 1.0:
        assert ac == pytest.approx(2.0 * uniform_state_threshold(0.0, rho, n))
    cfg = DynamicsConfig(snapshot_schedule="none")
    fd = FixedDensityConfig(rho=rho)
    below = evolve_fixed_density(_perturbed(n, rho, _rng(4)), reduced_to_raw(0.9 * ac, gh, n), cfg, fd)
    assert below.status is RunStatus.CONVERGED
    assert float(np.max(np.abs(below.w_final - rho * _jmi(n)))) < 1e-6
    above = evolve_fixed_density(_perturbed(n, rho, _rng(4)), reduced_to_raw(1.3 * ac, gh, n), cfg, fd)
    assert float(np.max(np.abs(above.w_final - rho * _jmi(n)))) > 0.1


def test_single_clique_above_threshold_without_gamma() -> None:
    n, rho = 40, 0.1
    ac = uniform_state_threshold(0.0, rho, n)
    p = reduced_to_raw(2.0 * ac, 0.0, n)
    tr = evolve_fixed_density(
        _perturbed(n, rho, _rng(5)), p, DynamicsConfig(snapshot_schedule="none"), FixedDensityConfig(rho=rho)
    )
    assert tr.status is RunStatus.CONVERGED
    w = tr.w_final
    nodes = np.where((w > 0.5).sum(axis=1) > 0)[0]
    rest = np.setdiff1d(np.arange(n), nodes)
    assert nodes.size >= 3
    inner = w[np.ix_(nodes, nodes)]
    off = ~np.eye(nodes.size, dtype=bool)
    assert float(inner[off].mean()) >= 0.99
    assert float(w[np.ix_(rest, np.arange(n))].max(initial=0.0)) <= 1e-6  # resto aislado
    assert abs(float(upper_triangle(w).sum()) - density_to_total(rho, n)) <= 1e-9 * 78.0


def test_evolve_structure_and_conservation() -> None:
    n, rho = 12, 0.3
    p = reduced_to_raw(2.0, 0.5, n)
    cfg = DynamicsConfig(max_steps=200)
    w0 = from_upper_triangle(_rng(6).random(n * (n - 1) // 2), n)
    tr = evolve_fixed_density(w0, p, cfg, FixedDensityConfig(rho=rho))
    total = density_to_total(rho, n)
    assert tr.dt == step_size(n, p, cfg)
    assert "multiplier" in tr.scalars and "action" in tr.scalars and "max_dw" in tr.scalars
    assert all(len(v) == tr.steps + 1 for v in tr.scalars.values())
    assert float(tr.scalars["max_dw"][0]) == 0.0
    sums = tr.snapshots.sum(axis=1)
    assert float(np.max(np.abs(sums - total))) <= 1e-9 * total
    assert tr.snapshots.shape[0] == tr.snapshot_steps.shape[0]
    # rho=None: densidad = media de W0.
    tr2 = evolve_fixed_density(w0, p, cfg, FixedDensityConfig())
    assert abs(float(upper_triangle(tr2.w_final).sum()) - float(upper_triangle(w0).sum())) <= 1e-9 * 66.0


def test_evolve_validation() -> None:
    n = 6
    p = reduced_to_raw(1.0, 0.0, n)
    w0 = 0.3 * _jmi(n)
    fd = FixedDensityConfig(rho=0.3)
    with pytest.raises(ValueError):
        evolve_fixed_density(w0, p, DynamicsConfig(integrator="sigmoid"), fd)
    with pytest.raises(TypeError):
        evolve_fixed_density(w0, p, DynamicsConfig(), None)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        evolve_fixed_density(w0, None, DynamicsConfig(), fd)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        fixed_density_step(w0, np.zeros((n, n)), 0.0, 4.0)
    with pytest.raises(ValueError):
        fixed_density_step(w0, np.zeros((n + 1, n + 1)), 0.1, 4.0)
    assert FunctionalParams(alpha=0.0).alpha == 0.0
