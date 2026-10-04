"""Pruebas de la ablacion del funcional (diseno §3.1; P§10)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.config.settings11 import FixedDensityConfig
from omega.dynamics.evolution import evolve
from omega.dynamics.fixed_density import evolve_fixed_density
from omega.dynamics.functional import action
from omega.dynamics.gradient import grad_action, lipschitz_bound
from omega.network.weights import from_upper_triangle
from omega.phases.ablation import (
    AblatedParams,
    AblationSpec,
    ablated_action,
    ablated_gradient,
    ablated_lipschitz,
    ablation_set,
    effective_coefficients,
    evolve_ablated,
)
from omega.types import FloatArray, RunStatus


def _w0(n: int, seed: int = 7) -> FloatArray:
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))
    return from_upper_triangle(rng.random(n * (n - 1) // 2), n)


def _jmi(n: int) -> FloatArray:
    return np.ones((n, n)) - np.eye(n)


def _get(n: int, ah: float, gh: float, spec: AblationSpec) -> AblatedParams:
    return next(a for a in ablation_set(ah, gh, n) if a.spec is spec)


def test_effective_coefficients() -> None:
    ap = AblatedParams(0.3, 1.0, 0.2, AblationSpec.FULL)
    for spec, exp in (
        (AblationSpec.FULL, (0.3, 1.0, 0.2)),
        (AblationSpec.NO_TRIANGLES, (0.0, 1.0, 0.2)),
        (AblationSpec.NO_DENSITY, (0.3, 0.0, 0.2)),
        (AblationSpec.NO_DEGREE, (0.3, 1.0, 0.0)),
        (AblationSpec.TRIANGLES_ONLY, (0.3, 0.0, 0.0)),
    ):
        assert effective_coefficients(AblatedParams(0.3, 1.0, 0.2, spec)) == exp
    assert ap.spec is AblationSpec.FULL


def test_params_validation() -> None:
    with pytest.raises(ValueError):
        AblatedParams(0.1, -1.0, 0.0, AblationSpec.FULL)
    with pytest.raises(ValueError):
        AblatedParams(0.1, 1.0, -0.1, AblationSpec.FULL)
    with pytest.raises(ValueError):
        AblatedParams(float("nan"), 1.0, 0.0, AblationSpec.FULL)
    with pytest.raises(TypeError):
        AblatedParams(0.1, 1.0, 0.0, "full")  # type: ignore[arg-type]
    AblatedParams(0.1, 0.0, 0.0, AblationSpec.TRIANGLES_ONLY)  # beta = 0 permitido


def test_full_action_gradient_match_omega10() -> None:
    n = 14
    w = _w0(n)
    p = reduced_to_raw(1.3, 0.7, n)
    ap = AblatedParams(p.alpha, p.beta, p.gamma, AblationSpec.FULL)
    assert abs(ablated_action(w, ap) - action(w, p)) <= 1e-12 * max(1.0, abs(action(w, p)))
    assert float(np.max(np.abs(ablated_gradient(w, ap) - grad_action(w, p)))) <= 1e-12
    assert ablated_lipschitz(n, ap) == pytest.approx(lipschitz_bound(n, p), rel=1e-14)


def test_ablated_gradient_is_derivative_of_action() -> None:
    n, h = 8, 1e-6
    w = (0.2 + 0.5 * _w0(n)) * _jmi(n)
    iu = np.triu_indices(n, k=1)
    for spec in AblationSpec:
        ap = _get(n, 1.5, 1.0, spec)
        g = ablated_gradient(w, ap)
        for e in (0, 5, 17):
            i, j = int(iu[0][e]), int(iu[1][e])
            wp, wm = w.copy(), w.copy()
            wp[i, j] = wp[j, i] = w[i, j] + h
            wm[i, j] = wm[j, i] = w[i, j] - h
            fd = (ablated_action(wp, ap) - ablated_action(wm, ap)) / (2 * h)
            assert fd == pytest.approx(g[i, j], rel=1e-6, abs=1e-8)


def test_lipschitz_per_variant() -> None:
    n = 20
    ap = AblatedParams(0.3, 1.0, 0.2, AblationSpec.FULL)
    assert ablated_lipschitz(n, ap) == pytest.approx(2.0 + 2.0 * (0.2 + 0.3) * 18)
    assert ablated_lipschitz(n, AblatedParams(0.3, 1.0, 0.2, AblationSpec.NO_TRIANGLES)) == pytest.approx(
        2.0 + 2.0 * 0.2 * 18
    )
    assert ablated_lipschitz(n, AblatedParams(0.3, 1.0, 0.2, AblationSpec.TRIANGLES_ONLY)) == pytest.approx(
        2.0 * 0.3 * 18
    )
    with pytest.raises(ValueError):
        ablated_lipschitz(n, AblatedParams(0.0, 1.0, 0.0, AblationSpec.TRIANGLES_ONLY))
    with pytest.raises(ValueError):
        ablated_lipschitz(1, ap)


def test_ablation_set() -> None:
    s = ablation_set(1.5, 1.0, 40)
    assert tuple(a.spec for a in s) == tuple(AblationSpec)
    p = reduced_to_raw(1.5, 1.0, 40)
    assert all((a.alpha, a.beta, a.gamma) == (p.alpha, p.beta, p.gamma) for a in s)
    assert ablation_set(1.5, 1.0, 40, beta=2.0)[0].beta == 2.0


@pytest.mark.parametrize("ah,gh", [(0.5, 0.0), (1.5, 1.0), (2.5, 0.3)])
def test_full_equals_evolve(ah: float, gh: float) -> None:
    n = 16
    w0 = _w0(n)
    p = reduced_to_raw(ah, gh, n)
    cfg = DynamicsConfig(max_steps=3000)
    ref = evolve(w0, p, cfg)
    out = evolve_ablated(w0, AblatedParams(p.alpha, p.beta, p.gamma, AblationSpec.FULL), cfg)
    assert out.status is ref.status and out.steps == ref.steps
    assert out.dt == pytest.approx(ref.dt, rel=1e-14)
    assert float(np.max(np.abs(out.w_final - ref.w_final))) <= 1e-12
    assert np.array_equal(out.snapshot_steps, ref.snapshot_steps)
    assert float(np.max(np.abs(out.snapshots - ref.snapshots))) <= 1e-12
    assert set(out.scalars) == set(ref.scalars)
    for k, v in ref.scalars.items():
        assert float(np.max(np.abs(out.scalars[k] - v))) <= 1e-12 * max(1.0, float(np.max(np.abs(v))))


def test_triangles_only_reaches_complete_graph() -> None:
    n = 12
    ap = _get(n, 1.5, 0.0, AblationSpec.TRIANGLES_ONLY)
    tr = evolve_ablated(_w0(n), ap, DynamicsConfig())
    assert tr.status is RunStatus.CONVERGED
    assert float(np.max(np.abs(tr.w_final - _jmi(n)))) <= 1e-6
    # Ninguna arista decrece (dS/dw = -alpha (W^2)_ab <= 0).
    assert float(np.min(np.diff(tr.snapshots, axis=0))) >= -1e-15


def test_no_triangles_reaches_zero() -> None:
    n = 20
    for gh in (0.0, 1.0):
        ap = _get(n, 3.0, gh, AblationSpec.NO_TRIANGLES)
        tr = evolve_ablated(_w0(n), ap, DynamicsConfig())
        assert tr.status is RunStatus.CONVERGED
        assert float(tr.w_final.max()) <= 1e-6
        assert tr.scalars["action"][-1] <= tr.scalars["action"][0]


def test_no_density_and_no_degree_run() -> None:
    n = 12
    for spec in (AblationSpec.NO_DENSITY, AblationSpec.NO_DEGREE):
        tr = evolve_ablated(_w0(n), _get(n, 1.5, 1.0, spec), DynamicsConfig(max_steps=500))
        assert tr.status in (RunStatus.CONVERGED, RunStatus.MAX_STEPS)
        assert float(tr.w_final.min()) >= 0.0 and float(tr.w_final.max()) <= 1.0


def test_no_triangles_with_fixed_density_gives_uniform() -> None:
    n, rho = 14, 0.25
    for gh in (0.0, 1.0):
        ap = _get(n, 3.0, gh, AblationSpec.NO_TRIANGLES)
        a, b, g = effective_coefficients(ap)
        tr = evolve_fixed_density(
            _w0(n),
            FunctionalParams(alpha=a, beta=b, gamma=g),
            DynamicsConfig(),
            FixedDensityConfig(rho=rho),
        )
        assert tr.status is RunStatus.CONVERGED
        assert float(np.max(np.abs(tr.w_final - rho * _jmi(n)))) <= 1e-6


def test_evolve_ablated_validation() -> None:
    n = 6
    ap = _get(n, 1.0, 0.0, AblationSpec.FULL)
    w0 = 0.3 * _jmi(n)
    with pytest.raises(ValueError):
        evolve_ablated(w0, ap, DynamicsConfig(integrator="sigmoid"))
    with pytest.raises(TypeError):
        evolve_ablated(w0, None, DynamicsConfig())  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        evolve_ablated(w0, ap, None)  # type: ignore[arg-type]
    zero = AblatedParams(0.0, 0.0, 0.0, AblationSpec.FULL)
    with pytest.raises(ValueError):
        evolve_ablated(w0, zero, DynamicsConfig())
    tr = evolve_ablated(w0, zero, DynamicsConfig(dt_mode="fixed", dt_fixed=0.1, max_steps=50))
    assert tr.status is RunStatus.CONVERGED and np.array_equal(tr.w_final, w0)
