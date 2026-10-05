"""L-3a (Consejo Omega Rev. 2, R2.5): lemas KKT de S0 y de Omega-B verificados en codigo.

Gradiente: G = -alpha W^2 + 2 beta W + 2 gamma (u_i + u_j) (u = k - k_medio). Alcance: demuestra
condiciones necesarias de primer orden sobre estados concretos; no afirma nada sobre minimos globales.
"""

from __future__ import annotations

from math import comb

import numpy as np
import pytest

from omega.config.convert import reduced_to_raw
from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.config.settings11 import FixedDensityConfig
from omega.dynamics.evolution import evolve, step_scalars
from omega.dynamics.fixed_density import evolve_fixed_density, uniform_state_threshold
from omega.dynamics.gradient import grad_action
from omega.experiments.reference_graphs import periodic_lattice
from omega.controls.random_geometric import rgg_torus
from omega.landscape.kkt import KKTResiduals, codegree, kkt_residuals_fixed_density, kkt_residuals_s0, open_p3_count
from omega.network.initialization import random_uniform_weights
from omega.types import FloatArray, RunStatus

GRID_ALPHA_S0 = (0.5, 1.0, 1.5, 1.9, 2.1, 2.5, 3.0)
GRID_GAMMA = (0.0, 1.0, 10.0, 100.0)
SUPPORT_THR = 1e-6
CFG = DynamicsConfig(snapshot_schedule="none")


def _rng(*key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(20261005, spawn_key=key)))


def _mixed_state(n: int, seed: int) -> FloatArray:
    """Estado con ceros y unos exactos y entradas interiores."""
    w = random_uniform_weights(n, _rng(9, seed))
    w[w < 0.3] = 0.0
    w[w > 0.8] = 1.0
    return w


# ------------------------------------------------------------------ punto 1


def _brute_s0(w: FloatArray, p: FunctionalParams, tol: float) -> tuple[float, float, float]:
    """Referencia por bucles, independiente de kkt.py, a partir de G = grad_action (la verdad del codigo)."""
    g = grad_action(w, p)
    n = w.shape[0]
    v0 = v1 = vi = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            x, gij = w[i, j], g[i, j]
            if x <= tol:
                v0 = max(v0, -gij)
            elif x >= 1.0 - tol:
                v1 = max(v1, gij)
            else:
                vi = max(vi, abs(gij))
    return v0, v1, vi


def test_residuals_match_gradient_conditions_s0() -> None:
    for seed in range(4):
        w = _mixed_state(14, seed)
        for ah, gh in ((0.5, 0.0), (2.0, 3.0), (4.0, 20.0)):
            p = reduced_to_raw(ah, gh, 14)
            r = kkt_residuals_s0(w, p, tol=1e-6)
            v0, v1, vi = _brute_s0(w, p, 1e-6)
            assert r.max_violation_zero == pytest.approx(max(0.0, v0), abs=1e-13)
            assert r.max_violation_one == pytest.approx(max(0.0, v1), abs=1e-13)
            assert r.max_violation_interior == pytest.approx(vi, abs=1e-13)
            assert r.n_zero + r.n_one + r.n_interior == 14 * 13 // 2
            # coincide con la norma de gradiente proyectado del integrador (R5) cuando tol = 0
            sc = step_scalars(w, w, p)
            assert kkt_residuals_s0(w, p, tol=0.0).max_violation == pytest.approx(sc["grad_proj_inf"], abs=1e-13)


def test_multiplier_matches_code_multiplier_omega_b() -> None:
    n = 30
    for gh, f in ((0.0, 0.5), (10.0, 0.5), (10.0, 1.5)):
        rho = 0.3
        p = reduced_to_raw(f * uniform_state_threshold(gh, rho, n), gh, n)
        tr = evolve_fixed_density(random_uniform_weights(n, _rng(8, int(gh), int(10 * f))), p, CFG,
                                  FixedDensityConfig(rho=rho))
        assert tr.status is RunStatus.CONVERGED
        r = kkt_residuals_fixed_density(tr.w_final, p)
        assert r.multiplier is not None
        assert r.multiplier == pytest.approx(float(tr.scalars["multiplier"][-1]), rel=1e-9, abs=1e-9)
        assert r.relative_violation <= 1e-7


def test_multiplier_without_interior_minimises_violation() -> None:
    # K_m disjuntas en N=12 (sin interiores): G+lam factible en [a,b]; con lam=punto medio la violacion es (a-b)/2
    n = 12
    w = np.zeros((n, n))
    w[:6, :6] = 1.0
    w[6:, 6:] = 1.0
    np.fill_diagonal(w, 0.0)
    p = FunctionalParams(alpha=0.9, beta=1.0)  # unos: G = 2-0.9*4 = -1.6 ; ceros: G = 0
    r = kkt_residuals_fixed_density(w, p)
    assert r.n_interior == 0 and r.multiplier is not None
    assert r.max_violation == 0.0 and r.multiplier == pytest.approx(0.8)  # punto medio de [0, 1.6]
    p2 = FunctionalParams(alpha=0.2, beta=1.0)  # unos: G = 2 - 0.8 = 1.2 > 0 (viola); ceros: G = 0
    r2 = kkt_residuals_fixed_density(w, p2)
    assert r2.max_violation == pytest.approx(0.6)  # lam = 0.5*(0 + (-1.2)) = -0.6


# ------------------------------------------------------------------ punto 2 (P-KKT) y 5 (lema de escala)


def _pkkt_finals() -> list[tuple[float, int, FloatArray, FunctionalParams, RunStatus]]:
    out = []
    n = 30
    for ia, ah in enumerate((0.5, 1.5, 2.5, 4.0)):
        p = reduced_to_raw(ah, 0.0, n)
        for s in range(20):
            tr = evolve(random_uniform_weights(n, _rng(ia, s)), p, CFG)
            out.append((ah, s, tr.w_final, p, tr.status))
    return out


@pytest.fixture(scope="module")
def pkkt_finals() -> list[tuple[float, int, FloatArray, FunctionalParams, RunStatus]]:
    return _pkkt_finals()


def test_p_kkt_s0_gamma_zero_converged_finals(
    pkkt_finals: list[tuple[float, int, FloatArray, FunctionalParams, RunStatus]],
) -> None:
    assert len(pkkt_finals) == 80
    for ah, s, w, p, status in pkkt_finals:
        assert status is RunStatus.CONVERGED, (ah, s)
        assert open_p3_count(w, SUPPORT_THR) == 0, (ah, s)
        assert kkt_residuals_s0(w, p).relative_violation <= 1e-7, (ah, s)


def test_p_kkt_nontrivial_union_of_cliques() -> None:
    """Contenido no vacio de P-KKT/lema de escala: K_10 x 3 (N=30, gamma=0) es KKT sii alpha_hat >= (N-2)/(m-2)."""
    n, m = 30, 10
    w = np.zeros((n, n))
    for b in range(3):
        w[b * m:(b + 1) * m, b * m:(b + 1) * m] = 1.0
    np.fill_diagonal(w, 0.0)
    assert open_p3_count(w, SUPPORT_THR) == 0
    ah_min = (n - 2) / (m - 2)  # = 3.5
    assert kkt_residuals_s0(w, reduced_to_raw(ah_min + 0.5, 0.0, n)).max_violation == 0.0
    assert kkt_residuals_s0(w, reduced_to_raw(ah_min - 0.5, 0.0, n)).max_violation > 0.0


def test_scale_lemma_on_converged_finals(
    pkkt_finals: list[tuple[float, int, FloatArray, FunctionalParams, RunStatus]],
) -> None:
    """Toda arista con W=1: alpha c_ij >= 2 beta + 2 gamma (u_i+u_j) - tol (finales del punto 2)."""
    n_ones = 0
    for ah, s, w, p, status in pkkt_finals:
        assert status is RunStatus.CONVERGED
        ones = np.triu(w == 1.0, 1)
        if not ones.any():
            continue
        k = w.sum(axis=1)
        u = k - k.mean()
        rhs = 2.0 * p.beta + 2.0 * p.gamma * (u[:, None] + u[None, :])
        lhs = p.alpha * codegree(w)
        n_ones += int(ones.sum())
        assert np.all(lhs[ones] >= rhs[ones] - 1e-6), (ah, s, float(np.min((lhs - rhs)[ones])))
    assert n_ones > 0


# ------------------------------------------------------------------ puntos 3 y 4


@pytest.mark.parametrize("side", [6, 8])
def test_p_kkt_reg_torus_binary_never_kkt_in_s0(side: int) -> None:
    w = periodic_lattice((side, side, side))
    n = w.shape[0]
    for ah in GRID_ALPHA_S0:
        for gh in GRID_GAMMA:
            r = kkt_residuals_s0(w, reduced_to_raw(ah, gh, n))
            assert r.max_violation > 0.0, (side, ah, gh)


@pytest.mark.parametrize("n", [216, 512])
def test_rgg3_k12_binary_never_kkt_in_s0(n: int) -> None:
    for s in range(3):
        w = rgg_torus(n, 3, 12.0, _rng(7, n, s), weights="binary")
        for ah in GRID_ALPHA_S0:
            for gh in GRID_GAMMA:
                r = kkt_residuals_s0(w, reduced_to_raw(ah, gh, n))
                assert r.max_violation > 0.0, (n, s, ah, gh)


# ------------------------------------------------------------------ punto 6


def _k22_halo(rho: float) -> tuple[FloatArray, float]:
    n, k = 100, 22
    w0 = rho * comb(n, 2)
    h = (w0 - comb(k, 2)) / comb(n - k, 2)
    w = np.zeros((n, n))
    w[:k, :k] = 1.0
    w[k:, k:] = h
    np.fill_diagonal(w, 0.0)
    return w, h


@pytest.mark.parametrize("rho", [0.05, 0.1])
def test_k22_halo_is_kkt_with_gamma_and_not_without(rho: float) -> None:
    """rho=0.05 reproduce el halo h~0.0055 de O-04; alpha_hat = 0.5*alpha_c(gamma_hat) con su propio gamma_hat."""
    n = 100
    w, h = _k22_halo(rho)
    assert float(w.sum() / 2.0) == pytest.approx(rho * comb(n, 2), rel=1e-12)
    assert 0.0 < h < 0.5

    def res(gh_state: float, gh: float) -> KKTResiduals:
        p = reduced_to_raw(0.5 * uniform_state_threshold(gh_state, rho, n), gh, n)
        return kkt_residuals_fixed_density(w, p)

    on = res(10.0, 10.0)
    off = res(0.0, 0.0)
    assert on.relative_violation <= 1e-7
    assert on.multiplier is not None and on.multiplier > 0.0
    assert off.relative_violation > 1e-7
    assert off.max_violation_zero > 0.0


@pytest.mark.parametrize(("rho", "kkt_without_gamma"), [(0.05, False), (0.1, True)])
def test_k22_halo_gamma_role_at_fixed_alpha(rho: float, kkt_without_gamma: bool) -> None:
    """Matiz registrado (resultados del bloque L): con el MISMO alpha_hat = 0.5*alpha_c(10) y gamma_hat=0, el halo
    de O-04 (rho=0.05, h~0.0055) deja de ser KKT (gamma es necesario), pero con rho=0.1 (h~0.088) sigue siendo KKT.
    "Halo KKT exige gamma>0" no es una afirmacion general: depende de rho y del nivel de alpha_hat."""
    n = 100
    w, _ = _k22_halo(rho)
    p = reduced_to_raw(0.5 * uniform_state_threshold(10.0, rho, n), 0.0, n)
    assert (kkt_residuals_fixed_density(w, p).relative_violation <= 1e-7) is kkt_without_gamma


def test_binary_torus_not_kkt_in_omega_b() -> None:
    w = periodic_lattice((6, 6, 6))
    n = w.shape[0]
    for ah in (0.5, 1.0, 2.0, 3.0):
        for gh in GRID_GAMMA:
            r = kkt_residuals_fixed_density(w, reduced_to_raw(ah, gh, n))
            assert r.max_violation > 0.0, (ah, gh)
