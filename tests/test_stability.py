"""Pruebas de stability.py: resumen por semillas, permutacion, nulos, histeresis, F-confirmada."""

from __future__ import annotations

import dataclasses
import math
from typing import Any

import numpy as np
import pytest

from omega.config.seeds import seed_key
from omega.config.settings import DynamicsConfig, OmegaConfig, PhaseThresholds
from omega.experiments.reference_graphs import periodic_lattice
from omega.network.initialization import random_uniform_weights
from omega.network.topology import adjacency, binary_degree
from omega.network.weights import upper_triangle, validate_weight_matrix
from omega.phases.classification import classify_run
from omega.phases.scan import default_config, point_params, simulate
from omega.phases.stability import (
    confirm_geometric,
    degree_preserving_null,
    hysteresis_gap,
    hysteresis_sweep,
    permutation_invariance,
    random_permutation,
    seed_summary,
    shuffled_weight_null,
    wmin_candidate_flags,
)
from omega.types import (
    DimensionEstimate,
    FloatArray,
    GeometryObservables,
    Observables,
    PhaseLabel,
    RunResult,
    RunStatus,
    TopologyObservables,
    Trajectory,
)
from omega.config.settings import FunctionalParams

N = 16
T = PhaseThresholds()


def cfg_small() -> OmegaConfig:
    cfg = default_config(N, master_entropy=5, experiment_id=5, replicates=3)
    return dataclasses.replace(cfg, dynamics=DynamicsConfig(max_steps=20_000))


def rng(seed: int = 1) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


# ---------------------------------------------------------------- resultados sinteticos
_Z = np.zeros(0, dtype=np.float64)


def _dim(v: float) -> DimensionEstimate:
    return DimensionEstimate(v, 0.0, "ok", (0, 2), _Z, _Z, _Z, True, "shell")


def fake_run(d_eff: float, candidate: bool = True, alpha_hat: float = 1.0) -> RunResult:
    """RunResult sintetico: F-candidata (D_s=D_eff) o C (D_eff sin ventana)."""
    top = TopologyObservables(
        n=100, mean_strength=5.0, std_strength=0.2, cv_strength=0.04, mean_binary_degree=6.0,
        binary_density=0.03, mean_weight=0.03, giant_size=100, giant_fraction=1.0, n_components=1,
        n_large_components=1, clustering_weighted=0.05, clustering_binary=0.05, frac_at_zero=0.9,
        frac_at_one=0.0,
    )
    nw = DimensionEstimate(math.nan, math.nan, "no_window", None, _Z, _Z, _Z, False, "shell")
    geo = GeometryObservables(5.0, 5.0, 9.0, _dim(d_eff) if candidate else nw, nw, nw,
                              _dim(d_eff) if candidate else nw)
    o = Observables(top, geo)
    traj = Trajectory(np.zeros((2, 2)), RunStatus.CONVERGED, 0, 0.0, 0.0, {}, np.zeros(1, dtype=np.int64),
                      np.zeros((1, 1)))
    return RunResult(seed=seed_key(default_config(10).seeds, 0, 0), params=FunctionalParams(alpha=0.0),
                     alpha_hat=alpha_hat, gamma_hat=0.0, w0=np.zeros((2, 2)), trajectory=traj,
                     observables=o, assessment=classify_run(o, RunStatus.CONVERGED, T))


def deff_key(r: RunResult) -> float:
    return r.observables.geometry.d_eff.value


# ---------------------------------------------------------------- seed_summary
def test_seed_summary() -> None:
    rs = [fake_run(2.9), fake_run(3.1), fake_run(3.0)]
    m, s, n = seed_summary(rs, deff_key)
    assert (m, n) == (pytest.approx(3.0), 3) and s == pytest.approx(0.1)
    m1, s1, n1 = seed_summary([fake_run(3.0), fake_run(0.0, candidate=False)], deff_key)
    assert (m1, s1, n1) == (3.0, 0.0, 1)  # NaN excluido
    mm, ss, nn = seed_summary([], deff_key)
    assert math.isnan(mm) and math.isnan(ss) and nn == 0


# ---------------------------------------------------------------- permutacion
@pytest.mark.parametrize("kind", ["uniform", "ring", "evolved"])
def test_permutation_invariance_exact(kind: str) -> None:
    cfg = cfg_small()
    r = rng(3)
    w: FloatArray
    if kind == "uniform":
        w = random_uniform_weights(N, r)
    elif kind == "ring":
        w = periodic_lattice((N,))
    else:
        w = simulate(cfg, point_params(1.0, 1.0, N), seed_key(cfg.seeds, 0, 0)).trajectory.w_final
    perm = random_permutation(N, r)
    out = permutation_invariance(w, perm, cfg)
    assert out["status_mismatches"] == 0.0
    assert out["max_abs_diff"] <= 1e-9
    assert out["topology.mean_strength"] <= 1e-9


def test_permutation_invariance_larger_uniform_and_input_unchanged() -> None:
    cfg = dataclasses.replace(cfg_small(), init=dataclasses.replace(cfg_small().init, n=60))
    w = random_uniform_weights(60, rng(8))
    w0 = w.copy()
    out = permutation_invariance(w, random_permutation(60, rng(9)), cfg)
    assert out["max_abs_diff"] <= 1e-9
    assert np.array_equal(w, w0)
    with pytest.raises(ValueError):
        permutation_invariance(w, np.zeros(60, dtype=np.int64), cfg)


# ---------------------------------------------------------------- nulos
def test_shuffled_weight_null_conserves_weights() -> None:
    w = random_uniform_weights(30, rng(1))
    w0 = w.copy()
    nl = shuffled_weight_null(w, 0.1, rng(2))
    validate_weight_matrix(nl)
    assert np.array_equal(np.sort(upper_triangle(nl)), np.sort(upper_triangle(w)))
    assert not np.array_equal(nl, w)
    assert np.array_equal(w, w0)
    assert np.array_equal(nl, shuffled_weight_null(w, 0.1, rng(2)))  # determinista
    with pytest.raises(ValueError):
        shuffled_weight_null(w, 1.5, rng(2))


def _structured(n: int = 30) -> FloatArray:
    r = rng(4)
    w = random_uniform_weights(n, r)
    return np.where(w > 0.7, w, 0.0)  # grafo disperso con pesos variables


def test_degree_preserving_null_conserves_degrees_and_weights() -> None:
    w = _structured()
    w0 = w.copy()
    nl = degree_preserving_null(w, 0.1, 2000, rng(5))
    validate_weight_matrix(nl)
    a0, a1 = adjacency(w, 0.1), adjacency(nl, 0.1)
    assert np.array_equal(binary_degree(a0), binary_degree(a1))
    assert np.array_equal(np.sort(upper_triangle(nl)), np.sort(upper_triangle(w)))
    assert not np.array_equal(a0, a1)  # el recableado cambio la topologia
    assert np.array_equal(w, w0)
    assert np.array_equal(nl, degree_preserving_null(w, 0.1, 2000, rng(5)))
    assert np.array_equal(degree_preserving_null(w, 0.1, 0, rng(5)), w)
    with pytest.raises(ValueError):
        degree_preserving_null(w, 0.1, -1, rng(5))


def test_degree_preserving_null_degenerate() -> None:
    z = np.zeros((6, 6))
    assert np.array_equal(degree_preserving_null(z, 0.1, 50, rng(1)), z)


# ---------------------------------------------------------------- histeresis
def test_hysteresis_sweep_shows_loop() -> None:
    cfg = cfg_small()
    path = (0.0, 1.5, 3.0)
    w_rand = random_uniform_weights(N, rng(6))
    fwd = hysteresis_sweep(w_rand, path, 0.0, cfg)
    ones = np.ones((N, N)) - np.eye(N)
    bwd = hysteresis_sweep(ones, path[::-1], 0.0, cfg)
    assert [r.alpha_hat for r in fwd] == pytest.approx(list(path))
    assert [r.alpha_hat for r in bwd] == pytest.approx(list(path[::-1]))
    assert [r.assessment.label for r in fwd] == [PhaseLabel.A] * 3  # W=0 es punto fijo
    assert [r.assessment.label for r in bwd] == [PhaseLabel.E, PhaseLabel.E, PhaseLabel.A]
    gap = hysteresis_gap(fwd, bwd, lambda r: r.observables.topology.mean_weight)
    assert gap.shape == (3,)
    assert gap == pytest.approx([0.0, -1.0, -1.0], abs=1e-6)
    # continuacion: el estado inicial del 2.o paso es el final del 1.o
    assert np.array_equal(fwd[1].w0, fwd[0].trajectory.w_final)
    assert np.array_equal(bwd[1].w0, bwd[0].trajectory.w_final)


def test_hysteresis_validation() -> None:
    cfg = cfg_small()
    with pytest.raises(ValueError):
        hysteresis_sweep(np.zeros((N, N)), (), 0.0, cfg)
    a, b = fake_run(3.0, alpha_hat=1.0), fake_run(3.0, alpha_hat=2.0)
    with pytest.raises(ValueError):
        hysteresis_gap([a], [b], deff_key)
    with pytest.raises(ValueError):
        hysteresis_gap([a, a], [a], deff_key)


# ---------------------------------------------------------------- confirm_geometric
def _group(vals: list[float], cand: bool = True) -> list[RunResult]:
    return [fake_run(v, cand) for v in vals]


def _confirm(**over: Any) -> Any:
    args: dict[str, Any] = dict(
        point=_group([3.0, 3.05, 2.95, 3.0, 3.1]),
        by_size={100: _group([3.0, 3.0]), 200: _group([3.1, 3.1]), 300: _group([3.0, 3.0])},
        neighbors=[_group([3.0, 3.0]), _group([3.1, 3.0])],
        t=T,
        permutation_max_diff=1e-12,
        wmin_flags=(True, True, True, True, False),
        null_candidates=(False, False, False),
    )
    args.update(over)
    return confirm_geometric(**args)


def test_confirm_geometric_all_criteria() -> None:
    a = _confirm()
    assert a.label is PhaseLabel.F and a.flags["F_confirmed"]
    assert all(a.flags.values())


@pytest.mark.parametrize(
    ("over", "failed"),
    [
        ({"point": _group([3.0, 3.0, 3.0, 3.0, 3.0], True)[:3] + _group([3.0, 3.0], False)}, "seed_fraction"),
        ({"point": _group([2.0, 3.0, 4.0, 2.5, 3.5])}, "seed_std"),
        ({"by_size": {100: _group([3.0]), 200: _group([3.5])}}, "size_stability"),
        ({"by_size": {100: _group([2.9]), 200: _group([3.0]), 300: _group([3.1])}}, "size_stability"),  # deriva
        ({"by_size": {100: _group([3.0])}}, "size_stability"),
        ({"neighbors": [_group([3.0, 3.0]), _group([3.0, 3.0], False)]}, "param_stability"),
        ({"neighbors": []}, "param_stability"),
        ({"permutation_max_diff": 1e-6}, "permutation"),
        ({"permutation_max_diff": None}, "permutation"),
        ({"wmin_flags": (True, True, True, False, False)}, "wmin_stability"),
        ({"wmin_flags": None}, "wmin_stability"),
        ({"null_candidates": (True, True, False)}, "null_distinct"),
        ({"null_candidates": None}, "null_distinct"),
    ],
)
def test_confirm_geometric_each_criterion_can_fail(over: dict[str, Any], failed: str) -> None:
    a = _confirm(**over)
    assert a.flags[failed] is False
    assert a.flags["F_confirmed"] is False
    assert a.label is not PhaseLabel.F


def test_confirm_geometric_label_when_not_confirmed() -> None:
    a = _confirm(permutation_max_diff=None)
    assert a.label is PhaseLabel.U  # F-candidata no confirmada
    only_c = _confirm(point=_group([3.0] * 4, False))
    assert only_c.label is PhaseLabel.C
    with pytest.raises(ValueError):
        confirm_geometric([], {}, [], T)


def test_wmin_candidate_flags_shape() -> None:
    cfg = cfg_small()
    w = random_uniform_weights(N, rng(2))
    f = wmin_candidate_flags(w, cfg)
    assert len(f) == len(cfg.graph.w_min_sensitivity) == 5
    assert not any(f)  # U(0,1) es hiperdensa/sin ventana: nunca F-candidata
