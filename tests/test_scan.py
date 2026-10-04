"""Pruebas de scan.py: malla, parametros, simulate, barrido y probabilidades de fase."""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

from omega.config.convert import raw_to_reduced
from omega.config.seeds import SeedKey, seed_key
from omega.config.settings import DynamicsConfig, FunctionalParams, OmegaConfig, ScanConfig
from omega.phases.scan import (
    default_config,
    observe,
    parameter_grid,
    phase_probabilities,
    point_params,
    run_summary,
    scan_phase_diagram,
    simulate,
    simulate_from,
    static_run,
    to_jsonable,
    with_experiment,
    with_params,
)
from omega.statistics.summary import proportion_ci
from omega.types import PhaseLabel, RunStatus

N = 16


def small_cfg(scan: ScanConfig | None = None, replicates: int = 2) -> OmegaConfig:
    cfg = default_config(N, master_entropy=99, experiment_id=3, replicates=replicates)
    return dataclasses.replace(cfg, dynamics=DynamicsConfig(max_steps=20_000), scan=scan)


def test_parameter_grid_order_and_indices() -> None:
    scan = ScanConfig(alpha_hat=(0.0, 1.0, 2.0), gamma_hat=(0.0, 5.0))
    g = parameter_grid(scan)
    assert len(g) == 6
    assert [x[0] for x in g] == list(range(6))
    assert g[0] == (0, 0.0, 0.0) and g[1] == (1, 0.0, 5.0) and g[5] == (5, 2.0, 5.0)
    with pytest.raises(TypeError):
        parameter_grid(None)  # type: ignore[arg-type]


def test_point_params_reduced_roundtrip() -> None:
    scan = ScanConfig(alpha_hat=(1.0,), gamma_hat=(1.0,))
    for n in (12, 100, 200):
        p = point_params(2.5, 3.0, n, scan)
        a, g = raw_to_reduced(p, n)
        assert a == pytest.approx(2.5, abs=1e-12) and g == pytest.approx(3.0, abs=1e-12)
        assert p.alpha == pytest.approx(2.0 * 2.5 / (n - 2), abs=1e-15)
    assert point_params(2.5, 3.0, 50, None) == point_params(2.5, 3.0, 50, scan)
    p = point_params(1.0, 2.0, 50, scan, beta=2.0)
    assert p.beta == 2.0 and raw_to_reduced(p, 50) == pytest.approx((1.0, 2.0), abs=1e-12)


def test_point_params_raw_scaling() -> None:
    raw = ScanConfig(alpha_hat=(0.3,), gamma_hat=(0.1,), scaling="raw")
    p = point_params(0.3, 0.1, 200, raw)
    assert (p.alpha, p.beta, p.gamma) == (0.3, 1.0, 0.1)
    p2 = point_params(0.3, 0.1, 200, raw, beta=2.0)
    assert (p2.alpha, p2.gamma) == (0.6, 0.2)


def test_simulate_reproducible_by_seedkey() -> None:
    cfg = small_cfg()
    p = point_params(3.0, 1.0, N)
    key = seed_key(cfg.seeds, 4, 1)
    r1, r2 = simulate(cfg, p, key), simulate(cfg, p, key)
    assert np.array_equal(r1.w0, r2.w0)
    assert np.array_equal(r1.trajectory.w_final, r2.trajectory.w_final)
    assert r1.trajectory.steps == r2.trajectory.steps
    assert run_summary(r1) == run_summary(r2)
    assert np.array_equal(r1.observables.geometry.d_eff.profile, r2.observables.geometry.d_eff.profile, equal_nan=True)
    r3 = simulate(cfg, p, seed_key(cfg.seeds, 4, 0))
    assert not np.array_equal(r1.w0, r3.w0)
    assert r1.seed == key and isinstance(r1.seed, SeedKey)


def test_simulate_phases_follow_bistability() -> None:
    cfg = small_cfg()
    key = seed_key(cfg.seeds, 0, 0)
    a = simulate(cfg, point_params(0.0, 0.0, N), key)
    e = simulate(cfg, point_params(3.0, 0.0, N), key)
    assert a.assessment.label is PhaseLabel.A and a.trajectory.status is RunStatus.CONVERGED
    assert e.assessment.label is PhaseLabel.E
    assert e.observables.topology.mean_weight == pytest.approx(1.0)
    assert a.alpha_hat == pytest.approx(0.0) and e.alpha_hat == pytest.approx(3.0)


def test_max_steps_gives_U() -> None:
    cfg = dataclasses.replace(small_cfg(), dynamics=DynamicsConfig(max_steps=2))
    r = simulate(cfg, point_params(3.0, 0.0, N), seed_key(cfg.seeds, 0, 0))
    assert r.trajectory.status is RunStatus.MAX_STEPS
    assert r.assessment.label is PhaseLabel.U


def test_simulate_from_checks_and_does_not_mutate() -> None:
    cfg = small_cfg()
    key = seed_key(cfg.seeds, 0, 0)
    w0 = simulate(cfg, point_params(0.0, 0.0, N), key).w0
    w0c = w0.copy()
    r = simulate_from(cfg, point_params(3.0, 0.0, N), key, w0)
    assert np.array_equal(w0, w0c) and np.array_equal(r.w0, w0c)
    with pytest.raises(ValueError):
        simulate_from(cfg, point_params(3.0, 0.0, N), key, np.zeros((N + 1, N + 1)))


def test_scan_phase_diagram_matches_direct_simulation_and_order_independent() -> None:
    scan = ScanConfig(alpha_hat=(0.0, 3.0), gamma_hat=(0.0, 1.0))
    cfg = small_cfg(scan, replicates=2)
    res = scan_phase_diagram(cfg)
    assert len(res) == 8
    # el punto 3 (alpha_hat=3, gamma_hat=1), replica 1 se reproduce sin ejecutar los demas
    direct = simulate(cfg, point_params(3.0, 1.0, N, scan), seed_key(cfg.seeds, 3, 1))
    got = res[3 * 2 + 1]
    assert np.array_equal(direct.trajectory.w_final, got.trajectory.w_final)
    assert got.seed.spawn_key == (3, 3, 1)
    assert [r.assessment.label for r in res[:4]] == [PhaseLabel.A] * 4
    assert [r.assessment.label for r in res[4:]] == [PhaseLabel.E] * 4
    with pytest.raises(ValueError):
        scan_phase_diagram(small_cfg(None))


def test_phase_probabilities_wilson() -> None:
    scan = ScanConfig(alpha_hat=(0.0, 3.0), gamma_hat=(0.0,))
    res = scan_phase_diagram(small_cfg(scan, replicates=3))
    pr = phase_probabilities(res)
    assert set(pr) == {(0.0, 0.0), (3.0, 0.0)}
    a, e = pr[(0.0, 0.0)], pr[(3.0, 0.0)]
    assert a["n"] == 3 and a["counts"]["A"] == 3 and a["probabilities"]["A"] == 1.0
    assert e["counts"]["E"] == 3 and e["probabilities"]["F"] == 0.0
    assert a["ci"]["A"] == proportion_ci(3, 3)
    assert set(a["counts"]) == {lab.value for lab in PhaseLabel}
    lo, hi = a["ci"]["A"]
    assert 0.0 < lo < 1.0 and hi == pytest.approx(1.0)


def test_static_run_and_observe() -> None:
    cfg = small_cfg()
    key = seed_key(cfg.seeds, 0, 0)
    w = simulate(cfg, point_params(0.0, 0.0, N), key).w0
    r = static_run(cfg, FunctionalParams(alpha=0.0), key, w)
    assert r.trajectory.steps == 0 and r.trajectory.dt == 0.0
    assert r.trajectory.snapshots.shape == (1, N * (N - 1) // 2)
    assert run_summary(r)["mean_strength"] == observe(w, cfg).topology.mean_strength
    assert r.assessment.label is PhaseLabel.E  # U(0,1): rho(W>0.1) ~ 0.9
    with pytest.raises(ValueError):
        static_run(cfg, FunctionalParams(alpha=0.0), key, np.zeros((N + 2, N + 2)))


def test_with_helpers() -> None:
    cfg = small_cfg()
    assert with_experiment(cfg, 9).seeds.experiment_id == 9 and cfg.seeds.experiment_id == 3
    p = FunctionalParams(alpha=0.1, gamma=0.2)
    c2 = with_params(cfg, p, n=20)
    assert c2.functional == p and c2.init.n == 20 and cfg.init.n == N


def test_to_jsonable_and_run_summary_serializable() -> None:
    assert to_jsonable({"a": float("nan"), 1: np.float64(2.0), "b": (np.int64(3),), "l": PhaseLabel.F}) == {
        "a": None, "1": 2.0, "b": [3], "l": "F"
    }
    cfg = small_cfg()
    r = simulate(cfg, point_params(0.0, 0.0, N), seed_key(cfg.seeds, 0, 0))
    s = run_summary(r)
    json.dumps(s, allow_nan=False)
    assert s["label"] == "A" and s["status"] == "CONVERGED"
