"""Clasificacion de fases con Observables sinteticos (ANALYSIS §4.1, §6.2)."""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import pytest

from omega.config.settings import PhaseThresholds
from omega.phases.classification import (
    classify_run,
    is_clustered,
    is_connected_homogeneous,
    is_dispersed,
    is_fragmented,
    is_geometric_candidate,
    is_hyperdense,
    phase_flags,
)
from omega.types import (
    DimensionEstimate,
    GeometryObservables,
    Observables,
    PhaseLabel,
    RunStatus,
    TopologyObservables,
)

T = PhaseThresholds()
_Z = np.zeros(0, dtype=np.float64)


def dim(value: float = 3.0, status: str = "ok", plateau: bool = True) -> DimensionEstimate:
    return DimensionEstimate(value, 0.01, status, (0, 2), _Z, _Z, _Z, plateau, "shell")  # type: ignore[arg-type]


def obs(topo: dict[str, Any] | None = None, d_eff: DimensionEstimate | None = None,
        d_s: DimensionEstimate | None = None) -> Observables:
    """Observables de una fase C 'limpia' (G=1, cv bajo, bajo clustering, rho moderada)."""
    base: dict[str, Any] = dict(
        n=200, mean_strength=5.0, std_strength=0.2, cv_strength=0.04, mean_binary_degree=6.0,
        binary_density=0.03, mean_weight=0.03, giant_size=200, giant_fraction=1.0, n_components=1,
        n_large_components=1, clustering_weighted=0.05, clustering_binary=0.05,
        frac_at_zero=0.9, frac_at_one=0.01,
    )
    base.update(topo or {})
    nan = float("nan")
    geo = GeometryObservables(
        path_length=5.0, path_length_hops=5.0, diameter=9.0,
        d_eff=d_eff or dim(nan, "no_window", False), d_eff_ball=dim(nan, "no_window", False),
        d_eff_hops=dim(nan, "no_window", False), d_s=d_s or dim(nan, "no_window", False),
    )
    return Observables(TopologyObservables(**base), geo)


def label(o: Observables, status: RunStatus = RunStatus.CONVERGED) -> PhaseLabel:
    return classify_run(o, status, T).label


def test_phase_E_by_density_or_mean_weight() -> None:
    assert label(obs({"binary_density": 0.5})) is PhaseLabel.E
    assert label(obs({"mean_weight": 0.5})) is PhaseLabel.E
    assert is_hyperdense(obs({"binary_density": 0.49, "mean_weight": 0.49}), T) is False


def test_phase_A() -> None:
    o = obs({"giant_fraction": 0.005, "mean_binary_degree": 0.0, "binary_density": 0.0,
             "mean_weight": 0.0, "n_large_components": 0})
    assert is_dispersed(o, T) and label(o) is PhaseLabel.A
    assert not is_dispersed(obs({"giant_fraction": 0.05, "mean_binary_degree": 1.0}), T)


def test_phase_B() -> None:
    o = obs({"giant_fraction": 0.5, "n_large_components": 2})
    assert is_fragmented(o, T) and label(o) is PhaseLabel.B
    assert not is_fragmented(obs({"giant_fraction": 0.5, "n_large_components": 1}), T)
    assert not is_fragmented(obs({"giant_fraction": 0.95, "n_large_components": 2}), T)


def test_phase_F_candidate_and_its_conditions() -> None:
    o = obs(d_eff=dim(3.0), d_s=dim(3.2))
    assert is_geometric_candidate(o, T) and label(o) is PhaseLabel.F
    assert not is_geometric_candidate(obs(d_eff=dim(3.0, plateau=False), d_s=dim(3.0)), T)
    assert not is_geometric_candidate(obs(d_eff=dim(3.0), d_s=dim(3.6)), T)  # |dD| > 0.5
    assert not is_geometric_candidate(obs(d_eff=dim(3.0, "no_window", True), d_s=dim(3.0)), T)
    assert not is_geometric_candidate(obs(d_eff=dim(3.0), d_s=dim(3.0, "no_window")), T)
    assert not is_geometric_candidate(obs({"giant_fraction": 0.93}, dim(3.0), dim(3.0)), T)
    assert not is_geometric_candidate(obs({"cv_strength": 0.25}, dim(3.0), dim(3.0)), T)
    assert not is_geometric_candidate(obs({"binary_density": 0.6}, dim(3.0), dim(3.0)), T)


def test_phase_D() -> None:
    o = obs({"clustering_binary": 0.5, "binary_density": 0.05})
    assert is_clustered(o, T) and label(o) is PhaseLabel.D
    assert not is_clustered(obs({"clustering_binary": 0.5, "binary_density": 0.3}), T)  # ratio < 3
    assert not is_clustered(obs({"clustering_binary": 0.2, "binary_density": 0.01}), T)  # C < 0.3
    assert not is_clustered(obs({"clustering_binary": 0.5, "binary_density": 0.0}), T)


def test_phase_C() -> None:
    o = obs()
    assert is_connected_homogeneous(o, T) and label(o) is PhaseLabel.C
    assert not is_connected_homogeneous(obs({"cv_strength": 0.3}), T)
    assert not is_connected_homogeneous(obs({"giant_fraction": 0.85}), T)


def test_phase_U_when_nothing_matches() -> None:
    assert label(obs({"cv_strength": 0.5})) is PhaseLabel.U  # conexa pero heterogenea
    assert label(obs({"giant_fraction": 0.5, "n_large_components": 1, "mean_binary_degree": 3.0})) is PhaseLabel.U


def test_precedence_E_over_everything() -> None:
    o = obs({"binary_density": 0.9, "mean_weight": 0.9, "clustering_binary": 1.0}, dim(3.0), dim(3.0))
    f = phase_flags(o, T)
    assert f["E"] and f["C"]
    assert label(o) is PhaseLabel.E


def test_precedence_A_over_B_and_F_over_D_over_C() -> None:
    # F > D > C: con banderas F, D y C simultaneas gana F
    o = obs({"clustering_binary": 0.5, "binary_density": 0.05}, dim(3.0), dim(3.0))
    f = phase_flags(o, T)
    assert f["F_candidate"] and f["D"] and f["C"]
    assert label(o) is PhaseLabel.F
    # sin F: D gana a C
    o2 = obs({"clustering_binary": 0.5, "binary_density": 0.05})
    assert phase_flags(o2, T)["C"] and label(o2) is PhaseLabel.D
    # B gana a F/D/C cuando aplica y no hay E/A
    o3 = obs({"giant_fraction": 0.5, "n_large_components": 2})
    assert label(o3) is PhaseLabel.B
    # E gana a A (A exige G<0.1; forzamos ambas con densidad alta)
    o4 = obs({"giant_fraction": 0.05, "mean_binary_degree": 0.5, "binary_density": 0.7})
    f4 = phase_flags(o4, T)
    assert f4["E"] and f4["A"] and label(o4) is PhaseLabel.E


@pytest.mark.parametrize("status", [RunStatus.MAX_STEPS, RunStatus.NONFINITE])
def test_nonconverged_gives_U_but_keeps_flags(status: RunStatus) -> None:
    o = obs({"binary_density": 0.9, "mean_weight": 0.9})
    a = classify_run(o, status, T)
    assert a.label is PhaseLabel.U
    assert a.flags["E"] is True and a.flags["converged"] is False
    assert set(a.flags) == {"E", "A", "B", "F_candidate", "D", "C", "converged"}


def test_converged_flag_true_and_thresholds_are_used() -> None:
    a = classify_run(obs(), RunStatus.CONVERGED, T)
    assert a.flags["converged"] is True
    strict = dataclasses.replace(T, cv_homogeneous=0.01)
    assert classify_run(obs(), RunStatus.CONVERGED, strict).label is PhaseLabel.U


def test_classify_run_rejects_bad_status() -> None:
    with pytest.raises(TypeError):
        classify_run(obs(), "CONVERGED", T)  # type: ignore[arg-type]
