"""Tests de la taxonomia Omega-F0..F10 (WP-F, OMEGA_1_1_DESIGN §1.11)."""

from __future__ import annotations

import dataclasses
import math
from typing import Any

import pytest

from omega.certificate.taxonomy import (
    RUN_FLAG_KEYS,
    assess_run,
    consensus_dimension,
    primary_code,
    resistance_consistent,
    run_failure_codes,
    run_flags,
)
from omega.config.settings11 import DistanceSuiteConfig
from omega.contracts import FAILURE_PRECEDENCE, FailureCode, IsotropyReport, LocalityReport
from omega.types import RunStatus
from tests.factories_v11 import make_cfg11, make_dim, make_evidence, make_topology, with_dimensions

CFG = make_cfg11()
F = FailureCode


def test_all_pass_state() -> None:
    a = assess_run(make_evidence(), CFG)
    assert a.codes == () and a.primary is None and a.passes
    flags = run_flags(make_evidence(), CFG)
    assert tuple(flags) == RUN_FLAG_KEYS
    assert all(type(v) is bool for v in flags.values())
    assert flags["nontrivial"] and not flags["empty"] and not flags["dense_or_uniform"] and not flags["small_world"]


def _isolated() -> list[tuple[str, dict[str, Any], FailureCode]]:
    base = make_evidence()
    return [
        ("converged", dict(status=RunStatus.MAX_STEPS), F.F10),
        ("nonfinite", dict(status=RunStatus.NONFINITE), F.F10),
        ("dense_rho", dict(topology=make_topology(binary_density=0.6)), F.F1),
        ("dense_meanw", dict(weight_mean=0.6), F.F1),
        ("uniform", dict(weight_cv=0.0, weight_mean=0.3), F.F1),
        ("connected", dict(topology=make_topology(giant_fraction=0.9)), F.F2),
        ("locality", dict(locality=LocalityReport(detour_fraction=0.9, ok=False)), F.F3),
        ("small_world", dict(clustering_ratio=5.0, curvature=dataclasses.replace(base.curvature, tail_fraction=0.1)), F.F3),
        ("d_volume", dict(geometry=dataclasses.replace(base.geometry, d_eff=make_dim(3.5))), F.F4),
        ("d_spectral", dict(geometry=dataclasses.replace(base.geometry, d_s=make_dim(2.5))), F.F4),
        ("d_weyl", dict(d_weyl=make_dim(3.0, plateau=False)), F.F4),
        ("d_status", dict(d_weyl=make_dim(3.0, status="no_window")), F.F4),
        ("metric", dict(suite=dataclasses.replace(base.suite, distance_sensitive=True)), F.F5),
        ("zeta", dict(suite=dataclasses.replace(base.suite, resistance_exponent=0.6)), F.F5),
        ("isotropy", dict(isotropy=IsotropyReport(0.2, 0.1, 4, 3, 64, False)), F.F9),
        ("homogeneity", dict(homogeneity=dataclasses.replace(base.homogeneity, ok=False)), F.F9),
        ("topology", dict(topo=dataclasses.replace(base.topo, stable=False)), F.F9),
        ("manifold", dict(annulus=dataclasses.replace(base.annulus, ok=False)), F.F9),
        ("curvature", dict(curvature=dataclasses.replace(base.curvature, ok=False)), F.F9),
    ]


@pytest.mark.parametrize("name,over,code", _isolated(), ids=[x[0] for x in _isolated()])
def test_each_flag_isolated_produces_its_code(name: str, over: dict[str, Any], code: FailureCode) -> None:
    a = assess_run(make_evidence(**over), CFG)
    expected = (F.F3, F.F9) if name == "small_world" else (code,)  # cola de curvatura rompe tambien manifold_proxy
    assert a.codes == expected, (name, a.codes)
    assert a.primary is code and not a.passes


def test_empty_state_is_f0_with_full_code_set() -> None:
    ev = make_evidence(topology=make_topology(giant_fraction=0.05, mean_binary_degree=0.5, binary_density=0.0005))
    a = assess_run(ev, CFG)
    assert F.F0 in a.codes and F.F2 in a.codes  # se informan todos
    assert a.primary is F.F0
    assert a.flags["empty"] and not a.flags["nontrivial"]


def test_uniform_state_is_f1() -> None:
    ev = make_evidence(weight_cv=0.0, weight_mean=0.3, topology=make_topology(binary_density=1.0, mean_binary_degree=799.0))
    a = assess_run(ev, CFG)
    assert F.F1 in a.codes and a.flags["dense_or_uniform"] and not a.flags["nontrivial"]
    assert a.primary is F.F1


def test_not_converged_is_f10_and_dominates() -> None:
    ev = make_evidence(status=RunStatus.MAX_STEPS, topology=make_topology(giant_fraction=0.9))
    a = assess_run(ev, CFG)
    assert F.F10 in a.codes and F.F2 in a.codes
    assert a.primary is F.F10


def test_class_none_fails_metric_robust() -> None:
    ev = make_evidence(dimension_class=None)
    flags = run_flags(ev, CFG)
    assert flags["metric_robust"] is False
    assert assess_run(ev, CFG).primary is F.F5


def test_nonfinite_evidence_counts_as_false() -> None:
    nan = math.nan
    ev = make_evidence(consensus_dimension=nan)
    assert not run_flags(ev, CFG)["d_volume"]
    ev2 = make_evidence(clustering_ratio=nan)
    assert run_flags(ev2, CFG)["locality"] is False
    ev3 = make_evidence(topology=make_topology(giant_fraction=nan))
    f3 = run_flags(ev3, CFG)
    assert not f3["connected"] and not f3["nontrivial"] and not f3["empty"]
    assert F.F1 in assess_run(ev3, CFG).codes  # conservador: no puede pasar


def test_precedence_order_and_primary() -> None:
    assert FAILURE_PRECEDENCE == (F.F10, F.F0, F.F1, F.F2, F.F3, F.F4, F.F5, F.F7, F.F6, F.F8, F.F9)
    assert primary_code([]) is None
    remaining = list(FAILURE_PRECEDENCE)
    while remaining:
        assert primary_code(reversed(remaining)) is remaining[0]
        remaining.pop(0)
    assert primary_code([F.F6, F.F7]) is F.F7  # F7 > F6
    assert primary_code([F.F8, F.F9]) is F.F8
    assert primary_code([F.F9, F.F6]) is F.F6


def test_run_failure_codes_all_failed_ordered() -> None:
    flags = {k: False for k in RUN_FLAG_KEYS}
    flags["empty"] = True
    flags["dense_or_uniform"] = True
    codes = run_failure_codes(flags)
    assert codes == (F.F10, F.F0, F.F1, F.F2, F.F3, F.F4, F.F5, F.F9)
    assert run_failure_codes({}) != ()  # bandera ausente = no cumplida


def test_run_failure_codes_rejects_np_bool() -> None:
    import numpy as np

    flags: dict[str, bool] = {k: True for k in RUN_FLAG_KEYS}
    flags["empty"] = False
    flags["dense_or_uniform"] = False
    flags["small_world"] = False
    assert run_failure_codes(flags) == ()
    flags["connected"] = np.True_  # type: ignore[assignment]
    with pytest.raises(TypeError):
        run_failure_codes(flags)


def test_consensus_dimension() -> None:
    est = [make_dim(3.06), make_dim(2.60), make_dim(2.91)]
    d, cls = consensus_dimension(est, 0.35)
    assert d == pytest.approx(2.91) and cls == 3
    d, cls = consensus_dimension([make_dim(3.17), make_dim(1.85), make_dim(1.55)], 0.35)
    assert d == pytest.approx(1.85) and cls == 2
    d, cls = consensus_dimension([make_dim(2.5)] * 3, 0.35)
    assert d == pytest.approx(2.5) and cls is None
    d, cls = consensus_dimension([make_dim(3.0), make_dim(3.0), make_dim(3.0, status="no_window")], 0.35)
    assert math.isnan(d) and cls is None
    d, cls = consensus_dimension([make_dim(3.0), make_dim(math.nan), make_dim(3.0)], 0.35)
    assert math.isnan(d) and cls is None
    d, cls = consensus_dimension([make_dim(0.1)] * 3, 0.35)
    assert cls is None  # clase 0 no existe
    assert consensus_dimension([], 0.35)[1] is None


def test_resistance_consistent() -> None:
    s = DistanceSuiteConfig()
    assert resistance_consistent(0.9, 1, s) and not resistance_consistent(0.7, 1, s)
    assert resistance_consistent(0.5, 2, s) and not resistance_consistent(0.7, 2, s) and not resistance_consistent(0.2, 2, s)
    assert resistance_consistent(0.2, 3, s) and resistance_consistent(0.0, 4, s) and not resistance_consistent(0.35, 3, s)
    assert not resistance_consistent(0.2, None, s)
    assert not resistance_consistent(math.nan, 3, s)


def test_with_dimensions_helper_consistent() -> None:
    ev = with_dimensions(make_evidence(), 2.0)
    assert assess_run(ev, dataclasses.replace(CFG)).codes[:1] in ((), (F.F5,))  # zeta 0.2 no coherente con clase 2
