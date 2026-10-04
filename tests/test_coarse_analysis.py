"""omega.coarse_graining.analysis: jerarquia de evidencia y consistencia de clase (anillo y 2D conservan clase)."""

from __future__ import annotations

import dataclasses

import pytest

from omega.certificate.evidence import collect_run_evidence
from omega.coarse_graining.analysis import class_consistency, hierarchy_evidence
from omega.config.seeds import SeedKey, make_rng
from omega.config.settings11 import CoarseGrainConfig, Omega11Config
from omega.contracts import AnnulusReport, FailureCode, LocalityReport
from omega.experiments.reference_graphs import periodic_lattice
from omega.phases.scan import default_config
from omega.types import RunStatus
from tests.factories_v11 import make_evidence, with_dimensions

KEY = SeedKey(20240901, (12, 0, 0))


def _cfg(replicates: int = 1) -> Omega11Config:
    cfg = Omega11Config(base=default_config(200, master_entropy=20240901, replicates=3))
    return dataclasses.replace(cfg, coarse=CoarseGrainConfig(replicates=replicates))


def test_ring_conserves_class() -> None:
    cfg = _cfg(2)
    w = periodic_lattice((200,))
    levels = hierarchy_evidence(w, cfg, make_rng(KEY))
    assert levels[0].n == 200 and len(levels) >= 5  # nivel 0 + 2 replicas x (>= 2 niveles)
    cc = class_consistency(levels, cfg)
    measured = [c for c in cc.classes if c is not None]
    assert measured and set(measured) == {1}
    assert cc.measured_levels >= 3 and cc.consistent
    assert cc.measured_levels + cc.undetermined_levels == len(levels)


def test_torus_2d_conserves_class() -> None:
    cfg = _cfg(1)
    levels = hierarchy_evidence(periodic_lattice((20, 20)), cfg, make_rng(KEY))
    cc = class_consistency(levels, cfg)
    assert {c for c in cc.classes if c is not None} == {2}
    assert cc.measured_levels >= 2 and cc.consistent


def test_level0_is_plain_evidence_and_deterministic() -> None:
    cfg = _cfg(1)
    w = periodic_lattice((120,))
    a = hierarchy_evidence(w, cfg, make_rng(KEY))
    b = hierarchy_evidence(w, cfg, make_rng(KEY))
    assert [e.n for e in a] == [e.n for e in b] and all(x.n < 120 for x in a[1:])
    ref = collect_run_evidence(w, RunStatus.CONVERGED, cfg, make_rng(KEY))
    assert a[0].dimension_class == ref.dimension_class == 1
    assert a[0].consensus_dimension == ref.consensus_dimension


def test_hierarchy_stops_below_n_min() -> None:
    cfg = _cfg(1)
    levels = hierarchy_evidence(periodic_lattice((60,)), cfg, make_rng(KEY))
    assert len(levels) == 1  # 60 -> ~35 < n_min=50: sin niveles medibles


def test_class_change_is_inconsistent() -> None:
    cfg = _cfg()
    lv = [with_dimensions(make_evidence(), 3.0), with_dimensions(make_evidence(), 3.0), with_dimensions(make_evidence(), 2.0)]
    cc = class_consistency(lv, cfg)
    assert cc.classes == (3, 3, 2) and cc.measured_levels == 3 and not cc.consistent


def test_new_failure_code_is_inconsistent() -> None:
    cfg = _cfg()
    bad = make_evidence(locality=LocalityReport(detour_fraction=0.5, ok=False))
    cc = class_consistency([make_evidence(), bad], cfg)
    assert cc.codes[0] is None and cc.codes[1] is FailureCode.F3 and not cc.consistent


def test_preexisting_code_is_not_new() -> None:
    cfg = _cfg()
    ann = AnnulusReport(fractions={2: 0.2}, ok=False)
    cc = class_consistency([make_evidence(annulus=ann), make_evidence(annulus=ann)], cfg)
    assert cc.codes == (None, None) and cc.consistent


def test_undetermined_levels_do_not_count_as_failure() -> None:
    cfg = _cfg()
    und = with_dimensions(make_evidence(), 2.5, cls=None)
    ok = with_dimensions(make_evidence(), 3.0)
    cc = class_consistency([ok, und, ok], cfg)
    assert cc.classes == (3, None, 3) and cc.undetermined_levels == 1 and cc.measured_levels == 2 and cc.consistent
    only = class_consistency([ok, und], cfg)  # un solo nivel medible: nada que comparar
    assert not only.consistent and only.measured_levels == 1


def test_errors() -> None:
    cfg = _cfg()
    with pytest.raises(ValueError):
        class_consistency([], cfg)
    with pytest.raises(TypeError):
        hierarchy_evidence(periodic_lattice((60,)), cfg, "rng")  # type: ignore[arg-type]
