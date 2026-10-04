"""Tests del GeometryCertificate y del veredicto por punto (WP-F, OMEGA_1_1_DESIGN §1.12)."""

from __future__ import annotations

import dataclasses
import itertools
import math
from typing import Any

import numpy as np
import pytest

from omega.certificate.certificate import (
    RULE_D3_NEVER_SUFFICIENT,
    build_certificate,
    finalize_verdict,
    null_separated,
    point_verdict,
    seed_robust,
    size_robust,
)
from omega.config.settings11 import NullModel
from omega.contracts import (
    CERTIFICATE_FIELDS,
    EnsembleSummary,
    FailureCode,
    GeometryCertificate,
    LocalityReport,
    RunEvidence,
    Verdict,
)
from omega.types import RunStatus
from tests.factories_v11 import make_cfg11, make_evidence, make_series, with_dimensions

CFG = make_cfg11()
F = FailureCode


def _cert(values: tuple[bool, ...], cls: int | None) -> GeometryCertificate:
    return GeometryCertificate(
        **dict(zip(CERTIFICATE_FIELDS, values)), consensus_dimension=3.0, dimension_class=cls, n_runs=5, evidence={}
    )


def test_rule_constant_is_text() -> None:
    assert isinstance(RULE_D3_NEVER_SUFFICIENT, str) and "jamás" in RULE_D3_NEVER_SUFFICIENT


def test_exhaustive_2_pow_14_candidate_only_if_all_true() -> None:
    assert len(CERTIFICATE_FIELDS) == 14
    n_candidates = 0
    for values in itertools.product((False, True), repeat=14):
        for cls in (None, 3):
            cert = _cert(values, cls)
            for codes in ((), (F.F7,), (F.F10,)):
                v = finalize_verdict(cert, codes)
                expected = all(values) and cls is not None and codes == ()
                assert (v is Verdict.GEOMETRIC_CANDIDATE) == expected
                n_candidates += expected
            assert (cert.verdict is Verdict.GEOMETRIC_CANDIDATE) == (all(values) and cls is not None)
    assert n_candidates == 1


def test_np_bool_rejected() -> None:
    kw: dict[str, Any] = dict(zip(CERTIFICATE_FIELDS, (True,) * 14))
    kw["connected"] = np.True_
    with pytest.raises(TypeError):
        GeometryCertificate(**kw, consensus_dimension=3.0, dimension_class=3, n_runs=1, evidence={})


# ------------------------------------------------------------------ escenario completo


def _nulls(d_null: float = 1.0) -> dict[NullModel, list[RunEvidence]]:
    # Replicas nulas que fallan (no convergen) y con D* distinto.
    bad = with_dimensions(make_evidence(status=RunStatus.MAX_STEPS), d_null)
    return {m: [bad] * 4 for m in CFG.null_models}


def _by_size(series: tuple[float, ...] = (3.06, 3.11, 3.16)) -> dict[int, list[RunEvidence]]:
    return {n: make_series(d, 4) for n, d in zip((100, 200, 400), series)}


def _good_point() -> tuple[list[RunEvidence], dict[int, list[RunEvidence]], dict[NullModel, list[RunEvidence]]]:
    return make_series(3.0, 5), _by_size(), _nulls()


def test_full_pass_is_candidate() -> None:
    runs, by_size, nulls = _good_point()
    pv = point_verdict(runs, by_size, nulls, CFG)
    assert pv.certificate.satisfied_all, pv.certificate.failed_fields
    assert pv.certificate.dimension_class == 3
    assert pv.codes == () and pv.primary is None
    assert pv.verdict is Verdict.GEOMETRIC_CANDIDATE
    assert pv.run_outcome_fractions == {"PASS": 1.0}
    assert "missing" not in pv.certificate.evidence


@pytest.mark.parametrize(
    "field,over",
    [
        ("connected", dict(topology=None)),
        ("locality", dict(locality=LocalityReport(0.5, False))),
        ("metric_robust", dict(dimension_class=None)),
    ],
)
def test_run_level_field_failures_block_candidate(field: str, over: dict[str, Any]) -> None:
    runs, by_size, nulls = _good_point()
    from tests.factories_v11 import make_topology

    if over.get("topology", 0) is None:
        over = dict(topology=make_topology(giant_fraction=0.5))
    bad = dataclasses.replace(runs[0], **over)
    pv = point_verdict([bad] * 5, by_size, nulls, CFG)
    assert getattr(pv.certificate, field) is False
    assert pv.verdict is Verdict.NOT_CANDIDATE
    assert pv.primary is not None


def test_eighty_percent_rule_at_point_level() -> None:
    runs, by_size, nulls = _good_point()
    bad = dataclasses.replace(runs[0], locality=LocalityReport(0.5, False))
    ok4 = point_verdict(runs[:4] + [bad], by_size, nulls, CFG)  # 4/5 = 80% -> campo True
    assert ok4.certificate.locality is True
    ok3 = point_verdict(runs[:3] + [bad, bad], by_size, nulls, CFG)  # 60% -> False
    assert ok3.certificate.locality is False
    assert ok3.verdict is Verdict.NOT_CANDIDATE


def test_exact_d3_with_locality_false_is_not_candidate_primary_f3() -> None:
    runs = make_series(3.0, 5, locality=LocalityReport(0.5, False))
    assert all(r.consensus_dimension == 3.0 and r.d_weyl.value == 3.0 for r in runs)
    pv = point_verdict(runs, _by_size(), _nulls(), CFG)
    assert pv.certificate.dimension_class == 3
    assert pv.certificate.d_volume and pv.certificate.d_spectral and pv.certificate.d_weyl
    assert not pv.certificate.locality
    assert pv.verdict is Verdict.NOT_CANDIDATE
    assert pv.primary is F.F3


def test_missing_evidence_is_false_not_candidate() -> None:
    runs = make_series(3.0, 5)
    pv = point_verdict(runs, {}, {}, CFG)
    c = pv.certificate
    assert not c.size_robust and not c.null_separated
    assert pv.verdict is Verdict.NOT_CANDIDATE
    assert F.F6 in pv.codes and F.F8 in pv.codes
    assert "sizes" in str(c.evidence["missing"]) and "nulls" in str(c.evidence["missing"])
    # sin corridas
    pv0 = point_verdict([], _by_size(), _nulls(), CFG)
    assert pv0.verdict is Verdict.NOT_CANDIDATE and F.F10 in pv0.codes
    assert pv0.certificate.n_runs == 0 and not pv0.certificate.satisfied_all
    # menos de min_sizes tamanos
    pv2 = point_verdict(runs, dict(list(_by_size().items())[:2]), _nulls(), CFG)
    assert not pv2.certificate.size_robust and pv2.verdict is Verdict.NOT_CANDIDATE


def test_missing_null_model_rejects() -> None:
    runs = make_series(3.0, 5)
    nulls = _nulls()
    nulls.pop(NullModel.ERDOS_RENYI)
    ok, ev = null_separated(runs, nulls, CFG)
    assert not ok and ev["null_erdos_renyi_missing"] == 1.0


# ------------------------------------------------------------------ size_robust


def test_size_robust_accepts_rgg_series_and_rejects_trend() -> None:
    ok, ev = size_robust(_by_size((3.06, 3.11, 3.16)), CFG)
    assert ok, ev
    assert ev["size_range"] == pytest.approx(0.10)
    assert abs(ev["size_slope"]) <= 0.1
    bad, _ = size_robust(_by_size((2.6, 3.0, 3.4)), CFG)
    assert not bad
    bad2, ev2 = size_robust(_by_size((2.75, 3.0, 3.25)), CFG)  # corridas limpias pero rango 0.5 > 0.3
    assert not bad2 and ev2["size_range"] > 0.3


def test_size_robust_slope_and_cv_rules() -> None:
    # Pendiente grande con rango pequeno: tamanos muy juntos en ln N.
    by_size = {100: make_series(2.95, 3), 110: make_series(3.05, 3), 120: make_series(3.15, 3)}
    ok, ev = size_robust(by_size, CFG)
    assert not ok and abs(ev["size_slope"]) > 0.1 and ev["size_range"] <= 0.3
    # cv de homogeneidad crece con N
    by_cv = {100: make_series(3.06, 3, cv=0.05), 200: make_series(3.11, 3, cv=0.06), 400: make_series(3.16, 3, cv=0.10)}
    ok, ev = size_robust(by_cv, CFG)
    assert not ok and ev["size_cv_large"] > ev["size_cv_small"] + 0.02
    by_cv_ok = {100: make_series(3.06, 3, cv=0.05), 200: make_series(3.11, 3, cv=0.06), 400: make_series(3.16, 3, cv=0.07)}
    assert size_robust(by_cv_ok, CFG)[0]


def test_size_robust_requires_all_sizes_clean() -> None:
    by_size = _by_size()
    by_size[400] = [dataclasses.replace(r, status=RunStatus.MAX_STEPS) for r in by_size[400]]
    ok, ev = size_robust(by_size, CFG)
    assert not ok and ev["size_400_frac_clean"] == 0.0


# ------------------------------------------------------------------ seed_robust


def test_seed_robust() -> None:
    ok, ev = seed_robust(make_series(3.0, 5), CFG)
    assert ok and ev["seed_std_dstar"] == 0.0
    assert not seed_robust(make_series(3.0, 1), CFG)[0]  # una semilla no basta
    assert not seed_robust([], CFG)[0]
    wide = make_series(2.8, 3) + make_series(3.3, 3)  # std 0.25 > 0.2
    assert not seed_robust(wide, CFG)[0]
    split = make_series(3.0, 4) + make_series(2.0, 4)  # clases distintas
    assert not seed_robust(split, CFG)[0]


def test_seed_dependent_point_is_f7() -> None:
    runs = make_series(3.0, 3) + [dataclasses.replace(make_evidence(), locality=LocalityReport(0.5, False))] * 2
    pv = point_verdict(runs, _by_size(), _nulls(), CFG)
    assert F.F7 in pv.codes and pv.verdict is Verdict.NOT_CANDIDATE
    assert pv.run_outcome_fractions["PASS"] == pytest.approx(0.6)
    assert pv.primary is F.F3 or pv.primary is F.F7


def test_no_modal_outcome_reaches_80_gives_f7_primary() -> None:
    pass_run = make_evidence()
    f2 = make_evidence(topology=dataclasses.replace(make_evidence().topology, giant_fraction=0.9, giant_size=720))
    runs = [pass_run] * 3 + [f2] * 2
    pv = point_verdict(runs, _by_size(), _nulls(), CFG)
    assert F.F7 in pv.codes and F.F2 not in pv.codes
    assert pv.primary is F.F7


def test_modal_failure_code_is_reported() -> None:
    f2 = make_evidence(topology=dataclasses.replace(make_evidence().topology, giant_fraction=0.9, giant_size=720))
    pv = point_verdict([f2] * 5, _by_size(), _nulls(), CFG)
    assert F.F2 in pv.codes and pv.primary is F.F2


# ------------------------------------------------------------------ F10 / ensemble


def test_f10_from_convergence_fraction_and_ensemble() -> None:
    runs = make_series(3.0, 3) + [make_evidence(status=RunStatus.MAX_STEPS)] * 2  # 60% converged
    pv = point_verdict(runs, _by_size(), _nulls(), CFG)
    assert pv.primary is F.F10 and pv.verdict is Verdict.NOT_CANDIDATE
    runs, by_size, nulls = _good_point()
    ens = EnsembleSummary(means={}, ses={}, rhat={}, tau_int={}, ess={}, geweke_ok=True, equilibrated=False)
    pv = point_verdict(runs, by_size, nulls, CFG, ensemble=ens)
    assert F.F10 in pv.codes and pv.verdict is Verdict.NOT_CANDIDATE
    ens_ok = dataclasses.replace(ens, equilibrated=True)
    assert point_verdict(runs, by_size, nulls, CFG, ensemble=ens_ok).verdict is Verdict.GEOMETRIC_CANDIDATE


# ------------------------------------------------------------------ null_separated


def test_null_separated_rejects_passing_null() -> None:
    runs = make_series(3.0, 5)
    ok, _ = null_separated(runs, _nulls(), CFG)
    assert ok
    nulls = _nulls()
    nulls[NullModel.ERDOS_RENYI] = make_series(3.0, 4)  # nulo que pasa todo
    ok, ev = null_separated(runs, nulls, CFG)
    assert not ok and ev["null_erdos_renyi_frac_pass"] == 1.0


def test_null_separated_dimension_gap_rules() -> None:
    runs = make_series(3.0, 5)
    near = with_dimensions(make_evidence(status=RunStatus.MAX_STEPS), 2.8)  # falla pero D* cercano (sigma=0, gap .2 < .35)
    nulls = _nulls()
    nulls[NullModel.SHUFFLED_WEIGHTS] = [near] * 4
    ok, ev = null_separated(runs, nulls, CFG)
    assert not ok and ev["null_shuffled_weights_gap"] == pytest.approx(0.2)
    # sigma > 0: se exige gap >= 3 sigma
    spread = [with_dimensions(make_evidence(status=RunStatus.MAX_STEPS), d) for d in (1.5, 2.0, 2.5, 3.0)]
    nulls[NullModel.SHUFFLED_WEIGHTS] = spread
    sigma = float(np.std([1.5, 2.0, 2.5, 3.0], ddof=1))
    mu = 2.25
    assert (3.0 - mu) < 3 * sigma
    assert not null_separated(runs, nulls, CFG)[0]
    far = [with_dimensions(make_evidence(status=RunStatus.MAX_STEPS), d) for d in (1.0, 1.05, 0.95, 1.0)]
    nulls[NullModel.SHUFFLED_WEIGHTS] = far
    assert null_separated(runs, nulls, CFG)[0]


def test_null_with_nan_dimension_only_needs_failing() -> None:
    runs = make_series(3.0, 5)
    nan_null = make_evidence(status=RunStatus.MAX_STEPS, consensus_dimension=math.nan, dimension_class=None)
    nulls = {m: [nan_null] * 3 for m in CFG.null_models}
    assert null_separated(runs, nulls, CFG)[0]


def test_candidate_with_nonfinite_dimension_never_separated() -> None:
    bad = make_evidence(consensus_dimension=math.nan, dimension_class=None)
    assert not null_separated([bad] * 3, _nulls(), CFG)[0]
    cert = build_certificate([bad] * 3, _by_size(), _nulls(), CFG)
    assert cert.dimension_class is None and math.isnan(cert.consensus_dimension)
    assert cert.verdict is Verdict.NOT_CANDIDATE


def test_certificate_fields_are_strict_bool() -> None:
    runs, by_size, nulls = _good_point()
    cert = build_certificate(runs, by_size, nulls, CFG)
    for name in CERTIFICATE_FIELDS:
        assert type(getattr(cert, name)) is bool
    assert cert.n_runs == 5
