"""Pruebas minimas de WP-0: settings11, contratos y regla RULE_D3_NEVER_SUFFICIENT."""

from __future__ import annotations

import dataclasses
import itertools
from typing import Any

import numpy as np
import pytest

import omega
from omega.baseline import BASELINE_COMMIT, BASELINE_NAME, BASELINE_SCHEMA_VERSION, MODEL_VERSION
from omega.config.settings import FunctionalParams, OmegaConfig
from omega.config.settings11 import (
    GEODESIC_MODES,
    CertificateThresholds,
    CoarseGrainConfig,
    DistanceMode,
    DistanceSuiteConfig,
    Engine,
    EnsembleConfig,
    LangevinConfig,
    MetropolisConfig,
    NullModel,
    Omega11Config,
    TopologyConfig,
    WeylConfig,
)
from omega.contracts import (
    CERTIFICATE_FIELDS,
    FAILURE_MEANING,
    FAILURE_PRECEDENCE,
    FIELD_CODE,
    RULE_D3_NEVER_SUFFICIENT,
    FailureCode,
    GeometryCertificate,
    PointVerdict,
    Verdict,
)
from omega.phases.scan import default_config


def _cert(**over: Any) -> GeometryCertificate:
    kw: dict[str, Any] = {name: True for name in CERTIFICATE_FIELDS}
    kw.update(consensus_dimension=2.9, dimension_class=3, n_runs=10, evidence={"k": 1.0})
    kw.update(over)
    return GeometryCertificate(**kw)


def test_versions_and_baseline_constants() -> None:
    assert omega.__version__ == "1.1.0"
    assert BASELINE_NAME == "OMEGA_1_0_BASELINE"
    assert BASELINE_COMMIT == "9dfbea7d6be28a8f17fbf336579c23628ac5669d"
    assert BASELINE_SCHEMA_VERSION == "1.0"
    assert MODEL_VERSION == "Ω-1.1"


def test_defaults_construct_and_are_frozen() -> None:
    cfg = Omega11Config(base=default_config(40))
    assert cfg.schema_version == "1.1" and cfg.engine is Engine.GRADIENT
    assert len(cfg.null_models) == 4 and NullModel.SMALL_WORLD not in cfg.null_models
    assert GEODESIC_MODES == (DistanceMode.HOP, DistanceMode.WEIGHTED_INVERSE, DistanceMode.WEIGHTED_LOG)
    assert DistanceSuiteConfig().modes[-1] is DistanceMode.RESISTANCE
    with pytest.raises(dataclasses.FrozenInstanceError):
        cfg.engine = Engine.LANGEVIN  # type: ignore[misc]
    assert not hasattr(cfg, "__dict__")


@pytest.mark.parametrize(
    "factory",
    [
        lambda: DistanceSuiteConfig(metric_tol=0.0),
        lambda: DistanceSuiteConfig(modes=()),
        lambda: DistanceSuiteConfig(zeta_2d_min=0.7),
        lambda: WeylConfig(laplacian="x"),  # type: ignore[arg-type]
        lambda: WeylConfig(count_max_frac=0.0),
        lambda: TopologyConfig(thetas=(0.1, 0.5)),
        lambda: TopologyConfig(thetas=(1.5,)),
        lambda: CertificateThresholds(g_connected=1.5),
        lambda: CoarseGrainConfig(n_min=1),
        lambda: LangevinConfig(theta_hat=0.0),
        lambda: MetropolisConfig(theta_hat=1.0, acceptance_low=0.6, acceptance_high=0.5),
        lambda: EnsembleConfig(n_chains=1),
        lambda: EnsembleConfig(rhat_max=float("nan")),
    ],
)
def test_settings11_validation_rejects(factory: Any) -> None:
    with pytest.raises((ValueError, TypeError)):
        factory()


def test_omega11_engine_requirements() -> None:
    base = default_config(40)
    with pytest.raises(ValueError):
        Omega11Config(base=base, engine=Engine.LANGEVIN)
    with pytest.raises(ValueError):
        Omega11Config(base=base, engine=Engine.METROPOLIS)
    with pytest.raises(ValueError):
        Omega11Config(base=base, engine=Engine.FIXED_DENSITY)
    ok = Omega11Config(base=base, engine=Engine.LANGEVIN, langevin=LangevinConfig(theta_hat=1.0))
    assert ok.langevin is not None
    with pytest.raises(ValueError, match="eta"):  # Langevin exige eta == 0
        Omega11Config(
            base=dataclasses.replace(
                base,
                functional=FunctionalParams(alpha=0.0, eta=0.1),
                dynamics=dataclasses.replace(base.dynamics, dt_mode="fixed", dt_fixed=0.01),
            ),
            engine=Engine.LANGEVIN,
            langevin=LangevinConfig(theta_hat=1.0),
        )
    with pytest.raises(TypeError):
        Omega11Config(base="x")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        Omega11Config(base=base, schema_version="1.0")
    with pytest.raises(ValueError):
        Omega11Config(base=dataclasses.replace(base, schema_version="9.9"))
    assert isinstance(base, OmegaConfig)
    with pytest.raises(ValueError):
        Omega11Config(base=base, null_models=())


def test_certificate_verdict_requires_all_14_and_dimension_class() -> None:
    full = _cert()
    assert full.satisfied_all and full.failed_fields == ()
    assert full.verdict is Verdict.GEOMETRIC_CANDIDATE
    assert _cert(dimension_class=None).verdict is Verdict.NOT_CANDIDATE
    for name in CERTIFICATE_FIELDS:
        c = _cert(**{name: False})
        assert not c.satisfied_all and c.failed_fields == (name,)
        assert c.verdict is Verdict.NOT_CANDIDATE


def test_rule_d3_never_sufficient_exhaustive() -> None:
    """2^14 combinaciones: el veredicto solo es candidato con los 14 True y clase no None; D~3 no basta."""
    assert "jamás será suficiente por sí solo" in RULE_D3_NEVER_SUFFICIENT
    n_candidate = 0
    for bits in itertools.product((True, False), repeat=14):
        kw = dict(zip(CERTIFICATE_FIELDS, bits))
        for cls in (3, None):
            c = _cert(dimension_class=cls, **kw)
            expect = all(bits) and cls is not None
            assert (c.verdict is Verdict.GEOMETRIC_CANDIDATE) == expect
            n_candidate += expect
    assert n_candidate == 1


def test_certificate_rejects_numpy_bool_and_non_bool() -> None:
    for bad in (np.bool_(True), 1, 0, "True", None):
        with pytest.raises(TypeError):
            _cert(connected=bad)
    for name in CERTIFICATE_FIELDS:
        with pytest.raises(TypeError):
            _cert(**{name: np.True_})
    with pytest.raises(TypeError):
        _cert(dimension_class=True)
    with pytest.raises(ValueError):
        _cert(n_runs=-1)
    with pytest.raises(TypeError):
        _cert(evidence={"k": [1]})


def test_certificate_nan_dimension_allowed() -> None:
    c = _cert(consensus_dimension=float("nan"), dimension_class=None)
    assert c.verdict is Verdict.NOT_CANDIDATE


def test_taxonomy_tables() -> None:
    assert set(FAILURE_MEANING) == set(FailureCode) and len(FAILURE_MEANING) == 11
    assert FAILURE_PRECEDENCE[0] is FailureCode.F10 and set(FAILURE_PRECEDENCE) == set(FailureCode)
    assert FAILURE_PRECEDENCE.index(FailureCode.F7) < FAILURE_PRECEDENCE.index(FailureCode.F6)
    assert set(FIELD_CODE) == set(CERTIFICATE_FIELDS) and len(CERTIFICATE_FIELDS) == 14
    assert FIELD_CODE["nontrivial"] is FailureCode.F1 and FIELD_CODE["d_weyl"] is FailureCode.F4


def test_point_verdict_consistency() -> None:
    ok = PointVerdict(_cert(), (), None, Verdict.GEOMETRIC_CANDIDATE, {"A": 1.0})
    assert ok.verdict is Verdict.GEOMETRIC_CANDIDATE
    with pytest.raises(ValueError):
        PointVerdict(_cert(connected=False), (), None, Verdict.GEOMETRIC_CANDIDATE, {})
    with pytest.raises(ValueError):
        PointVerdict(_cert(), (FailureCode.F10,), FailureCode.F10, Verdict.GEOMETRIC_CANDIDATE, {})
