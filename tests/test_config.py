"""Pruebas de config/settings y config/convert (ANALYSIS §6.2 config)."""

from __future__ import annotations

import dataclasses
import json

import pytest

from omega.config.convert import config_from_dict, config_to_dict, raw_to_reduced, reduced_to_raw
from omega.config.settings import (
    DimensionConfig,
    DynamicsConfig,
    FunctionalParams,
    GraphConfig,
    InitConfig,
    OmegaConfig,
    PhaseThresholds,
    ScanConfig,
    SeedConfig,
    SpectralConfig,
)
from tests.conftest import make_config


def test_roundtrip_dict_and_json() -> None:
    for scan in (True, False):
        cfg = make_config(scan=scan)
        d = config_to_dict(cfg)
        assert config_from_dict(d) == cfg
        assert config_from_dict(json.loads(json.dumps(d))) == cfg


def test_to_dict_does_not_alias_and_from_dict_does_not_mutate() -> None:
    d = config_to_dict(make_config())
    snapshot = json.dumps(d, sort_keys=True)
    config_from_dict(d)
    assert json.dumps(d, sort_keys=True) == snapshot


def test_unknown_keys_and_schema_rejected() -> None:
    d = config_to_dict(make_config())
    with pytest.raises(ValueError, match="desconocidas"):
        config_from_dict({**d, "extra": 1})
    bad = json.loads(json.dumps(d))
    bad["init"]["foo"] = 1
    with pytest.raises(ValueError, match="desconocidas"):
        config_from_dict(bad)
    with pytest.raises(ValueError, match="schema_version"):
        config_from_dict({**d, "schema_version": "2.0"})
    missing = {k: v for k, v in d.items() if k != "graph"}
    with pytest.raises(ValueError, match="faltan"):
        config_from_dict(missing)


@pytest.mark.parametrize(
    "factory",
    [
        lambda: FunctionalParams(alpha=1.0, beta=0.0),
        lambda: FunctionalParams(alpha=1.0, gamma=-0.1),
        lambda: FunctionalParams(alpha=1.0, eta=-1.0),
        lambda: FunctionalParams(alpha=float("nan")),
        lambda: FunctionalParams(alpha=float("inf")),
        lambda: InitConfig(n=1),
        lambda: DynamicsConfig(dt_mode="fixed"),
        lambda: DynamicsConfig(dt_safety=1.0),
        lambda: DynamicsConfig(max_steps=0),
        lambda: DynamicsConfig(tol_step=0.0),
        lambda: DynamicsConfig(integrator="rk4"),  # type: ignore[arg-type]
        lambda: GraphConfig(w_min=1.5),
        lambda: GraphConfig(epsilon=0.0),
        lambda: GraphConfig(w_min_sensitivity=()),
        lambda: DimensionConfig(saturation=0.0),
        lambda: DimensionConfig(min_points=1),
        lambda: SpectralConfig(laziness=0.0),
        lambda: PhaseThresholds(g_dispersed=2.0),
        lambda: PhaseThresholds(g_connected=0.95, g_geometric=0.9),
        lambda: SeedConfig(master_entropy=-1, experiment_id=0),
        lambda: SeedConfig(master_entropy=1, experiment_id=0, replicates=0),
        lambda: ScanConfig(alpha_hat=(), gamma_hat=(1.0,)),
        lambda: ScanConfig(alpha_hat=(1.0,), gamma_hat=(-1.0,)),
    ],
)
def test_invalid_values_rejected(factory) -> None:  # type: ignore[no-untyped-def]
    with pytest.raises(ValueError):
        factory()


def test_type_errors() -> None:
    with pytest.raises(TypeError):
        FunctionalParams(alpha="1")  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        InitConfig(n=3.0)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        ScanConfig(alpha_hat=[1.0], gamma_hat=(1.0,))  # type: ignore[arg-type]


def test_frozen_and_ints_coerced_to_float() -> None:
    p = FunctionalParams(alpha=1, beta=2)
    assert isinstance(p.alpha, float) and isinstance(p.beta, float)
    with pytest.raises(dataclasses.FrozenInstanceError):
        p.alpha = 3.0  # type: ignore[misc]


def test_eta_requires_fixed_dt() -> None:
    base = make_config()
    with pytest.raises(ValueError, match="eta"):
        dataclasses.replace(base, functional=FunctionalParams(alpha=1.0, eta=0.1))
    ok = dataclasses.replace(
        base,
        functional=FunctionalParams(alpha=1.0, eta=0.1),
        dynamics=DynamicsConfig(dt_mode="fixed", dt_fixed=1e-3),
    )
    assert ok.functional.eta == 0.1


@pytest.mark.parametrize("n", [3, 12, 200])
@pytest.mark.parametrize("ah,gh", [(0.0, 0.0), (2.0, 10.0), (0.37, 3.3)])
def test_reduced_raw_identity(n: int, ah: float, gh: float) -> None:
    p = reduced_to_raw(ah, gh, n, beta=1.7)
    a2, g2 = raw_to_reduced(p, n)
    assert abs(a2 - ah) < 1e-12 and abs(g2 - gh) < 1e-12
    assert p.alpha == pytest.approx(2 * 1.7 * ah / (n - 2))


def test_reduced_requires_n_ge_3() -> None:
    with pytest.raises(ValueError):
        reduced_to_raw(1.0, 1.0, 2)
