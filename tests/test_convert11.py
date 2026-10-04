"""Pruebas de convert11: ida y vuelta y hash canónico estable (WP-E)."""

from __future__ import annotations

import dataclasses
import json

import pytest

from omega.config.convert11 import canonical_json, config11_from_dict, config11_to_dict, config_hash
from omega.config.settings11 import (
    CurvatureConfig,
    DistanceMode,
    Engine,
    LangevinConfig,
    NullModel,
    Omega11Config,
    TopologyConfig,
)
from tests.conftest import make_config


def _cfg() -> Omega11Config:
    return Omega11Config(base=make_config(12))


def test_roundtrip_default() -> None:
    cfg = _cfg()
    d = config11_to_dict(cfg)
    json.dumps(d, allow_nan=False)
    assert config11_from_dict(json.loads(json.dumps(d))) == cfg


def test_roundtrip_enums_and_optionals() -> None:
    base = make_config(12)
    cfg = Omega11Config(
        base=base,
        curvature=CurvatureConfig(distance_mode=DistanceMode.WEIGHTED_LOG),
        engine=Engine.LANGEVIN,
        langevin=LangevinConfig(theta_hat=0.5),
        null_models=(NullModel.SMALL_WORLD, NullModel.RANDOM_GEOMETRIC),
    )
    d = json.loads(canonical_json(config11_to_dict(cfg)))
    assert d["engine"] == "langevin"
    assert config11_from_dict(d) == cfg


def test_rejects_unknown_keys() -> None:
    d = config11_to_dict(_cfg())
    d["foo"] = 1
    with pytest.raises(ValueError):
        config11_from_dict(d)
    d = config11_to_dict(_cfg())
    d["weyl"]["bogus"] = 1
    with pytest.raises(ValueError):
        config11_from_dict(d)
    d = config11_to_dict(_cfg())
    d["schema_version"] = "1.0"
    with pytest.raises(ValueError):
        config11_from_dict(d)


def test_hash_stable_under_key_order_and_sensitive() -> None:
    cfg = _cfg()
    d = config11_to_dict(cfg)
    rev = dict(reversed(list(d.items())))
    assert canonical_json(d) == canonical_json(rev)
    h = config_hash(cfg)
    assert h == config_hash(config11_from_dict(rev))
    assert len(h) == 64
    changed = dataclasses.replace(cfg, topology=TopologyConfig(b1_density_max=0.04))
    assert config_hash(changed) != h
