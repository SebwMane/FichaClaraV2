"""Conversión de Omega11Config a/desde dict JSON-compatible y hash canónico (§1.13, §3.5)."""

from __future__ import annotations

import dataclasses
import hashlib
import json
from collections.abc import Mapping
from enum import Enum
from typing import Any

from omega.config.convert import config_from_dict, config_to_dict
from omega.config.settings11 import (
    CertificateThresholds,
    CoarseGrainConfig,
    CurvatureConfig,
    DistanceMode,
    DistanceSuiteConfig,
    EnsembleConfig,
    Engine,
    FiniteSizeConfig,
    FixedDensityConfig,
    LangevinConfig,
    LocalStructureConfig,
    MetropolisConfig,
    NullModel,
    Omega11Config,
    TopologyConfig,
    WeylConfig,
)

__all__ = ["config11_to_dict", "config11_from_dict", "canonical_json", "config_hash"]

_SECTIONS: dict[str, type] = {
    "distance": DistanceSuiteConfig,
    "weyl": WeylConfig,
    "topology": TopologyConfig,
    "curvature": CurvatureConfig,
    "local": LocalStructureConfig,
    "certificate": CertificateThresholds,
    "coarse": CoarseGrainConfig,
    "ensemble": EnsembleConfig,
    "finite_size": FiniteSizeConfig,
}
_OPTIONAL: dict[str, type] = {
    "langevin": LangevinConfig,
    "metropolis": MetropolisConfig,
    "fixed_density": FixedDensityConfig,
}
_ENUM_FIELDS: dict[tuple[str, str], type[Enum]] = {
    ("distance", "modes"): DistanceMode,
    ("curvature", "distance_mode"): DistanceMode,
}
_TOP = {"base", "engine", "null_models", "schema_version"} | set(_SECTIONS) | set(_OPTIONAL)


def _plain(obj: Any) -> Any:
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {f.name: _plain(getattr(obj, f.name)) for f in dataclasses.fields(obj)}
    if isinstance(obj, Enum):
        return obj.value
    if isinstance(obj, (tuple, list)):
        return [_plain(v) for v in obj]
    return obj


def config11_to_dict(cfg: Omega11Config) -> dict[str, Any]:
    """Serializa a dict JSON-compatible (enums -> value, tuplas -> listas); base vía Ω-1.0."""
    if not isinstance(cfg, Omega11Config):
        raise TypeError("cfg debe ser Omega11Config")
    out: dict[str, Any] = {"base": config_to_dict(cfg.base)}
    for f in dataclasses.fields(cfg):
        if f.name != "base":
            out[f.name] = _plain(getattr(cfg, f.name))
    return out


def _build(cls: type, d: Any, where: str) -> Any:
    if not isinstance(d, Mapping):
        raise TypeError(f"{where} debe ser dict, recibido {type(d).__name__}")
    names = {f.name for f in dataclasses.fields(cls)}
    unknown = set(d) - names
    if unknown:
        raise ValueError(f"claves desconocidas en {where}: {sorted(unknown)}")
    kwargs: dict[str, Any] = {}
    for k, v in d.items():
        enum_cls = _ENUM_FIELDS.get((where, k))
        if enum_cls is not None:
            v = tuple(enum_cls(x) for x in v) if isinstance(v, list) else enum_cls(v)
        elif isinstance(v, list):
            v = tuple(v)
        kwargs[k] = v
    return cls(**kwargs)


def config11_from_dict(d: Mapping[str, Any]) -> Omega11Config:
    """Reconstruye Omega11Config; rechaza claves desconocidas y schema_version != '1.1'."""
    if not isinstance(d, Mapping):
        raise TypeError("d debe ser dict")
    unknown = set(d) - _TOP
    if unknown:
        raise ValueError(f"claves desconocidas en la configuración 1.1: {sorted(unknown)}")
    if "base" not in d:
        raise ValueError("falta la sección 'base'")
    if d.get("schema_version", "1.1") != "1.1":
        raise ValueError(f"schema_version no soportado: {d.get('schema_version')!r}")
    kwargs: dict[str, Any] = {"base": config_from_dict(dict(d["base"]))}
    for name, cls in _SECTIONS.items():
        if name in d:
            kwargs[name] = _build(cls, d[name], name)
    for name, cls in _OPTIONAL.items():
        if d.get(name) is not None:
            kwargs[name] = _build(cls, d[name], name)
    if "engine" in d:
        kwargs["engine"] = Engine(d["engine"])
    if "null_models" in d:
        kwargs["null_models"] = tuple(NullModel(x) for x in d["null_models"])
    return Omega11Config(**kwargs)


def canonical_json(obj: Any) -> str:
    """JSON canónico: claves ordenadas, sin espacios, UTF-8 literal, sin NaN/inf."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)


def config_hash(cfg: Omega11Config) -> str:
    """sha256 hexadecimal de canonical_json(config11_to_dict(cfg))."""
    return hashlib.sha256(canonical_json(config11_to_dict(cfg)).encode("utf-8")).hexdigest()
