"""Conversion de configuracion a/desde dict (JSON-compatible) y parametros reducidos (D-15)."""

from __future__ import annotations

import dataclasses
from typing import Any

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

__all__ = ["config_to_dict", "config_from_dict", "reduced_to_raw", "raw_to_reduced"]

_SECTIONS: dict[str, type] = {
    "init": InitConfig,
    "functional": FunctionalParams,
    "dynamics": DynamicsConfig,
    "graph": GraphConfig,
    "dimension": DimensionConfig,
    "spectral": SpectralConfig,
    "phases": PhaseThresholds,
    "seeds": SeedConfig,
}


def _to_plain(obj: Any) -> Any:
    """Dataclasses -> dict; tuplas -> listas (compatible con JSON)."""
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {f.name: _to_plain(getattr(obj, f.name)) for f in dataclasses.fields(obj)}
    if isinstance(obj, (tuple, list)):
        return [_to_plain(v) for v in obj]
    return obj


def _build(cls: type, d: Any, where: str) -> Any:
    """Construye `cls` desde dict rechazando claves desconocidas; listas -> tuplas."""
    if not isinstance(d, dict):
        raise TypeError(f"{where} debe ser dict, recibido {type(d).__name__}")
    names = {f.name for f in dataclasses.fields(cls)}
    unknown = set(d) - names
    if unknown:
        raise ValueError(f"claves desconocidas en {where}: {sorted(unknown)}")
    kwargs = {k: tuple(v) if isinstance(v, list) else v for k, v in d.items()}
    return cls(**kwargs)


def config_to_dict(cfg: OmegaConfig) -> dict[str, Any]:
    """Serializa la configuracion completa a un dict JSON-compatible."""
    if not isinstance(cfg, OmegaConfig):
        raise TypeError("cfg debe ser OmegaConfig")
    out = _to_plain(cfg)
    assert isinstance(out, dict)
    return out


def config_from_dict(d: dict[str, Any]) -> OmegaConfig:
    """Reconstruye OmegaConfig; rechaza claves desconocidas y schema_version distinto de '1.0'."""
    if not isinstance(d, dict):
        raise TypeError("d debe ser dict")
    allowed = set(_SECTIONS) | {"scan", "schema_version"}
    unknown = set(d) - allowed
    if unknown:
        raise ValueError(f"claves desconocidas en la configuracion: {sorted(unknown)}")
    missing = set(_SECTIONS) - set(d)
    if missing:
        raise ValueError(f"faltan secciones: {sorted(missing)}")
    schema = d.get("schema_version", "1.0")
    if schema != "1.0":
        raise ValueError(f"schema_version no soportado: {schema!r}")
    sections = {name: _build(cls, d[name], name) for name, cls in _SECTIONS.items()}
    scan_d = d.get("scan")
    scan = None if scan_d is None else _build(ScanConfig, scan_d, "scan")
    return OmegaConfig(scan=scan, schema_version=schema, **sections)


def _check_n(n: int) -> None:
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n debe ser entero")
    if n < 3:
        raise ValueError(f"n debe ser >= 3 para la reparametrizacion, recibido {n}")


def reduced_to_raw(alpha_hat: float, gamma_hat: float, n: int, beta: float = 1.0) -> FunctionalParams:
    """alpha = 2*beta*alpha_hat/(n-2), gamma = beta*gamma_hat/(n-2) (D-15); n>=3, hats>=0."""
    _check_n(n)
    if alpha_hat < 0.0 or gamma_hat < 0.0:
        raise ValueError("alpha_hat y gamma_hat deben ser >= 0")
    return FunctionalParams(
        alpha=2.0 * beta * alpha_hat / (n - 2),
        beta=beta,
        gamma=beta * gamma_hat / (n - 2),
    )


def raw_to_reduced(p: FunctionalParams, n: int) -> tuple[float, float]:
    """Inversa de reduced_to_raw: (alpha_hat, gamma_hat) = (alpha(n-2)/(2 beta), gamma(n-2)/beta)."""
    _check_n(n)
    return (p.alpha * (n - 2) / (2.0 * p.beta), p.gamma * (n - 2) / p.beta)
