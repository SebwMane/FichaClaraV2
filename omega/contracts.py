"""Contratos de datos de Omega-1.1 (OMEGA_1_1_DESIGN §3.0 y §1.12).

Dataclasses frozen+slots sin logica salvo validacion y propiedades triviales. Congelado tras WP-0.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Final, Literal

import numpy as np

from omega.config.settings11 import DistanceMode, Engine
from omega.types import (
    DimensionEstimate,
    FloatArray,
    GeometryObservables,
    IntArray,
    RunStatus,
    TopologyObservables,
)

__all__ = [
    "RULE_D3_NEVER_SUFFICIENT",
    "FailureCode",
    "FAILURE_MEANING",
    "FAILURE_PRECEDENCE",
    "Verdict",
    "CERTIFICATE_FIELDS",
    "FIELD_CODE",
    "DistanceSuiteResult",
    "HomogeneityReport",
    "IsotropyReport",
    "LocalityReport",
    "AnnulusReport",
    "BettiResult",
    "ShortCycleBetti",
    "FiltrationLevel",
    "BettiCurves",
    "PersistenceH0",
    "TopologySummary",
    "CurvatureSummary",
    "QRCProfile",
    "RunEvidence",
    "RunAssessment11",
    "GeometryCertificate",
    "PointVerdict",
    "ChainResult",
    "EnsembleSummary",
    "TransitionSummary",
    "CoarseLevel",
    "CoarseConsistency",
]

RULE_D3_NEVER_SUFFICIENT: Final = "D≈3 jamás será suficiente por sí solo para declarar geometría emergente"
"""D≈3 jamás será suficiente por sí solo para declarar geometría emergente"""


class FailureCode(Enum):
    F0 = "Ω-F0"
    F1 = "Ω-F1"
    F2 = "Ω-F2"
    F3 = "Ω-F3"
    F4 = "Ω-F4"
    F5 = "Ω-F5"
    F6 = "Ω-F6"
    F7 = "Ω-F7"
    F8 = "Ω-F8"
    F9 = "Ω-F9"
    F10 = "Ω-F10"


FAILURE_MEANING: Final[Mapping[FailureCode, str]] = MappingProxyType(
    {
        FailureCode.F0: "trivial-empty",
        FailureCode.F1: "trivial-complete",
        FailureCode.F2: "fragmented",
        FailureCode.F3: "small-world",
        FailureCode.F4: "dimension-artifact",
        FailureCode.F5: "metric-dependent",
        FailureCode.F6: "finite-size artifact",
        FailureCode.F7: "seed-dependent",
        FailureCode.F8: "null-reproducible",
        FailureCode.F9: "non-manifold",
        FailureCode.F10: "no stable phase",
    }
)

FAILURE_PRECEDENCE: Final[tuple[FailureCode, ...]] = (
    FailureCode.F10,
    FailureCode.F0,
    FailureCode.F1,
    FailureCode.F2,
    FailureCode.F3,
    FailureCode.F4,
    FailureCode.F5,
    FailureCode.F7,
    FailureCode.F6,
    FailureCode.F8,
    FailureCode.F9,
)


class Verdict(Enum):
    GEOMETRIC_CANDIDATE = "Ω-CANDIDATE"
    NOT_CANDIDATE = "NOT_CANDIDATE"


CERTIFICATE_FIELDS: Final[tuple[str, ...]] = (
    "connected",
    "nontrivial",
    "locality",
    "d_volume",
    "d_spectral",
    "d_weyl",
    "metric_robust",
    "seed_robust",
    "size_robust",
    "null_separated",
    "isotropy_ok",
    "homogeneity_ok",
    "topology_stable",
    "manifold_proxy_ok",
)

FIELD_CODE: Final[Mapping[str, FailureCode]] = MappingProxyType(
    {
        "connected": FailureCode.F2,
        "nontrivial": FailureCode.F1,
        "locality": FailureCode.F3,
        "d_volume": FailureCode.F4,
        "d_spectral": FailureCode.F4,
        "d_weyl": FailureCode.F4,
        "metric_robust": FailureCode.F5,
        "size_robust": FailureCode.F6,
        "seed_robust": FailureCode.F7,
        "null_separated": FailureCode.F8,
        "isotropy_ok": FailureCode.F9,
        "homogeneity_ok": FailureCode.F9,
        "topology_stable": FailureCode.F9,
        "manifold_proxy_ok": FailureCode.F9,
    }
)


# ---------------------------------------------------------------- validacion

def _bool(obj: object, name: str) -> None:
    if type(getattr(obj, name)) is not bool:
        raise TypeError(f"{name} debe ser bool estricto, recibido {type(getattr(obj, name)).__name__}")


def _int(obj: object, name: str, minimum: int | None = None) -> None:
    v = getattr(obj, name)
    if isinstance(v, bool) or not isinstance(v, int):
        raise TypeError(f"{name} debe ser entero, recibido {type(v).__name__}")
    if minimum is not None and v < minimum:
        raise ValueError(f"{name} debe ser >= {minimum}, recibido {v}")


def _opt_int(obj: object, name: str, minimum: int | None = None) -> None:
    if getattr(obj, name) is not None:
        _int(obj, name, minimum)


def _float(obj: object, name: str, *, finite: bool = False, minimum: float | None = None,
           maximum: float | None = None) -> None:
    """Real (no bool); NaN admitido salvo finite=True. Normaliza ints a float."""
    v = getattr(obj, name)
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise TypeError(f"{name} debe ser real, recibido {type(v).__name__}")
    x = float(v)
    if finite and not math.isfinite(x):
        raise ValueError(f"{name} debe ser finito, recibido {x}")
    if not math.isnan(x):
        if minimum is not None and x < minimum:
            raise ValueError(f"{name} debe ser >= {minimum}, recibido {x}")
        if maximum is not None and x > maximum:
            raise ValueError(f"{name} debe ser <= {maximum}, recibido {x}")
    object.__setattr__(obj, name, x)


def _isinst(obj: object, name: str, cls: type | tuple[type, ...]) -> None:
    if not isinstance(getattr(obj, name), cls):
        want = cls.__name__ if isinstance(cls, type) else "/".join(c.__name__ for c in cls)
        raise TypeError(f"{name} debe ser {want}, recibido {type(getattr(obj, name)).__name__}")


def _status(obj: object, name: str) -> None:
    if getattr(obj, name) not in ("ok", "over_budget"):
        raise ValueError(f"{name} debe ser 'ok' u 'over_budget', recibido {getattr(obj, name)!r}")


def _check_array(obj: object, name: str, kind: str, ndim: int = 1) -> None:
    a = getattr(obj, name)
    if not isinstance(a, np.ndarray):
        raise TypeError(f"{name} debe ser ndarray")
    if a.ndim != ndim:
        raise ValueError(f"{name} debe tener ndim={ndim}, recibido {a.ndim}")
    if kind == "float" and a.dtype != np.float64:
        raise TypeError(f"{name} debe ser float64, recibido {a.dtype}")
    if kind == "int" and a.dtype != np.int64:
        raise TypeError(f"{name} debe ser int64, recibido {a.dtype}")


def _mapping(obj: object, name: str, key_type: type | tuple[type, ...]) -> None:
    m = getattr(obj, name)
    if not isinstance(m, Mapping):
        raise TypeError(f"{name} debe ser Mapping")
    for k in m:
        if not isinstance(k, key_type):
            raise TypeError(f"claves de {name} con tipo invalido: {k!r}")


def _all_float_mapping(obj: object, name: str) -> None:
    _mapping(obj, name, str)
    for k, v in getattr(obj, name).items():
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            raise TypeError(f"{name}[{k!r}] debe ser real")


# ---------------------------------------------------------------- geometria

@dataclass(frozen=True, slots=True)
class DistanceSuiteResult:
    estimates: Mapping[DistanceMode, DimensionEstimate]
    fallback_used: Mapping[DistanceMode, bool]
    resistance_exponent: float
    metric_spread: float
    distance_sensitive: bool
    reason: str

    def __post_init__(self) -> None:
        _mapping(self, "estimates", DistanceMode)
        _mapping(self, "fallback_used", DistanceMode)
        if any(not isinstance(v, DimensionEstimate) for v in self.estimates.values()):
            raise TypeError("estimates debe mapear a DimensionEstimate")
        if any(type(v) is not bool for v in self.fallback_used.values()):
            raise TypeError("fallback_used debe mapear a bool")
        _float(self, "resistance_exponent")
        _float(self, "metric_spread")
        _bool(self, "distance_sensitive")
        _isinst(self, "reason", str)


@dataclass(frozen=True, slots=True)
class HomogeneityReport:
    mu: float
    sigma: float
    cv: float
    n_valid: int
    ok: bool

    def __post_init__(self) -> None:
        _float(self, "mu")
        _float(self, "sigma", minimum=0.0)
        _float(self, "cv", minimum=0.0)
        _int(self, "n_valid", 0)
        _bool(self, "ok")


@dataclass(frozen=True, slots=True)
class IsotropyReport:
    median_ratio: float
    p10_ratio: float
    radius: int
    k: int
    n_sources: int
    ok: bool

    def __post_init__(self) -> None:
        _float(self, "median_ratio", minimum=0.0)
        _float(self, "p10_ratio", minimum=0.0)
        _int(self, "radius", 0)
        _int(self, "k", 0)
        _int(self, "n_sources", 0)
        _bool(self, "ok")


@dataclass(frozen=True, slots=True)
class LocalityReport:
    detour_fraction: float
    ok: bool

    def __post_init__(self) -> None:
        _float(self, "detour_fraction", minimum=0.0, maximum=1.0)
        _bool(self, "ok")


@dataclass(frozen=True, slots=True)
class AnnulusReport:
    """Enmienda A-17 (D1/B4): `excluded_radius` = r_hi (radio contaminado por saturación), None si no se excluyó.

    `sensitivity` es diagnóstico NO decisorio (nunca entra en taxonomía/certificado).
    """

    fractions: Mapping[int, float]
    ok: bool
    excluded_radius: int | None = None
    evaluated_radii: tuple[int, ...] = ()
    sensitivity: Mapping[str, float | int | bool | None] = field(default_factory=dict)

    def __post_init__(self) -> None:
        _mapping(self, "fractions", int)
        _opt_int(self, "excluded_radius")
        if not isinstance(self.evaluated_radii, tuple) or not all(
            isinstance(r, int) and not isinstance(r, bool) for r in self.evaluated_radii
        ):
            raise TypeError("evaluated_radii debe ser tuple[int, ...]")
        _mapping(self, "sensitivity", str)
        for k, v in self.fractions.items():
            if isinstance(k, bool) or isinstance(v, bool) or not isinstance(v, (int, float)):
                raise TypeError("fractions debe ser Mapping[int, float]")
        _bool(self, "ok")


# ---------------------------------------------------------------- topologia

@dataclass(frozen=True, slots=True)
class BettiResult:
    betti: tuple[int, ...]
    counts: tuple[int, ...]
    euler_betti: int
    euler_counts: int
    status: Literal["ok", "over_budget"]
    field: str = "Z2"

    def __post_init__(self) -> None:
        for name in ("betti", "counts"):
            v = getattr(self, name)
            if not isinstance(v, tuple) or any(isinstance(x, bool) or not isinstance(x, int) or x < 0 for x in v):
                raise ValueError(f"{name} debe ser tuple de enteros >= 0")
        _int(self, "euler_betti")
        _int(self, "euler_counts")
        _status(self, "status")
        _isinst(self, "field", str)


@dataclass(frozen=True, slots=True)
class ShortCycleBetti:
    b0: int
    b1: int
    n_edges: int
    n_faces: int
    b1_density: float
    max_length: int
    status: Literal["ok", "over_budget"]

    def __post_init__(self) -> None:
        for name in ("b0", "b1", "n_edges", "n_faces", "max_length"):
            _int(self, name, 0)
        _float(self, "b1_density", minimum=0.0)
        _status(self, "status")


@dataclass(frozen=True, slots=True)
class FiltrationLevel:
    theta: float
    n_edges: int
    n_components: int
    giant_fraction: float
    mean_degree: float

    def __post_init__(self) -> None:
        _float(self, "theta", finite=True, minimum=0.0, maximum=1.0)
        _int(self, "n_edges", 0)
        _int(self, "n_components", 0)
        _float(self, "giant_fraction", minimum=0.0, maximum=1.0)
        _float(self, "mean_degree", minimum=0.0)


@dataclass(frozen=True, slots=True)
class BettiCurves:
    thetas: FloatArray
    beta0: IntArray
    beta1_short: IntArray
    b1_density: FloatArray
    giant_fraction: FloatArray
    status: tuple[str, ...]

    def __post_init__(self) -> None:
        _check_array(self, "thetas", "float")
        _check_array(self, "beta0", "int")
        _check_array(self, "beta1_short", "int")
        _check_array(self, "b1_density", "float")
        _check_array(self, "giant_fraction", "float")
        n = self.thetas.shape[0]
        for name in ("beta0", "beta1_short", "b1_density", "giant_fraction"):
            if getattr(self, name).shape[0] != n:
                raise ValueError(f"{name} debe tener la misma longitud que thetas")
        if not isinstance(self.status, tuple) or len(self.status) != n:
            raise ValueError("status debe ser tuple con la longitud de thetas")


@dataclass(frozen=True, slots=True)
class PersistenceH0:
    births: FloatArray
    deaths: FloatArray
    n_essential: int

    def __post_init__(self) -> None:
        _check_array(self, "births", "float")
        _check_array(self, "deaths", "float")
        if self.births.shape != self.deaths.shape:
            raise ValueError("births y deaths deben tener la misma forma")
        _int(self, "n_essential", 0)


@dataclass(frozen=True, slots=True)
class TopologySummary:
    curves: BettiCurves
    h0: PersistenceH0
    at_w_min: ShortCycleBetti
    clique_at_w_min: BettiResult
    stable_flags: tuple[bool, ...]
    stable: bool

    def __post_init__(self) -> None:
        _isinst(self, "curves", BettiCurves)
        _isinst(self, "h0", PersistenceH0)
        _isinst(self, "at_w_min", ShortCycleBetti)
        _isinst(self, "clique_at_w_min", BettiResult)
        if not isinstance(self.stable_flags, tuple) or any(type(f) is not bool for f in self.stable_flags):
            raise TypeError("stable_flags debe ser tuple de bool")
        _bool(self, "stable")


@dataclass(frozen=True, slots=True)
class CurvatureSummary:
    mean: float
    std: float
    se: float
    tail_fraction: float
    n_edges: int
    sampled: bool
    edge_values: FloatArray
    ok: bool

    def __post_init__(self) -> None:
        _float(self, "mean")
        _float(self, "std", minimum=0.0)
        _float(self, "se", minimum=0.0)
        _float(self, "tail_fraction", minimum=0.0, maximum=1.0)
        _int(self, "n_edges", 0)
        _bool(self, "sampled")
        _check_array(self, "edge_values", "float")
        _bool(self, "ok")


@dataclass(frozen=True, slots=True)
class QRCProfile:
    deltas: IntArray
    ratio_mean: FloatArray
    ratio_se: FloatArray
    n_pairs: IntArray

    def __post_init__(self) -> None:
        _check_array(self, "deltas", "int")
        _check_array(self, "ratio_mean", "float")
        _check_array(self, "ratio_se", "float")
        _check_array(self, "n_pairs", "int")
        n = self.deltas.shape[0]
        for name in ("ratio_mean", "ratio_se", "n_pairs"):
            if getattr(self, name).shape[0] != n:
                raise ValueError(f"{name} debe tener la misma longitud que deltas")


# ---------------------------------------------------------------- corrida

@dataclass(frozen=True, slots=True)
class RunEvidence:
    status: RunStatus
    topology: TopologyObservables
    geometry: GeometryObservables
    suite: DistanceSuiteResult
    d_weyl: DimensionEstimate
    homogeneity: HomogeneityReport
    isotropy: IsotropyReport
    locality: LocalityReport
    annulus: AnnulusReport
    topo: TopologySummary
    curvature: CurvatureSummary
    qrc: QRCProfile | None
    weight_mean: float
    weight_cv: float
    consensus_dimension: float
    dimension_class: int | None
    clustering_ratio: float
    fiedler_length: float
    n: int

    def __post_init__(self) -> None:
        for name, cls in (
            ("status", RunStatus), ("topology", TopologyObservables), ("geometry", GeometryObservables),
            ("suite", DistanceSuiteResult), ("d_weyl", DimensionEstimate), ("homogeneity", HomogeneityReport),
            ("isotropy", IsotropyReport), ("locality", LocalityReport), ("annulus", AnnulusReport),
            ("topo", TopologySummary), ("curvature", CurvatureSummary),
        ):
            _isinst(self, name, cls)
        if self.qrc is not None:
            _isinst(self, "qrc", QRCProfile)
        _float(self, "weight_mean")
        _float(self, "weight_cv")
        _float(self, "consensus_dimension")
        _opt_int(self, "dimension_class")
        _float(self, "clustering_ratio")
        _float(self, "fiedler_length")
        _int(self, "n", 2)


@dataclass(frozen=True, slots=True)
class RunAssessment11:
    codes: tuple[FailureCode, ...]
    primary: FailureCode | None
    flags: Mapping[str, bool]
    passes: bool

    def __post_init__(self) -> None:
        if not isinstance(self.codes, tuple) or any(not isinstance(c, FailureCode) for c in self.codes):
            raise TypeError("codes debe ser tuple de FailureCode")
        if self.primary is not None:
            _isinst(self, "primary", FailureCode)
            if self.primary not in self.codes:
                raise ValueError("primary debe pertenecer a codes")
        elif self.codes:
            raise ValueError("primary None exige codes vacio")
        _mapping(self, "flags", str)
        if any(type(v) is not bool for v in self.flags.values()):
            raise TypeError("flags debe mapear a bool")
        _bool(self, "passes")
        if self.passes and self.codes:
            raise ValueError("passes=True exige codes vacio")


@dataclass(frozen=True, slots=True)
class GeometryCertificate:
    connected: bool
    nontrivial: bool
    locality: bool
    d_volume: bool
    d_spectral: bool
    d_weyl: bool
    metric_robust: bool
    seed_robust: bool
    size_robust: bool
    null_separated: bool
    isotropy_ok: bool
    homogeneity_ok: bool
    topology_stable: bool
    manifold_proxy_ok: bool
    consensus_dimension: float
    dimension_class: int | None
    n_runs: int
    evidence: Mapping[str, float | int | str | bool | None]

    def __post_init__(self) -> None:
        for name in CERTIFICATE_FIELDS:
            _bool(self, name)  # type(x) is bool: rechaza np.bool_
        _float(self, "consensus_dimension")
        _opt_int(self, "dimension_class")
        _int(self, "n_runs", 0)
        _mapping(self, "evidence", str)
        for k, v in self.evidence.items():
            if v is not None and not isinstance(v, (float, int, str, bool)):
                raise TypeError(f"evidence[{k!r}] debe ser float|int|str|bool|None")

    @property
    def satisfied_all(self) -> bool:
        return all(getattr(self, name) for name in CERTIFICATE_FIELDS)

    @property
    def failed_fields(self) -> tuple[str, ...]:
        return tuple(name for name in CERTIFICATE_FIELDS if not getattr(self, name))

    @property
    def verdict(self) -> Verdict:
        if self.satisfied_all and self.dimension_class is not None:
            return Verdict.GEOMETRIC_CANDIDATE
        return Verdict.NOT_CANDIDATE


@dataclass(frozen=True, slots=True)
class PointVerdict:
    certificate: GeometryCertificate
    codes: tuple[FailureCode, ...]
    primary: FailureCode | None
    verdict: Verdict
    run_outcome_fractions: Mapping[str, float]

    def __post_init__(self) -> None:
        _isinst(self, "certificate", GeometryCertificate)
        if not isinstance(self.codes, tuple) or any(not isinstance(c, FailureCode) for c in self.codes):
            raise TypeError("codes debe ser tuple de FailureCode")
        if self.primary is not None and self.primary not in self.codes:
            raise ValueError("primary debe pertenecer a codes")
        _isinst(self, "verdict", Verdict)
        if self.verdict is Verdict.GEOMETRIC_CANDIDATE and (
            self.certificate.verdict is not Verdict.GEOMETRIC_CANDIDATE or FailureCode.F10 in self.codes
        ):
            raise ValueError("GEOMETRIC_CANDIDATE exige certificado completo y ausencia de F10")
        _all_float_mapping(self, "run_outcome_fractions")


# ---------------------------------------------------------------- estadistica y coarse-graining

@dataclass(frozen=True, slots=True)
class ChainResult:
    engine: Engine
    theta: float
    samples: Mapping[str, FloatArray]
    sample_steps: IntArray
    w_final: FloatArray
    states: tuple[FloatArray, ...]
    acceptance: float
    step_size: float
    n_steps: int

    def __post_init__(self) -> None:
        _isinst(self, "engine", Engine)
        _float(self, "theta", finite=True, minimum=0.0)
        _mapping(self, "samples", str)
        _check_array(self, "sample_steps", "int")
        _check_array(self, "w_final", "float", ndim=2)
        if not isinstance(self.states, tuple):
            raise TypeError("states debe ser tuple")
        _float(self, "acceptance", minimum=0.0, maximum=1.0)
        _float(self, "step_size", minimum=0.0)
        _int(self, "n_steps", 0)


@dataclass(frozen=True, slots=True)
class EnsembleSummary:
    means: Mapping[str, float]
    ses: Mapping[str, float]
    rhat: Mapping[str, float]
    tau_int: Mapping[str, float]
    ess: Mapping[str, float]
    geweke_ok: bool
    equilibrated: bool

    def __post_init__(self) -> None:
        for name in ("means", "ses", "rhat", "tau_int", "ess"):
            _all_float_mapping(self, name)
        _bool(self, "geweke_ok")
        _bool(self, "equilibrated")


@dataclass(frozen=True, slots=True)
class TransitionSummary:
    lambdas: FloatArray
    mean: FloatArray
    derivative: FloatArray
    susceptibility: FloatArray
    binder: FloatArray
    bimodality: FloatArray
    peak_lambda: float

    def __post_init__(self) -> None:
        names = ("lambdas", "mean", "derivative", "susceptibility", "binder", "bimodality")
        for name in names:
            _check_array(self, name, "float")
        n = self.lambdas.shape[0]
        for name in names[1:]:
            if getattr(self, name).shape[0] != n:
                raise ValueError(f"{name} debe tener la misma longitud que lambdas")
        _float(self, "peak_lambda")


@dataclass(frozen=True, slots=True)
class CoarseLevel:
    level: int
    w: FloatArray
    blocks: tuple[tuple[int, ...], ...]
    n: int

    def __post_init__(self) -> None:
        _int(self, "level", 0)
        _check_array(self, "w", "float", ndim=2)
        _int(self, "n", 1)
        if self.w.shape != (self.n, self.n):
            raise ValueError("w debe ser (n, n)")
        if not isinstance(self.blocks, tuple):
            raise TypeError("blocks debe ser tuple")


@dataclass(frozen=True, slots=True)
class CoarseConsistency:
    classes: tuple[int | None, ...]
    codes: tuple[FailureCode | None, ...]
    measured_levels: int
    undetermined_levels: int
    consistent: bool

    def __post_init__(self) -> None:
        if not isinstance(self.classes, tuple) or any(
            c is not None and (isinstance(c, bool) or not isinstance(c, int)) for c in self.classes
        ):
            raise TypeError("classes debe ser tuple de int|None")
        if not isinstance(self.codes, tuple) or any(c is not None and not isinstance(c, FailureCode) for c in self.codes):
            raise TypeError("codes debe ser tuple de FailureCode|None")
        _int(self, "measured_levels", 0)
        _int(self, "undetermined_levels", 0)
        _bool(self, "consistent")

