"""Dataclasses de configuracion (ANALYSIS §5.1). Frozen, slots y validadas en __post_init__."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

__all__ = [
    "FunctionalParams",
    "InitConfig",
    "DynamicsConfig",
    "GraphConfig",
    "DimensionConfig",
    "SpectralConfig",
    "PhaseThresholds",
    "SeedConfig",
    "ScanConfig",
    "OmegaConfig",
]


def _real(obj: object, name: str, value: object) -> float:
    """Valida que `value` sea real finito (no bool) y lo fija como float."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} debe ser un numero real, recibido {type(value).__name__}")
    x = float(value)
    if not math.isfinite(x):
        raise ValueError(f"{name} debe ser finito, recibido {x}")
    object.__setattr__(obj, name, x)
    return x


def _integer(name: str, value: object, minimum: int) -> int:
    """Valida entero (no bool) >= minimum."""
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} debe ser entero, recibido {type(value).__name__}")
    if value < minimum:
        raise ValueError(f"{name} debe ser >= {minimum}, recibido {value}")
    return value


def _choice(name: str, value: object, options: tuple[str, ...]) -> None:
    if value not in options:
        raise ValueError(f"{name} debe ser uno de {options}, recibido {value!r}")


def _float_tuple(obj: object, name: str, value: object, minimum: float = 0.0) -> None:
    """Valida tupla no vacia de reales finitos >= minimum; la normaliza a floats."""
    if not isinstance(value, tuple):
        raise TypeError(f"{name} debe ser tuple, recibido {type(value).__name__}")
    if len(value) == 0:
        raise ValueError(f"{name} no puede estar vacia")
    out: list[float] = []
    for i, v in enumerate(value):
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
            raise ValueError(f"{name}[{i}] debe ser real finito, recibido {v!r}")
        if v < minimum:
            raise ValueError(f"{name}[{i}] debe ser >= {minimum}, recibido {v}")
        out.append(float(v))
    object.__setattr__(obj, name, tuple(out))


@dataclass(frozen=True, slots=True)
class FunctionalParams:
    """Parametros crudos de S0 (M§13-M§15). beta>0, gamma>=0, eta>=0; alpha y mu reales."""

    alpha: float
    beta: float = 1.0
    gamma: float = 0.0
    eta: float = 0.0
    mu: float = 0.0

    def __post_init__(self) -> None:
        _real(self, "alpha", self.alpha)
        if _real(self, "beta", self.beta) <= 0.0:
            raise ValueError(f"beta debe ser > 0, recibido {self.beta}")
        if _real(self, "gamma", self.gamma) < 0.0:
            raise ValueError(f"gamma debe ser >= 0, recibido {self.gamma}")
        if _real(self, "eta", self.eta) < 0.0:
            raise ValueError(f"eta debe ser >= 0, recibido {self.eta}")
        _real(self, "mu", self.mu)


@dataclass(frozen=True, slots=True)
class InitConfig:
    """Inicializacion de W0 (D-0, D-9): n>=2 nodos."""

    n: int
    distribution: Literal["uniform"] = "uniform"
    symmetrization: Literal["upper_mirror"] = "upper_mirror"

    def __post_init__(self) -> None:
        _integer("n", self.n, 2)
        _choice("distribution", self.distribution, ("uniform",))
        _choice("symmetrization", self.symmetrization, ("upper_mirror",))


@dataclass(frozen=True, slots=True)
class DynamicsConfig:
    """Integrador y criterios de parada (D-3, D-4, D-13)."""

    integrator: Literal["clip", "sigmoid"] = "clip"
    dt_mode: Literal["auto", "fixed"] = "auto"
    dt_safety: float = 0.5
    dt_fixed: float | None = None
    max_steps: int = 200_000
    tol_step: float = 1e-10
    patience: int = 10
    snapshot_schedule: Literal["log2", "none"] = "log2"
    sigmoid_clip: float = 1e-6

    def __post_init__(self) -> None:
        _choice("integrator", self.integrator, ("clip", "sigmoid"))
        _choice("dt_mode", self.dt_mode, ("auto", "fixed"))
        _choice("snapshot_schedule", self.snapshot_schedule, ("log2", "none"))
        s = _real(self, "dt_safety", self.dt_safety)
        if not 0.0 < s < 1.0:
            raise ValueError(f"dt_safety debe estar en (0,1), recibido {s}")
        if self.dt_fixed is not None:
            if _real(self, "dt_fixed", self.dt_fixed) <= 0.0:
                raise ValueError(f"dt_fixed debe ser > 0, recibido {self.dt_fixed}")
        elif self.dt_mode == "fixed":
            raise ValueError("dt_mode='fixed' exige dt_fixed")
        _integer("max_steps", self.max_steps, 1)
        if _real(self, "tol_step", self.tol_step) <= 0.0:
            raise ValueError(f"tol_step debe ser > 0, recibido {self.tol_step}")
        _integer("patience", self.patience, 1)
        c = _real(self, "sigmoid_clip", self.sigmoid_clip)
        if not 0.0 < c < 0.5:
            raise ValueError(f"sigmoid_clip debe estar en (0,0.5), recibido {c}")


@dataclass(frozen=True, slots=True)
class GraphConfig:
    """Grafo derivado de W (D-1, D-2): umbral w_min, epsilon y barrido de sensibilidad."""

    w_min: float = 0.1
    epsilon: float = 1e-9
    w_min_sensitivity: tuple[float, ...] = (0.01, 0.05, 0.1, 0.2, 0.5)

    def __post_init__(self) -> None:
        w = _real(self, "w_min", self.w_min)
        if not 0.0 <= w <= 1.0:
            raise ValueError(f"w_min debe estar en [0,1], recibido {w}")
        if _real(self, "epsilon", self.epsilon) <= 0.0:
            raise ValueError(f"epsilon debe ser > 0, recibido {self.epsilon}")
        _float_tuple(self, "w_min_sensitivity", self.w_min_sensitivity)
        if any(v > 1.0 for v in self.w_min_sensitivity):
            raise ValueError("w_min_sensitivity debe estar en [0,1]")


@dataclass(frozen=True, slots=True)
class DimensionConfig:
    """Estimacion de D_eff (D-6, D-19, D-22)."""

    estimator: Literal["shell", "ball"] = "shell"
    n_radii: int = 40
    saturation: float = 0.2
    min_scale_ratio: float = 3.0
    min_points: int = 3
    plateau_tol: float = 0.3
    integer_tol: float = 1e-6

    def __post_init__(self) -> None:
        _choice("estimator", self.estimator, ("shell", "ball"))
        _integer("n_radii", self.n_radii, 2)
        s = _real(self, "saturation", self.saturation)
        if not 0.0 < s <= 1.0:
            raise ValueError(f"saturation debe estar en (0,1], recibido {s}")
        if _real(self, "min_scale_ratio", self.min_scale_ratio) <= 1.0:
            raise ValueError(f"min_scale_ratio debe ser > 1, recibido {self.min_scale_ratio}")
        _integer("min_points", self.min_points, 2)
        if _real(self, "plateau_tol", self.plateau_tol) <= 0.0:
            raise ValueError(f"plateau_tol debe ser > 0, recibido {self.plateau_tol}")
        if _real(self, "integer_tol", self.integer_tol) <= 0.0:
            raise ValueError(f"integer_tol debe ser > 0, recibido {self.integer_tol}")


@dataclass(frozen=True, slots=True)
class SpectralConfig:
    """Estimacion de D_s (D-7, D-12)."""

    method: Literal["lazy_walk", "heat_normalized"] = "lazy_walk"
    laziness: float = 0.5
    t_min: float = 4.0
    saturation_factor: float = 5.0
    min_scale_ratio: float = 3.0
    n_times: int = 200
    graph: Literal["thresholded_weighted", "weighted", "binary"] = "thresholded_weighted"

    def __post_init__(self) -> None:
        _choice("method", self.method, ("lazy_walk", "heat_normalized"))
        _choice("graph", self.graph, ("thresholded_weighted", "weighted", "binary"))
        q = _real(self, "laziness", self.laziness)
        if not 0.0 < q <= 1.0:
            raise ValueError(f"laziness debe estar en (0,1], recibido {q}")
        if _real(self, "t_min", self.t_min) <= 0.0:
            raise ValueError(f"t_min debe ser > 0, recibido {self.t_min}")
        if _real(self, "saturation_factor", self.saturation_factor) <= 0.0:
            raise ValueError(f"saturation_factor debe ser > 0, recibido {self.saturation_factor}")
        if _real(self, "min_scale_ratio", self.min_scale_ratio) <= 1.0:
            raise ValueError(f"min_scale_ratio debe ser > 1, recibido {self.min_scale_ratio}")
        _integer("n_times", self.n_times, 2)


@dataclass(frozen=True, slots=True)
class PhaseThresholds:
    """Umbrales preregistrados de clasificacion de fases (ANALYSIS §4.1, D-17..D-20)."""

    g_dispersed: float = 0.1
    kbin_dispersed: float = 1.0
    large_component_frac: float = 0.1
    g_connected: float = 0.9
    g_geometric: float = 0.95
    cv_homogeneous: float = 0.2
    rho_hyperdense: float = 0.5
    meanw_hyperdense: float = 0.5
    c_clustered: float = 0.3
    c_ratio_clustered: float = 3.0
    ds_deff_tol: float = 0.5
    f_seed_fraction: float = 0.8
    f_seed_std: float = 0.2
    f_size_tol: float = 0.3
    param_perturbation: float = 0.1

    def __post_init__(self) -> None:
        unit = (
            "g_dispersed", "large_component_frac", "g_connected", "g_geometric",
            "rho_hyperdense", "meanw_hyperdense", "c_clustered", "f_seed_fraction",
        )
        nonneg = (
            "kbin_dispersed", "cv_homogeneous", "c_ratio_clustered", "ds_deff_tol",
            "f_seed_std", "f_size_tol", "param_perturbation",
        )
        for name in unit:
            v = _real(self, name, getattr(self, name))
            if not 0.0 <= v <= 1.0:
                raise ValueError(f"{name} debe estar en [0,1], recibido {v}")
        for name in nonneg:
            v = _real(self, name, getattr(self, name))
            if v < 0.0:
                raise ValueError(f"{name} debe ser >= 0, recibido {v}")
        if self.g_geometric < self.g_connected:
            raise ValueError("g_geometric debe ser >= g_connected")


@dataclass(frozen=True, slots=True)
class SeedConfig:
    """Semillas (D-5, D-28): entropia maestra, identificador de experimento y replicas."""

    master_entropy: int
    experiment_id: int
    replicates: int = 10

    def __post_init__(self) -> None:
        _integer("master_entropy", self.master_entropy, 0)
        _integer("experiment_id", self.experiment_id, 0)
        _integer("replicates", self.replicates, 1)


@dataclass(frozen=True, slots=True)
class ScanConfig:
    """Malla de barrido en parametros reducidos (D-8, D-15)."""

    alpha_hat: tuple[float, ...]
    gamma_hat: tuple[float, ...]
    scaling: Literal["reduced", "raw"] = "reduced"

    def __post_init__(self) -> None:
        _float_tuple(self, "alpha_hat", self.alpha_hat)
        _float_tuple(self, "gamma_hat", self.gamma_hat)
        _choice("scaling", self.scaling, ("reduced", "raw"))


@dataclass(frozen=True, slots=True)
class OmegaConfig:
    """Configuracion completa, congelada y guardada en el pasaporte (ANALYSIS §4)."""

    init: InitConfig
    functional: FunctionalParams
    dynamics: DynamicsConfig
    graph: GraphConfig
    dimension: DimensionConfig
    spectral: SpectralConfig
    phases: PhaseThresholds
    seeds: SeedConfig
    scan: ScanConfig | None = None
    schema_version: str = "1.0"

    def __post_init__(self) -> None:
        expected = (
            ("init", InitConfig), ("functional", FunctionalParams),
            ("dynamics", DynamicsConfig), ("graph", GraphConfig),
            ("dimension", DimensionConfig), ("spectral", SpectralConfig),
            ("phases", PhaseThresholds), ("seeds", SeedConfig),
        )
        for name, cls in expected:
            if not isinstance(getattr(self, name), cls):
                raise TypeError(f"{name} debe ser {cls.__name__}")
        if self.scan is not None and not isinstance(self.scan, ScanConfig):
            raise TypeError("scan debe ser ScanConfig o None")
        if not isinstance(self.schema_version, str) or not self.schema_version:
            raise ValueError("schema_version debe ser una cadena no vacia")
        if self.functional.eta != 0.0 and self.dynamics.dt_mode != "fixed":
            raise ValueError("eta != 0 exige dynamics.dt_mode == 'fixed' (D-3)")
