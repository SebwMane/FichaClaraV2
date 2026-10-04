"""Configuracion de Omega-1.1 (OMEGA_1_1_DESIGN §3.0). Frozen, slots y validada en __post_init__.

Congelado tras WP-0: un cambio exige decision del orquestador y nueva version del documento.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Final, Literal

from omega.config.settings import OmegaConfig, _choice, _integer, _real

__all__ = [
    "DistanceMode",
    "GEODESIC_MODES",
    "NullModel",
    "Engine",
    "DistanceSuiteConfig",
    "WeylConfig",
    "TopologyConfig",
    "CurvatureConfig",
    "LocalStructureConfig",
    "CertificateThresholds",
    "CoarseGrainConfig",
    "LangevinConfig",
    "MetropolisConfig",
    "EnsembleConfig",
    "FixedDensityConfig",
    "FiniteSizeConfig",
    "Omega11Config",
]


class DistanceMode(Enum):
    HOP = "hop"
    WEIGHTED_INVERSE = "weighted_inverse"
    WEIGHTED_LOG = "weighted_log"
    RESISTANCE = "resistance"


GEODESIC_MODES: Final[tuple[DistanceMode, ...]] = (
    DistanceMode.HOP,
    DistanceMode.WEIGHTED_INVERSE,
    DistanceMode.WEIGHTED_LOG,
)


class NullModel(Enum):
    ERDOS_RENYI = "erdos_renyi"
    CONFIGURATION_MODEL = "configuration_model"
    DEGREE_PRESERVING_REWIRE = "degree_preserving_rewire"
    SHUFFLED_WEIGHTS = "shuffled_weights"
    SMALL_WORLD = "small_world"
    RANDOM_GEOMETRIC = "random_geometric"


class Engine(Enum):
    GRADIENT = "gradient"
    LANGEVIN = "langevin"
    METROPOLIS = "metropolis"
    FIXED_DENSITY = "fixed_density"


def _num(
    obj: object,
    name: str,
    lo: float | None = None,
    hi: float | None = None,
    *,
    lo_open: bool = False,
    hi_open: bool = False,
) -> float:
    """Real finito dentro de [lo, hi] (extremos abiertos si se indica)."""
    x = _real(obj, name, getattr(obj, name))
    if lo is not None and (x <= lo if lo_open else x < lo):
        raise ValueError(f"{name} debe ser {'>' if lo_open else '>='} {lo}, recibido {x}")
    if hi is not None and (x >= hi if hi_open else x > hi):
        raise ValueError(f"{name} debe ser {'<' if hi_open else '<='} {hi}, recibido {x}")
    return x


def _int(obj: object, name: str, minimum: int) -> int:
    return _integer(name, getattr(obj, name), minimum)


def _enum_tuple(obj: object, name: str, cls: type[Enum], allow_empty: bool = False) -> None:
    v = getattr(obj, name)
    if not isinstance(v, tuple):
        raise TypeError(f"{name} debe ser tuple, recibido {type(v).__name__}")
    if not allow_empty and len(v) == 0:
        raise ValueError(f"{name} no puede estar vacia")
    for item in v:
        if not isinstance(item, cls):
            raise TypeError(f"{name} debe contener {cls.__name__}, recibido {item!r}")
    if len(set(v)) != len(v):
        raise ValueError(f"{name} no puede tener duplicados")


def _tuple_of(obj: object, name: str, lo: float, hi: float, *, integer: bool = False) -> None:
    """Tupla no vacia de reales (o enteros) en [lo, hi]; normaliza a float (no a int)."""
    v = getattr(obj, name)
    if not isinstance(v, tuple):
        raise TypeError(f"{name} debe ser tuple, recibido {type(v).__name__}")
    if len(v) == 0:
        raise ValueError(f"{name} no puede estar vacia")
    for i, x in enumerate(v):
        if isinstance(x, bool) or not isinstance(x, (int, float)):
            raise TypeError(f"{name}[{i}] debe ser numerico, recibido {x!r}")
        if integer and not isinstance(x, int):
            raise TypeError(f"{name}[{i}] debe ser entero, recibido {x!r}")
        if not math.isfinite(x) or not lo <= x <= hi:
            raise ValueError(f"{name}[{i}] debe estar en [{lo}, {hi}], recibido {x!r}")
    if not integer:
        object.__setattr__(obj, name, tuple(float(x) for x in v))


@dataclass(frozen=True, slots=True)
class DistanceSuiteConfig:
    modes: tuple[DistanceMode, ...] = (
        DistanceMode.HOP,
        DistanceMode.WEIGHTED_INVERSE,
        DistanceMode.WEIGHTED_LOG,
        DistanceMode.RESISTANCE,
    )
    metric_tol: float = 0.5
    log_floor: float = 1e-12
    fallback_n_radii: int = 12
    resistance_max_nodes: int = 2000
    zeta_1d_min: float = 0.8
    zeta_2d_min: float = 0.35
    zeta_2d_max: float = 0.65
    zeta_3d_max: float = 0.35

    def __post_init__(self) -> None:
        _enum_tuple(self, "modes", DistanceMode)
        _num(self, "metric_tol", 0.0, lo_open=True)
        _num(self, "log_floor", 0.0, 1.0, lo_open=True, hi_open=True)
        _int(self, "fallback_n_radii", 2)
        _int(self, "resistance_max_nodes", 2)
        for name in ("zeta_1d_min", "zeta_2d_min", "zeta_2d_max", "zeta_3d_max"):
            _num(self, name, 0.0, 1.0)
        if self.zeta_2d_min > self.zeta_2d_max:
            raise ValueError("zeta_2d_min debe ser <= zeta_2d_max")


@dataclass(frozen=True, slots=True)
class WeylConfig:
    laplacian: Literal["combinatorial", "normalized"] = "combinatorial"
    graph: Literal["binary", "thresholded_weighted"] = "binary"
    count_min: int = 10
    count_max_frac: float = 0.2
    min_levels: int = 4
    degeneracy_rtol: float = 1e-9
    r2_min: float = 0.95
    min_nodes: int = 10

    def __post_init__(self) -> None:
        _choice("laplacian", self.laplacian, ("combinatorial", "normalized"))
        _choice("graph", self.graph, ("binary", "thresholded_weighted"))
        _int(self, "count_min", 1)
        _num(self, "count_max_frac", 0.0, 1.0, lo_open=True)
        _int(self, "min_levels", 2)
        _num(self, "degeneracy_rtol", 0.0, lo_open=True)
        _num(self, "r2_min", 0.0, 1.0)
        _int(self, "min_nodes", 2)


@dataclass(frozen=True, slots=True)
class TopologyConfig:
    thetas: tuple[float, ...] = (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01)
    clique_max_dim: int = 2
    max_simplices: int = 300_000
    short_cycle_length: int = 4
    max_faces: int = 500_000
    b1_density_max: float = 0.03
    stable_fraction: float = 0.8

    def __post_init__(self) -> None:
        _tuple_of(self, "thetas", 0.0, 1.0)
        if any(b >= a for a, b in zip(self.thetas, self.thetas[1:])):
            raise ValueError("thetas debe ser estrictamente decreciente")
        _int(self, "clique_max_dim", 1)
        _int(self, "max_simplices", 1)
        if _int(self, "short_cycle_length", 3) > 8:
            raise ValueError("short_cycle_length debe ser <= 8")
        _int(self, "max_faces", 1)
        _num(self, "b1_density_max", 0.0)
        _num(self, "stable_fraction", 0.0, 1.0, lo_open=True)


@dataclass(frozen=True, slots=True)
class CurvatureConfig:
    idleness: float = 0.5
    distance_mode: DistanceMode = DistanceMode.HOP
    max_edges: int = 1000
    mean_min: float = -0.1
    tail_cut: float = -0.5
    tail_frac_max: float = 0.03
    qrc_deltas: tuple[int, ...] = (1, 2, 3)
    qrc_max_pairs: int = 200

    def __post_init__(self) -> None:
        _num(self, "idleness", 0.0, 1.0, hi_open=True)
        if not isinstance(self.distance_mode, DistanceMode):
            raise TypeError("distance_mode debe ser DistanceMode")
        _int(self, "max_edges", 1)
        _num(self, "mean_min", -2.0, 1.0)
        _num(self, "tail_cut", -2.0, 1.0)
        _num(self, "tail_frac_max", 0.0, 1.0)
        _tuple_of(self, "qrc_deltas", 1, 1_000_000, integer=True)
        _int(self, "qrc_max_pairs", 1)


@dataclass(frozen=True, slots=True)
class LocalStructureConfig:
    homogeneity_cv_max: float = 0.2
    homogeneity_min_nodes: int = 10
    isotropy_sources: int = 64
    isotropy_radius_offset: int = 1
    isotropy_median_min: float = 0.5
    isotropy_p10_min: float = 0.3
    locality_min: float = 0.97
    annulus_sources: int = 64
    annulus_min: float = 0.9
    annulus_r_min: int = 2
    small_world_clustering_ratio: float = 3.0

    def __post_init__(self) -> None:
        _num(self, "homogeneity_cv_max", 0.0, lo_open=True)
        _int(self, "homogeneity_min_nodes", 2)
        _int(self, "isotropy_sources", 1)
        _int(self, "isotropy_radius_offset", 0)
        _num(self, "isotropy_median_min", 0.0, 1.0)
        _num(self, "isotropy_p10_min", 0.0, 1.0)
        if self.isotropy_p10_min > self.isotropy_median_min:
            raise ValueError("isotropy_p10_min debe ser <= isotropy_median_min")
        _num(self, "locality_min", 0.0, 1.0)
        _int(self, "annulus_sources", 1)
        _num(self, "annulus_min", 0.0, 1.0)
        _int(self, "annulus_r_min", 1)
        _num(self, "small_world_clustering_ratio", 1.0)


@dataclass(frozen=True, slots=True)
class CertificateThresholds:
    g_connected: float = 0.95
    dim_tol: float = 0.35
    class_tol: float = 0.35
    seed_fraction: float = 0.8
    seed_std: float = 0.2
    size_range_tol: float = 0.3
    size_slope_tol: float = 0.1
    min_sizes: int = 3
    homogeneity_size_slack: float = 0.02
    null_pass_max: float = 0.2
    null_sigma: float = 3.0
    uniform_cv_max: float = 1e-6
    empty_g: float = 0.1
    empty_kbin: float = 1.0
    dense_rho: float = 0.5
    dense_meanw: float = 0.5

    def __post_init__(self) -> None:
        for name in ("g_connected", "seed_fraction", "null_pass_max", "empty_g", "dense_rho", "dense_meanw"):
            _num(self, name, 0.0, 1.0)
        for name in (
            "dim_tol", "class_tol", "seed_std", "size_range_tol", "size_slope_tol",
            "homogeneity_size_slack", "null_sigma", "uniform_cv_max", "empty_kbin",
        ):
            _num(self, name, 0.0)
        _int(self, "min_sizes", 2)


@dataclass(frozen=True, slots=True)
class CoarseGrainConfig:
    rule: Literal["heavy_edge_matching"] = "heavy_edge_matching"
    aggregation: Literal["max", "mean"] = "max"
    n_min: int = 50
    max_levels: int = 4
    replicates: int = 3

    def __post_init__(self) -> None:
        _choice("rule", self.rule, ("heavy_edge_matching",))
        _choice("aggregation", self.aggregation, ("max", "mean"))
        _int(self, "n_min", 2)
        _int(self, "max_levels", 1)
        _int(self, "replicates", 1)


@dataclass(frozen=True, slots=True)
class LangevinConfig:
    theta_hat: float
    n_steps: int = 20_000
    dt_safety: float = 0.1
    boundary: Literal["reflect", "project"] = "reflect"
    thin: int = 10
    burn_in_fraction: float = 0.5

    def __post_init__(self) -> None:
        _num(self, "theta_hat", 0.0, lo_open=True)
        _int(self, "n_steps", 1)
        _num(self, "dt_safety", 0.0, 1.0, lo_open=True, hi_open=True)
        _choice("boundary", self.boundary, ("reflect", "project"))
        _int(self, "thin", 1)
        _num(self, "burn_in_fraction", 0.0, 1.0, hi_open=True)


@dataclass(frozen=True, slots=True)
class MetropolisConfig:
    theta_hat: float
    n_sweeps: int = 2_000
    step_init: float = 0.2
    acceptance_low: float = 0.3
    acceptance_high: float = 0.5
    adapt_every: int = 10
    thin: int = 1
    burn_in_fraction: float = 0.5

    def __post_init__(self) -> None:
        _num(self, "theta_hat", 0.0, lo_open=True)
        _int(self, "n_sweeps", 1)
        _num(self, "step_init", 0.0, 1.0, lo_open=True)
        _num(self, "acceptance_low", 0.0, 1.0, lo_open=True)
        _num(self, "acceptance_high", 0.0, 1.0, hi_open=True)
        if self.acceptance_low >= self.acceptance_high:
            raise ValueError("acceptance_low debe ser < acceptance_high")
        _int(self, "adapt_every", 1)
        _int(self, "thin", 1)
        _num(self, "burn_in_fraction", 0.0, 1.0, hi_open=True)


@dataclass(frozen=True, slots=True)
class EnsembleConfig:
    n_chains: int = 4
    autocorr_c: float = 5.0
    rhat_max: float = 1.05
    ess_min: float = 100.0
    geweke_z_max: float = 3.0
    phase_samples: int = 5

    def __post_init__(self) -> None:
        _int(self, "n_chains", 2)
        _num(self, "autocorr_c", 0.0, lo_open=True)
        _num(self, "rhat_max", 1.0)
        _num(self, "ess_min", 0.0, lo_open=True)
        _num(self, "geweke_z_max", 0.0, lo_open=True)
        _int(self, "phase_samples", 1)


@dataclass(frozen=True, slots=True)
class FixedDensityConfig:
    rho: float | None = None
    projection_tol: float = 1e-9

    def __post_init__(self) -> None:
        if self.rho is not None:
            _num(self, "rho", 0.0, 1.0, lo_open=True)
        _num(self, "projection_tol", 0.0, lo_open=True)


@dataclass(frozen=True, slots=True)
class FiniteSizeConfig:
    sizes: tuple[int, ...] = (64, 100, 150, 200, 300, 500, 800)
    xl_sizes: tuple[int, ...] = (1500, 3000)
    replicates: int = 10

    def __post_init__(self) -> None:
        _tuple_of(self, "sizes", 3, 10_000_000, integer=True)
        _tuple_of(self, "xl_sizes", 3, 10_000_000, integer=True)
        _int(self, "replicates", 1)


@dataclass(frozen=True, slots=True)
class Omega11Config:
    base: OmegaConfig
    distance: DistanceSuiteConfig = field(default_factory=DistanceSuiteConfig)
    weyl: WeylConfig = field(default_factory=WeylConfig)
    topology: TopologyConfig = field(default_factory=TopologyConfig)
    curvature: CurvatureConfig = field(default_factory=CurvatureConfig)
    local: LocalStructureConfig = field(default_factory=LocalStructureConfig)
    certificate: CertificateThresholds = field(default_factory=CertificateThresholds)
    coarse: CoarseGrainConfig = field(default_factory=CoarseGrainConfig)
    ensemble: EnsembleConfig = field(default_factory=EnsembleConfig)
    finite_size: FiniteSizeConfig = field(default_factory=FiniteSizeConfig)
    engine: Engine = Engine.GRADIENT
    langevin: LangevinConfig | None = None
    metropolis: MetropolisConfig | None = None
    fixed_density: FixedDensityConfig | None = None
    null_models: tuple[NullModel, ...] = (
        NullModel.DEGREE_PRESERVING_REWIRE,
        NullModel.CONFIGURATION_MODEL,
        NullModel.ERDOS_RENYI,
        NullModel.SHUFFLED_WEIGHTS,
    )
    schema_version: str = "1.1"

    def __post_init__(self) -> None:
        expected: tuple[tuple[str, type], ...] = (
            ("base", OmegaConfig), ("distance", DistanceSuiteConfig), ("weyl", WeylConfig),
            ("topology", TopologyConfig), ("curvature", CurvatureConfig),
            ("local", LocalStructureConfig), ("certificate", CertificateThresholds),
            ("coarse", CoarseGrainConfig), ("ensemble", EnsembleConfig),
            ("finite_size", FiniteSizeConfig),
        )
        for name, cls in expected:
            if not isinstance(getattr(self, name), cls):
                raise TypeError(f"{name} debe ser {cls.__name__}")
        if self.base.schema_version != "1.0":
            raise ValueError("base.schema_version debe ser '1.0'")
        if not isinstance(self.engine, Engine):
            raise TypeError("engine debe ser Engine")
        for name, cls in (
            ("langevin", LangevinConfig), ("metropolis", MetropolisConfig), ("fixed_density", FixedDensityConfig)
        ):
            val = getattr(self, name)
            if val is not None and not isinstance(val, cls):
                raise TypeError(f"{name} debe ser {cls.__name__} o None")
        required = {
            Engine.LANGEVIN: "langevin",
            Engine.METROPOLIS: "metropolis",
            Engine.FIXED_DENSITY: "fixed_density",
        }.get(self.engine)
        if required is not None and getattr(self, required) is None:
            raise ValueError(f"engine={self.engine.name} exige {required}")
        if self.engine in (Engine.LANGEVIN, Engine.METROPOLIS) and self.base.functional.eta != 0.0:
            raise ValueError("Langevin/Metropolis exigen base.functional.eta == 0")
        _enum_tuple(self, "null_models", NullModel)
        if self.schema_version != "1.1":
            raise ValueError("schema_version debe ser '1.1'")
