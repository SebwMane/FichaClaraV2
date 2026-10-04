"""Tipos compartidos del simulador (ANALYSIS §5.4). Dataclasses frozen sin logica."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from typing import Literal, TypeAlias

import numpy as np
import numpy.typing as npt

from omega.config.seeds import SeedKey
from omega.config.settings import FunctionalParams

FloatArray: TypeAlias = npt.NDArray[np.float64]
BoolArray: TypeAlias = npt.NDArray[np.bool_]
IntArray: TypeAlias = npt.NDArray[np.int64]

EstimateStatus: TypeAlias = Literal["ok", "no_window", "insufficient_component"]


class RunStatus(Enum):
    """Estado final de una evolucion (D-4)."""

    CONVERGED = "CONVERGED"
    MAX_STEPS = "MAX_STEPS"
    NONFINITE = "NONFINITE"


class PhaseLabel(Enum):
    """Fases preregistradas (ANALYSIS §4.1); U = no clasificada."""

    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"
    F = "F"
    U = "U"


@dataclass(frozen=True, slots=True)
class DimensionEstimate:
    """Estimacion de dimension (D_eff o D_s) con ventana y estado (D-6, D-7, D-22)."""

    value: float
    stderr: float
    status: EstimateStatus
    window: tuple[int, int] | None
    scales: FloatArray
    profile: FloatArray
    local_slopes: FloatArray
    plateau: bool
    method: str


@dataclass(frozen=True, slots=True)
class TopologyObservables:
    """Observables topologicos de una red de pesos (D-10, D-11)."""

    n: int
    mean_strength: float
    std_strength: float
    cv_strength: float
    mean_binary_degree: float
    binary_density: float
    mean_weight: float
    giant_size: int
    giant_fraction: float
    n_components: int
    n_large_components: int
    clustering_weighted: float
    clustering_binary: float
    frac_at_zero: float
    frac_at_one: float


@dataclass(frozen=True, slots=True)
class GeometryObservables:
    """Observables geometricos sobre la componente gigante (D-11, D-22)."""

    path_length: float
    path_length_hops: float
    diameter: float
    d_eff: DimensionEstimate
    d_eff_ball: DimensionEstimate
    d_eff_hops: DimensionEstimate
    d_s: DimensionEstimate


@dataclass(frozen=True, slots=True)
class Observables:
    """Conjunto completo de observables de un estado."""

    topology: TopologyObservables
    geometry: GeometryObservables


@dataclass(frozen=True, slots=True)
class Trajectory:
    """Registro de una evolucion. Las instantaneas son (n_snap, M), triangular superior."""

    w_final: FloatArray
    status: RunStatus
    steps: int
    dt: float
    tau: float
    scalars: Mapping[str, FloatArray]
    snapshot_steps: IntArray
    snapshots: FloatArray


@dataclass(frozen=True, slots=True)
class PhaseAssessment:
    """Etiqueta de fase y todas las banderas que la sustentan."""

    label: PhaseLabel
    flags: Mapping[str, bool]


@dataclass(frozen=True, slots=True)
class RunResult:
    """Resultado completo de una corrida."""

    seed: SeedKey
    params: FunctionalParams
    alpha_hat: float
    gamma_hat: float
    w0: FloatArray
    trajectory: Trajectory
    observables: Observables
    assessment: PhaseAssessment
