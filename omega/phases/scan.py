"""Simulacion de una corrida y barrido del diagrama de fases (ANALYSIS §5.8, D-8, D-15, D-28).

Todo es explicito: la aleatoriedad viene de `SeedKey`/`make_rng`; no hay estado global.
"""

from __future__ import annotations

import dataclasses
import itertools
import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from omega.config.convert import raw_to_reduced, reduced_to_raw
from omega.config.seeds import SeedKey, make_rng, seed_key
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
from omega.dynamics.evolution import evolve
from omega.geometry.observables import geometry_observables
from omega.network.initialization import random_uniform_weights
from omega.network.topology import topology_observables
from omega.network.weights import upper_triangle, validate_weight_matrix
from omega.phases.classification import classify_run
from omega.statistics.summary import proportion_ci
from omega.types import (
    FloatArray,
    Observables,
    PhaseLabel,
    RunResult,
    RunStatus,
    Trajectory,
)

__all__ = [
    "default_config",
    "with_experiment",
    "with_params",
    "parameter_grid",
    "point_params",
    "observe",
    "simulate_from",
    "simulate",
    "static_run",
    "scan_phase_diagram",
    "phase_probabilities",
    "run_summary",
    "to_jsonable",
]


def default_config(
    n: int = 200, *, master_entropy: int = 20240901, experiment_id: int = 0, replicates: int = 10
) -> OmegaConfig:
    """Configuracion preregistrada por defecto (umbrales y D-* de ANALYSIS §4 sin tocar)."""
    return OmegaConfig(
        init=InitConfig(n=n),
        functional=FunctionalParams(alpha=0.0),
        dynamics=DynamicsConfig(),
        graph=GraphConfig(),
        dimension=DimensionConfig(),
        spectral=SpectralConfig(),
        phases=PhaseThresholds(),
        seeds=SeedConfig(master_entropy=master_entropy, experiment_id=experiment_id, replicates=replicates),
    )


def with_experiment(cfg: OmegaConfig, experiment_id: int) -> OmegaConfig:
    """Copia de cfg con otro `seeds.experiment_id` (D-28)."""
    return dataclasses.replace(cfg, seeds=dataclasses.replace(cfg.seeds, experiment_id=experiment_id))


def with_params(cfg: OmegaConfig, p: FunctionalParams, n: int | None = None) -> OmegaConfig:
    """Copia de cfg con funcional `p` y, opcionalmente, otro numero de nodos."""
    init = cfg.init if n is None else dataclasses.replace(cfg.init, n=n)
    return dataclasses.replace(cfg, functional=p, init=init)


def parameter_grid(scan: ScanConfig) -> tuple[tuple[int, float, float], ...]:
    """Producto cartesiano (alpha_hat x gamma_hat) con indice de punto, alpha-mayor (D-8)."""
    if not isinstance(scan, ScanConfig):
        raise TypeError("scan debe ser ScanConfig")
    return tuple((i, a, g) for i, (a, g) in enumerate(itertools.product(scan.alpha_hat, scan.gamma_hat)))


def point_params(
    alpha_hat: float, gamma_hat: float, n: int, scan: ScanConfig | None = None, beta: float = 1.0
) -> FunctionalParams:
    """Parametros crudos de un punto de la malla.

    scaling 'reduced' (o scan None): (alpha_hat, gamma_hat) reducidos -> reduced_to_raw (D-15).
    scaling 'raw': los valores son directamente alpha/beta y gamma/beta (malla de M§25).
    """
    if scan is not None and scan.scaling == "raw":
        return FunctionalParams(alpha=alpha_hat * beta, beta=beta, gamma=gamma_hat * beta)
    return reduced_to_raw(alpha_hat, gamma_hat, n, beta)


def observe(w: FloatArray, cfg: OmegaConfig) -> Observables:
    """Observables topologicos y geometricos de un estado W (D-2, D-11, D-22)."""
    validate_weight_matrix(w)
    top = topology_observables(w, cfg.graph.w_min, cfg.phases.large_component_frac)
    geo = geometry_observables(w, cfg.graph, cfg.dimension, cfg.spectral)
    return Observables(topology=top, geometry=geo)


def simulate_from(cfg: OmegaConfig, p: FunctionalParams, key: SeedKey, w0: FloatArray) -> RunResult:
    """evolve -> observables -> classify desde un W0 dado (continuacion/histeresis)."""
    validate_weight_matrix(w0)
    if w0.shape[0] != cfg.init.n:
        raise ValueError(f"w0 tiene {w0.shape[0]} nodos, cfg.init.n={cfg.init.n}")
    traj = evolve(w0, p, cfg.dynamics)
    obs = observe(traj.w_final, cfg)
    a_hat, g_hat = raw_to_reduced(p, cfg.init.n)
    return RunResult(
        seed=key,
        params=p,
        alpha_hat=a_hat,
        gamma_hat=g_hat,
        w0=np.array(w0, dtype=np.float64, copy=True),
        trajectory=traj,
        observables=obs,
        assessment=classify_run(obs, traj.status, cfg.phases),
    )


def simulate(cfg: OmegaConfig, p: FunctionalParams, key: SeedKey) -> RunResult:
    """init (U(0,1), rng=make_rng(key)) -> evolve -> observables -> classify."""
    w0 = random_uniform_weights(cfg.init.n, make_rng(key))
    return simulate_from(cfg, p, key, w0)


def static_run(cfg: OmegaConfig, p: FunctionalParams, key: SeedKey, w: FloatArray) -> RunResult:
    """RunResult de un estado SIN dinamica (linea base, Exp 1): steps=0, dt=tau=0, estado CONVERGED
    (estado estatico), instantanea unica; `p` es solo un marcador de parametros."""
    validate_weight_matrix(w)
    if w.shape[0] != cfg.init.n:
        raise ValueError(f"w tiene {w.shape[0]} nodos, cfg.init.n={cfg.init.n}")
    edges = upper_triangle(w)
    traj = Trajectory(
        w_final=np.array(w, dtype=np.float64, copy=True),
        status=RunStatus.CONVERGED,
        steps=0,
        dt=0.0,
        tau=0.0,
        scalars={},
        snapshot_steps=np.zeros(1, dtype=np.int64),
        snapshots=edges.reshape(1, -1),
    )
    obs = observe(w, cfg)
    a_hat, g_hat = raw_to_reduced(p, cfg.init.n)
    return RunResult(
        seed=key,
        params=p,
        alpha_hat=a_hat,
        gamma_hat=g_hat,
        w0=np.array(w, dtype=np.float64, copy=True),
        trajectory=traj,
        observables=obs,
        assessment=classify_run(obs, traj.status, cfg.phases),
    )


def scan_phase_diagram(cfg: OmegaConfig) -> tuple[RunResult, ...]:
    """Barrido completo: para cada punto de cfg.scan y cada replica, simulate con
    seed_key(cfg.seeds, punto, replica). Independiente del orden de ejecucion (D-28)."""
    if cfg.scan is None:
        raise ValueError("cfg.scan es None: no hay malla que barrer")
    out: list[RunResult] = []
    for idx, a_hat, g_hat in parameter_grid(cfg.scan):
        p = point_params(a_hat, g_hat, cfg.init.n, cfg.scan)
        for rep in range(cfg.seeds.replicates):
            out.append(simulate(cfg, p, seed_key(cfg.seeds, idx, rep)))
    return tuple(out)


def phase_probabilities(
    results: Sequence[RunResult], level: float = 0.95
) -> dict[tuple[float, float], dict[str, Any]]:
    """Por punto (alpha_hat, gamma_hat) redondeado a 9 decimales: n, conteos, probabilidades e
    IC de Wilson por etiqueta (todas las etiquetas aparecen, tambien con 0)."""
    groups: dict[tuple[float, float], list[RunResult]] = {}
    for r in results:
        groups.setdefault((round(r.alpha_hat, 9), round(r.gamma_hat, 9)), []).append(r)
    out: dict[tuple[float, float], dict[str, Any]] = {}
    for pt, rs in groups.items():
        n = len(rs)
        counts = {lab.value: 0 for lab in PhaseLabel}
        for r in rs:
            counts[r.assessment.label.value] += 1
        out[pt] = {
            "n": n,
            "counts": counts,
            "probabilities": {k: v / n for k, v in counts.items()},
            "ci": {k: proportion_ci(v, n, level) for k, v in counts.items()},
        }
    return out


def to_jsonable(obj: Any) -> Any:
    """Convierte a tipos JSON puros: numpy -> python, no finitos -> None, claves -> str,
    enums -> valor."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return to_jsonable(obj.item())
    if isinstance(obj, (PhaseLabel, RunStatus)):
        return obj.value
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    return obj


def run_summary(r: RunResult) -> dict[str, Any]:
    """Fila resumen (JSON-compatible) de una corrida."""
    top, geo = r.observables.topology, r.observables.geometry
    out = to_jsonable(
        {
            "spawn_key": list(r.seed.spawn_key),
            "alpha_hat": r.alpha_hat,
            "gamma_hat": r.gamma_hat,
            "label": r.assessment.label,
            "flags": {k: bool(v) for k, v in r.assessment.flags.items()},
            "status": r.trajectory.status,
            "steps": r.trajectory.steps,
            "mean_strength": top.mean_strength,
            "cv_strength": top.cv_strength,
            "mean_binary_degree": top.mean_binary_degree,
            "binary_density": top.binary_density,
            "mean_weight": top.mean_weight,
            "giant_fraction": top.giant_fraction,
            "n_large_components": top.n_large_components,
            "clustering_binary": top.clustering_binary,
            "clustering_weighted": top.clustering_weighted,
            "frac_at_zero": top.frac_at_zero,
            "frac_at_one": top.frac_at_one,
            "d_eff": geo.d_eff.value,
            "d_eff_status": geo.d_eff.status,
            "d_eff_plateau": geo.d_eff.plateau,
            "d_s": geo.d_s.value,
            "d_s_status": geo.d_s.status,
        }
    )
    assert isinstance(out, dict)
    return out
