"""Ablacion del funcional S0 (P§10; diseno §3.1): quita terminos y evoluciona con clip.

Funciones puras; reutiliza funcional y gradiente de Omega-1.0.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

import numpy as np

from omega.config.convert import reduced_to_raw
from omega.config.settings import DynamicsConfig
from omega.dynamics.evolution import is_converged, snapshot_steps
from omega.dynamics.functional import degree_irregularity, density, triangles
from omega.dynamics.gradient import grad_degree_irregularity, grad_density, grad_triangles
from omega.network.weights import from_upper_triangle, upper_triangle, validate_weight_matrix
from omega.types import FloatArray, RunStatus, Trajectory

__all__ = [
    "AblationSpec",
    "AblatedParams",
    "effective_coefficients",
    "ablated_action",
    "ablated_gradient",
    "ablated_lipschitz",
    "evolve_ablated",
    "ablation_set",
]


class AblationSpec(Enum):
    """Variante del funcional: terminos que se conservan (P§10)."""

    FULL = "full"
    NO_TRIANGLES = "no_triangles"
    NO_DENSITY = "no_density"
    NO_DEGREE = "no_degree"
    TRIANGLES_ONLY = "triangles_only"


@dataclass(frozen=True, slots=True)
class AblatedParams:
    """Coeficientes crudos (alpha, beta, gamma) y variante; beta >= 0 (FunctionalParams exige > 0)."""

    alpha: float
    beta: float
    gamma: float
    spec: AblationSpec

    def __post_init__(self) -> None:
        for name in ("alpha", "beta", "gamma"):
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise TypeError(f"{name} debe ser real")
            if not math.isfinite(v):
                raise ValueError(f"{name} debe ser finito")
        if self.beta < 0.0:
            raise ValueError(f"beta debe ser >= 0, recibido {self.beta}")
        if self.gamma < 0.0:
            raise ValueError(f"gamma debe ser >= 0, recibido {self.gamma}")
        if not isinstance(self.spec, AblationSpec):
            raise TypeError("spec debe ser AblationSpec")


def effective_coefficients(ap: AblatedParams) -> tuple[float, float, float]:
    """(alpha_e, beta_e, gamma_e): pone a 0 los terminos ablacionados."""
    if not isinstance(ap, AblatedParams):
        raise TypeError("ap debe ser AblatedParams")
    a, b, g = float(ap.alpha), float(ap.beta), float(ap.gamma)
    if ap.spec is AblationSpec.FULL:
        return a, b, g
    if ap.spec is AblationSpec.NO_TRIANGLES:
        return 0.0, b, g
    if ap.spec is AblationSpec.NO_DENSITY:
        return a, 0.0, g
    if ap.spec is AblationSpec.NO_DEGREE:
        return a, b, 0.0
    return a, 0.0, 0.0  # TRIANGLES_ONLY


def ablated_action(w: FloatArray, ap: AblatedParams) -> float:
    """S = -alpha_e T + beta_e S_dens + gamma_e S_deg (M§13)."""
    validate_weight_matrix(w)
    a, b, g = effective_coefficients(ap)
    total = 0.0
    if a != 0.0:
        total -= a * triangles(w)
    if b != 0.0:
        total += b * density(w)
    if g != 0.0:
        total += g * degree_irregularity(w)
    return float(total)


def ablated_gradient(w: FloatArray, ap: AblatedParams) -> FloatArray:
    """G = -alpha_e G_T + beta_e G_dens + gamma_e G_deg (matriz simetrica, diagonal 0)."""
    validate_weight_matrix(w)
    a, b, g = effective_coefficients(ap)
    out = np.zeros_like(w)
    if a != 0.0:
        out = out - a * grad_triangles(w)
    if b != 0.0:
        out = out + b * grad_density(w)
    if g != 0.0:
        out = out + g * grad_degree_irregularity(w)
    np.fill_diagonal(out, 0.0)
    return np.asarray(out, dtype=np.float64)


def ablated_lipschitz(n: int, ap: AblatedParams) -> float:
    """L = 2 beta_e + 2 (gamma_e + |alpha_e|)(N-2); ValueError si L = 0 (D-3)."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError(f"n debe ser entero >= 2, recibido {n!r}")
    a, b, g = effective_coefficients(ap)
    lip = 2.0 * b + 2.0 * (g + abs(a)) * (n - 2)
    if lip <= 0.0:
        raise ValueError("cota de Lipschitz nula: use dt_fixed")
    return float(lip)


def _scalars(w: FloatArray, w_prev: FloatArray, g: FloatArray, ap: AblatedParams) -> dict[str, float]:
    """Mismas claves que `step_scalars` de Omega-1.0, con la accion ablacionada."""
    n = w.shape[0]
    edges = upper_triangle(w)
    blocked = ((w <= 0.0) & (g > 0.0)) | ((w >= 1.0) & (g < 0.0))
    return {
        "action": ablated_action(w, ap),
        "triangles": triangles(w),
        "max_dw": float(np.max(np.abs(w - w_prev))),
        "grad_proj_inf": float(np.max(np.abs(np.where(blocked, 0.0, g)))),
        "frac_zero": float(np.mean(edges == 0.0)),
        "frac_one": float(np.mean(edges == 1.0)),
        "mean_strength": float(w.sum() / n),
        "mean_weight": float(edges.mean()),
    }


def evolve_ablated(w0: FloatArray, ap: AblatedParams, cfg: DynamicsConfig) -> Trajectory:
    """Bucle clip identico a `evolve` de Omega-1.0 con el funcional ablacionado (FULL == evolve)."""
    validate_weight_matrix(w0)
    if not isinstance(ap, AblatedParams):
        raise TypeError("ap debe ser AblatedParams")
    if not isinstance(cfg, DynamicsConfig):
        raise TypeError("cfg debe ser DynamicsConfig")
    if cfg.integrator != "clip":
        raise ValueError("la ablacion solo admite integrator='clip'")
    n = w0.shape[0]
    if cfg.dt_mode == "fixed":
        assert cfg.dt_fixed is not None
        dt = float(cfg.dt_fixed)
    else:
        dt = float(cfg.dt_safety / ablated_lipschitz(n, ap))
    sched = {int(s) for s in snapshot_steps(cfg.max_steps, cfg.snapshot_schedule)}

    w = np.array(w0, dtype=np.float64, copy=True)
    g = ablated_gradient(w, ap)
    records: dict[str, list[float]] = {}

    def record(w_new: FloatArray, w_old: FloatArray, g_new: FloatArray) -> None:
        for key, val in _scalars(w_new, w_old, g_new, ap).items():
            records.setdefault(key, []).append(val)

    record(w, w, g)
    snaps: list[FloatArray] = []
    snap_idx: list[int] = []
    if 0 in sched:
        snaps.append(upper_triangle(w))
        snap_idx.append(0)

    status = RunStatus.MAX_STEPS
    steps = 0
    for t in range(1, cfg.max_steps + 1):
        with np.errstate(all="ignore"):
            w_new = np.clip(w - dt * g, 0.0, 1.0)
            np.fill_diagonal(w_new, 0.0)
            finite = bool(np.all(np.isfinite(w_new)))
            g_new = ablated_gradient(w_new, ap) if finite else g
            finite = finite and bool(np.all(np.isfinite(g_new)))
        if not finite:
            status = RunStatus.NONFINITE
            break
        record(w_new, w, g_new)
        w, g = w_new, g_new
        steps = t
        if t in sched:
            snaps.append(upper_triangle(w))
            snap_idx.append(t)
        if is_converged(records["max_dw"][1:][-cfg.patience :], cfg.tol_step, cfg.patience):
            status = RunStatus.CONVERGED
            break

    if not snap_idx or snap_idx[-1] != steps:
        snaps.append(upper_triangle(w))
        snap_idx.append(steps)
    return Trajectory(
        w_final=from_upper_triangle(upper_triangle(w), n),
        status=status,
        steps=steps,
        dt=dt,
        tau=steps * dt,
        scalars={k: np.asarray(v, dtype=np.float64) for k, v in records.items()},
        snapshot_steps=np.asarray(snap_idx, dtype=np.int64),
        snapshots=np.vstack(snaps).astype(np.float64),
    )


def ablation_set(
    alpha_hat: float, gamma_hat: float, n: int, beta: float = 1.0
) -> tuple[AblatedParams, ...]:
    """Las cinco variantes con parametros crudos de `reduced_to_raw` (D-15), en orden del Enum."""
    p = reduced_to_raw(alpha_hat, gamma_hat, n, beta)
    return tuple(AblatedParams(p.alpha, p.beta, p.gamma, spec) for spec in AblationSpec)
