"""Descenso de gradiente proyectado/sigmoide sobre S0 (M§16-M§17; ANALYSIS D-3, D-4, D-13, R5).

No calcula geometria: solo escalares baratos por paso e instantaneas.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.special import expit

from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.dynamics.functional import action, triangles
from omega.dynamics.gradient import grad_action, lipschitz_bound, theta_gradient
from omega.network.weights import from_upper_triangle, upper_triangle, validate_weight_matrix
from omega.types import FloatArray, IntArray, RunStatus, Trajectory

__all__ = [
    "step_size",
    "projected_step",
    "sigmoid",
    "logit",
    "sigmoid_step",
    "snapshot_steps",
    "step_scalars",
    "is_converged",
    "evolve",
]


def step_size(n: int, p: FunctionalParams, cfg: DynamicsConfig) -> float:
    """dt = dt_fixed (modo fixed) o dt_safety/L (modo auto; exige eta=0) (D-3)."""
    if not isinstance(cfg, DynamicsConfig):
        raise TypeError("cfg debe ser DynamicsConfig")
    if cfg.dt_mode == "fixed":
        assert cfg.dt_fixed is not None
        return float(cfg.dt_fixed)
    return float(cfg.dt_safety / lipschitz_bound(n, p))


def projected_step(w: FloatArray, g: FloatArray, dt: float) -> FloatArray:
    """W <- clip(W - dt G, 0, 1) con diagonal 0 (M§17, §1.5); no muta entradas."""
    validate_weight_matrix(w)
    g = np.asarray(g, dtype=np.float64)
    if g.shape != w.shape:
        raise ValueError(f"g y w deben tener la misma forma: {g.shape} vs {w.shape}")
    if not dt > 0.0:
        raise ValueError(f"dt debe ser > 0, recibido {dt}")
    out = np.clip(w - dt * g, 0.0, 1.0)
    np.fill_diagonal(out, 0.0)
    return np.asarray(out, dtype=np.float64)


def sigmoid(theta: FloatArray) -> FloatArray:
    """Logistica 1/(1+exp(-theta)) estable (D-13)."""
    return np.asarray(expit(np.asarray(theta, dtype=np.float64)), dtype=np.float64)


def logit(w: FloatArray, delta: float) -> FloatArray:
    """theta = log(x/(1-x)) con x = clip(w, delta, 1-delta) (D-13)."""
    if not 0.0 < delta < 0.5:
        raise ValueError(f"delta debe estar en (0,0.5), recibido {delta}")
    x = np.clip(np.asarray(w, dtype=np.float64), delta, 1.0 - delta)
    return np.asarray(np.log(x / (1.0 - x)), dtype=np.float64)


def _weights_from_theta(theta: FloatArray) -> FloatArray:
    w = sigmoid(theta)
    np.fill_diagonal(w, 0.0)
    return w


def sigmoid_step(theta: FloatArray, p: FunctionalParams, dt: float) -> FloatArray:
    """theta <- theta - dt * G W (1-W), con W = sigmoide(theta) y diag(W)=0 (D-13)."""
    theta = np.asarray(theta, dtype=np.float64)
    if theta.ndim != 2 or theta.shape[0] != theta.shape[1]:
        raise ValueError(f"theta debe ser cuadrada, forma {theta.shape}")
    if not dt > 0.0:
        raise ValueError(f"dt debe ser > 0, recibido {dt}")
    w = _weights_from_theta(theta)
    g = grad_action(w, p)
    out = theta - dt * theta_gradient(g, w)
    np.fill_diagonal(out, 0.0)
    return np.asarray(out, dtype=np.float64)


def snapshot_steps(max_steps: int, schedule: str) -> IntArray:
    """Pasos {0,1,2,4,...} <= max_steps para 'log2'; vacio para 'none' (D-4)."""
    if isinstance(max_steps, bool) or not isinstance(max_steps, int) or max_steps < 1:
        raise ValueError(f"max_steps debe ser entero >= 1, recibido {max_steps!r}")
    if schedule == "none":
        return np.zeros(0, dtype=np.int64)
    if schedule != "log2":
        raise ValueError(f"schedule invalido: {schedule!r}")
    steps = [0]
    s = 1
    while s <= max_steps:
        steps.append(s)
        s *= 2
    return np.asarray(steps, dtype=np.int64)


def _projected_inf_norm(w: FloatArray, g: FloatArray) -> float:
    """||proy_[0,1](G)||_inf: anula G donde empuja fuera de la caja (R5)."""
    blocked = ((w <= 0.0) & (g > 0.0)) | ((w >= 1.0) & (g < 0.0))
    return float(np.max(np.abs(np.where(blocked, 0.0, g))))


def step_scalars(
    w: FloatArray, w_prev: FloatArray, p: FunctionalParams, g: FloatArray | None = None
) -> dict[str, float]:
    """Escalares baratos de un paso (D-4, R5); `g` opcional evita recalcular el gradiente."""
    validate_weight_matrix(w)
    validate_weight_matrix(w_prev)
    if w.shape != w_prev.shape:
        raise ValueError("w y w_prev deben tener la misma forma")
    grad = grad_action(w, p) if g is None else g
    n = w.shape[0]
    edges = upper_triangle(w)
    return {
        "action": action(w, p),
        "triangles": triangles(w),
        "max_dw": float(np.max(np.abs(w - w_prev))),
        "grad_proj_inf": _projected_inf_norm(w, grad),
        "frac_zero": float(np.mean(edges == 0.0)),
        "frac_one": float(np.mean(edges == 1.0)),
        "mean_strength": float(w.sum() / n),
        "mean_weight": float(edges.mean()),
    }


def is_converged(max_dw_recent: Sequence[float] | FloatArray, tol: float, patience: int) -> bool:
    """True si hay >= patience valores y los ultimos `patience` son todos < tol (D-4)."""
    if isinstance(patience, bool) or not isinstance(patience, int) or patience < 1:
        raise ValueError(f"patience debe ser entero >= 1, recibido {patience!r}")
    if not tol > 0.0:
        raise ValueError(f"tol debe ser > 0, recibido {tol}")
    arr = np.asarray(max_dw_recent, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError("max_dw_recent debe ser unidimensional")
    if arr.shape[0] < patience:
        return False
    return bool(np.all(arr[-patience:] < tol))


def evolve(w0: FloatArray, p: FunctionalParams, cfg: DynamicsConfig) -> Trajectory:
    """Evoluciona W0 hasta convergencia, max_steps o valores no finitos (D-4, M§17).

    Escalares con longitud steps+1 (indice 0 = estado inicial, max_dw=0). Si hay NONFINITE
    se descarta el paso invalido: w_final es el ultimo estado finito.
    """
    validate_weight_matrix(w0)
    if not isinstance(p, FunctionalParams):
        raise TypeError("p debe ser FunctionalParams")
    if not isinstance(cfg, DynamicsConfig):
        raise TypeError("cfg debe ser DynamicsConfig")
    n = w0.shape[0]
    dt = step_size(n, p, cfg)
    use_sigmoid = cfg.integrator == "sigmoid"
    sched = {int(s) for s in snapshot_steps(cfg.max_steps, cfg.snapshot_schedule)}

    w = np.array(w0, dtype=np.float64, copy=True)
    theta = logit(w, cfg.sigmoid_clip) if use_sigmoid else w
    if use_sigmoid:
        np.fill_diagonal(theta, 0.0)
        w = _weights_from_theta(theta)

    records: dict[str, list[float]] = {}
    snaps: list[FloatArray] = []
    snap_idx: list[int] = []

    def record(w_new: FloatArray, w_old: FloatArray, g: FloatArray) -> None:
        for key, val in step_scalars(w_new, w_old, p, g).items():
            records.setdefault(key, []).append(val)

    g = grad_action(w, p)
    record(w, w, g)
    if 0 in sched:
        snaps.append(upper_triangle(w))
        snap_idx.append(0)

    status = RunStatus.MAX_STEPS
    steps = 0
    for t in range(1, cfg.max_steps + 1):
        if use_sigmoid:
            theta_new = theta - dt * theta_gradient(g, w)
            np.fill_diagonal(theta_new, 0.0)
            w_new = _weights_from_theta(theta_new)
        else:
            theta_new = theta
            w_new = np.clip(w - dt * g, 0.0, 1.0)
            np.fill_diagonal(w_new, 0.0)
        with np.errstate(all="ignore"):
            finite = bool(np.all(np.isfinite(w_new)) and np.all(np.isfinite(theta_new)))
            g_new = grad_action(w_new, p) if finite else g
            finite = finite and bool(np.all(np.isfinite(g_new)))
        if not finite:
            status = RunStatus.NONFINITE
            break
        record(w_new, w, g_new)
        w, theta, g = w_new, theta_new, g_new
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

    m = n * (n - 1) // 2
    snapshots = np.vstack(snaps).astype(np.float64) if snaps else np.zeros((0, m))
    return Trajectory(
        w_final=from_upper_triangle(upper_triangle(w), n),
        status=status,
        steps=steps,
        dt=dt,
        tau=steps * dt,
        scalars={k: np.asarray(v, dtype=np.float64) for k, v in records.items()},
        snapshot_steps=np.asarray(snap_idx, dtype=np.int64),
        snapshots=snapshots,
    )
