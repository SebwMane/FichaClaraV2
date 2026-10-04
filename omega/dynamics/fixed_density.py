"""Rama Omega-B: dinamica a densidad fija sobre {x in [0,1]^M : sum x = W0} (P§9; diseno §1.9).

Solo integrador clip. Todas las funciones son puras; no hay estado global ni RNG.
"""

from __future__ import annotations

import numpy as np

from omega.config.settings import DynamicsConfig, FunctionalParams
from omega.config.settings11 import FixedDensityConfig
from omega.dynamics.evolution import is_converged, snapshot_steps, step_scalars, step_size
from omega.dynamics.gradient import grad_action
from omega.network.weights import from_upper_triangle, upper_triangle, validate_weight_matrix
from omega.types import FloatArray, RunStatus, Trajectory

__all__ = [
    "density_to_total",
    "project_capped_simplex",
    "project_fixed_density",
    "fixed_density_step",
    "evolve_fixed_density",
    "uniform_state_threshold",
    "uniform_state_unstable",
]


def _check_n(n: int, minimum: int) -> None:
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError("n debe ser entero")
    if n < minimum:
        raise ValueError(f"n debe ser >= {minimum}, recibido {n}")


def _check_real(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.floating, np.integer)):
        raise TypeError(f"{name} debe ser real")
    v = float(value)
    if not np.isfinite(v):
        raise ValueError(f"{name} debe ser finito")
    return v


def density_to_total(rho: float, n: int) -> float:
    """W0 = rho N(N-1)/2 (diseno §1.9); rho en [0,1], n >= 2."""
    r = _check_real("rho", rho)
    _check_n(n, 2)
    if not 0.0 <= r <= 1.0:
        raise ValueError(f"rho debe estar en [0,1], recibido {r}")
    return r * n * (n - 1) / 2.0


def _project_with_shift(v: FloatArray, total: float) -> tuple[FloatArray, float]:
    """x = clip(v - nu, 0, 1) con sum x = total; devuelve (x, nu). Puntos de quiebre {v, v-1}."""
    v = np.asarray(v, dtype=np.float64)
    if v.ndim != 1 or v.shape[0] < 1:
        raise ValueError("v debe ser un vector no vacio")
    if not np.all(np.isfinite(v)):
        raise ValueError("v contiene valores no finitos")
    t = _check_real("total", total)
    m = v.shape[0]
    if t < 0.0 or t > float(m):
        raise ValueError(f"total debe estar en [0, {m}], recibido {t}")
    s = np.sort(v)
    pref = np.concatenate(([0.0], np.cumsum(s)))
    bp = np.sort(np.concatenate((s, s - 1.0)))

    def f(nu: FloatArray) -> FloatArray:
        lo = np.searchsorted(s, nu, side="right")  # v <= nu
        hi = np.searchsorted(s, nu + 1.0, side="right")  # v <= nu + 1
        mid = pref[hi] - pref[lo] - nu * (hi - lo)
        return np.asarray((m - hi) + mid, dtype=np.float64)

    fb = f(bp)
    k = int(np.argmax(fb <= t))  # primer indice con f <= total (f no creciente)
    if fb[k] > t:  # no ocurre salvo redondeo extremo: f(max v) = 0 <= total
        k = bp.shape[0] - 1
    if k == 0:
        nu = float(bp[0])
    else:
        x0, x1, f0, f1 = float(bp[k - 1]), float(bp[k]), float(fb[k - 1]), float(fb[k])
        nu = x1 if f0 == f1 else x0 + (f0 - t) * (x1 - x0) / (f0 - f1)
    return np.asarray(np.clip(v - nu, 0.0, 1.0), dtype=np.float64), nu


def project_capped_simplex(v: FloatArray, total: float) -> FloatArray:
    """Proyeccion euclidea exacta sobre {x in [0,1]^M : sum x = total} (diseno §1.9).

    Precision |sum x - total| <= 1e-9 max(1, total); ValueError si total no esta en [0, M].
    """
    return _project_with_shift(v, total)[0]


def project_fixed_density(w: FloatArray, total: float) -> FloatArray:
    """Proyecta la triangular superior de W sobre densidad fija y reconstruye la matriz simetrica."""
    validate_weight_matrix(w)
    x = project_capped_simplex(upper_triangle(w), total)
    return from_upper_triangle(x, w.shape[0])


def _step_with_shift(w: FloatArray, g: FloatArray, dt: float, total: float) -> tuple[FloatArray, float]:
    validate_weight_matrix(w)
    g = np.asarray(g, dtype=np.float64)
    if g.shape != w.shape:
        raise ValueError(f"g y w deben tener la misma forma: {g.shape} vs {w.shape}")
    if not dt > 0.0:
        raise ValueError(f"dt debe ser > 0, recibido {dt}")
    n = w.shape[0]
    iu = np.triu_indices(n, k=1)
    x, nu = _project_with_shift(w[iu] - dt * g[iu], total)
    return from_upper_triangle(x, n), nu


def fixed_density_step(w: FloatArray, g: FloatArray, dt: float, total: float) -> FloatArray:
    """W <- Proy_C(W - dt G); el gradiente se extrae con triu_indices (diseno §1.9)."""
    return _step_with_shift(w, g, dt, total)[0]


def _multiplier_estimate(w: FloatArray, g: FloatArray) -> float:
    """Multiplicador KKT en el estado inicial: media de G sobre aristas libres (o todas)."""
    iu = np.triu_indices(w.shape[0], k=1)
    ge, we = g[iu], w[iu]
    free = (we > 0.0) & (we < 1.0)
    return float(np.mean(ge[free]) if np.any(free) else np.mean(ge))


def evolve_fixed_density(
    w0: FloatArray, p: FunctionalParams, cfg: DynamicsConfig, fd: FixedDensityConfig
) -> Trajectory:
    """Evoluciona con densidad fija (P§9). Proyecta W0; dt de Omega-1.0; paradas de Omega-1.0.

    Escalares de `step_scalars` mas "multiplier" = nu/dt (nu = desplazamiento de la proyeccion;
    en el indice 0, estimacion KKT). Solo integrador clip. Si fd.rho es None, rho = media de W0.
    """
    validate_weight_matrix(w0)
    if not isinstance(p, FunctionalParams):
        raise TypeError("p debe ser FunctionalParams")
    if not isinstance(cfg, DynamicsConfig):
        raise TypeError("cfg debe ser DynamicsConfig")
    if not isinstance(fd, FixedDensityConfig):
        raise TypeError("fd debe ser FixedDensityConfig")
    if cfg.integrator != "clip":
        raise ValueError("densidad fija solo admite integrator='clip'")
    n = w0.shape[0]
    m = n * (n - 1) // 2
    total = float(upper_triangle(w0).sum()) if fd.rho is None else density_to_total(fd.rho, n)
    total = min(max(total, 0.0), float(m))
    dt = step_size(n, p, cfg)
    sched = {int(s) for s in snapshot_steps(cfg.max_steps, cfg.snapshot_schedule)}

    w = project_fixed_density(w0, total)
    g = grad_action(w, p)
    records: dict[str, list[float]] = {}

    def record(w_new: FloatArray, w_old: FloatArray, g_new: FloatArray, mult: float) -> None:
        for key, val in step_scalars(w_new, w_old, p, g_new).items():
            records.setdefault(key, []).append(val)
        records.setdefault("multiplier", []).append(mult)

    snaps: list[FloatArray] = []
    snap_idx: list[int] = []
    record(w, w, g, _multiplier_estimate(w, g))
    if 0 in sched:
        snaps.append(upper_triangle(w))
        snap_idx.append(0)

    status = RunStatus.MAX_STEPS
    steps = 0
    for t in range(1, cfg.max_steps + 1):
        with np.errstate(all="ignore"):
            if not np.all(np.isfinite(g)):
                status = RunStatus.NONFINITE
                break
            w_new, nu = _step_with_shift(w, g, dt, total)
            g_new = grad_action(w_new, p)
            finite = bool(np.all(np.isfinite(w_new)) and np.all(np.isfinite(g_new)))
        if not finite:
            status = RunStatus.NONFINITE
            break
        record(w_new, w, g_new, nu / dt)
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


def uniform_state_threshold(gamma_hat: float, rho: float, n: int) -> float:
    """alpha_c = (1+gamma_hat)(N-2)/(rho (N-4)): inestabilidad de W=rho (diseno §1.9); n >= 5."""
    _check_n(n, 5)
    g = _check_real("gamma_hat", gamma_hat)
    r = _check_real("rho", rho)
    if g < 0.0:
        raise ValueError("gamma_hat debe ser >= 0")
    if not 0.0 < r <= 1.0:
        raise ValueError("rho debe estar en (0,1]")
    return (1.0 + g) * (n - 2) / (r * (n - 4))


def uniform_state_unstable(alpha_hat: float, gamma_hat: float, rho: float, n: int) -> bool:
    """True si alpha_hat rho (N-4)/(N-2) > 1 + gamma_hat (diseno §1.9)."""
    a = _check_real("alpha_hat", alpha_hat)
    return bool(a > uniform_state_threshold(gamma_hat, rho, n))
