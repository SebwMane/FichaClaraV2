"""Langevin sobre el vector triangular u (Omega-1.1, DESIGN §1.8).

u <- B(u - dt*g + sqrt(2*Theta*dt)*xi), xi ~ N(0,1); densidad estacionaria ~ exp(-H/Theta).

IMPORTANTE: Theta = Theta_hat * beta es una TEMPERATURA ESTADISTICA. El paso dt es un
parametro numerico del integrador; nada de esto es tiempo fisico (D-1, P§8).
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from omega.config.settings import FunctionalParams
from omega.config.settings11 import Engine, LangevinConfig
from omega.contracts import ChainResult
from omega.dynamics.functional import action, triangles
from omega.dynamics.gradient import grad_action, lipschitz_bound
from omega.network.weights import from_upper_triangle, validate_weight_matrix
from omega.types import FloatArray

__all__ = ["reflect_unit", "temperature", "langevin_step", "chain_observables", "run_langevin"]

OBSERVABLE_KEYS = ("action", "mean_weight", "binary_density", "triangle_density", "strength_cv")


def reflect_unit(x: FloatArray) -> FloatArray:
    """Reflexion en [0,1]: y = x mod 2; B(x) = y si y<=1, 2-y si no. Pura (copia)."""
    y = np.mod(np.asarray(x, dtype=np.float64), 2.0)
    out = np.where(y <= 1.0, y, 2.0 - y)
    return np.asarray(np.clip(out, 0.0, 1.0), dtype=np.float64)


def temperature(theta_hat: float, beta: float) -> float:
    """Theta = Theta_hat * beta (temperatura estadistica, escala del termino beta*w^2)."""
    if not (theta_hat > 0.0) or not (beta > 0.0):
        raise ValueError("theta_hat y beta deben ser > 0")
    return float(theta_hat * beta)


def _require_eta0(p: FunctionalParams) -> None:
    if p.eta != 0.0:
        raise ValueError("Langevin/Metropolis exigen eta == 0")


def langevin_step(
    u: FloatArray,
    n: int,
    p: FunctionalParams,
    dt: float,
    theta: float,
    rng: np.random.Generator,
    boundary: Literal["reflect", "project"],
) -> FloatArray:
    """Un paso Euler-Maruyama sobre u (aristas a<b); devuelve un vector nuevo en [0,1]."""
    _require_eta0(p)
    if boundary not in ("reflect", "project"):
        raise ValueError("boundary debe ser 'reflect' o 'project'")
    if not (dt > 0.0) or theta < 0.0:
        raise ValueError("dt debe ser > 0 y theta >= 0")
    u = np.asarray(u, dtype=np.float64)
    w = from_upper_triangle(u, n)
    g = grad_action(w, p)[np.triu_indices(n, 1)]
    xi = rng.standard_normal(u.shape[0])
    y = u - dt * g + np.sqrt(2.0 * theta * dt) * xi
    if boundary == "reflect":
        return reflect_unit(y)
    return np.asarray(np.clip(y, 0.0, 1.0), dtype=np.float64)


def chain_observables(w: FloatArray, p: FunctionalParams, w_min: float) -> dict[str, float]:
    """action, mean_weight, binary_density (W>w_min), triangle_density=T/C(N,3), strength_cv."""
    validate_weight_matrix(w)
    n = w.shape[0]
    iu = np.triu_indices(n, 1)
    v = w[iu]
    k = w.sum(axis=1)
    kbar = float(k.mean())
    n3 = n * (n - 1) * (n - 2) / 6.0
    return {
        "action": float(action(w, p)),
        "mean_weight": float(v.mean()),
        "binary_density": float(np.mean(v > w_min)),
        "triangle_density": float(triangles(w) / n3) if n >= 3 else 0.0,
        "strength_cv": float(k.std() / kbar) if kbar > 0.0 else 0.0,
    }


def pick_state_indices(n_samples: int, n_states: int) -> set[int]:
    """Indices (de muestra) equiespaciados de los estados a guardar."""
    if n_states <= 0 or n_samples <= 0:
        return set()
    return {int(i) for i in np.linspace(0, n_samples - 1, min(n_states, n_samples)).astype(int)}


def run_langevin(
    w0: FloatArray,
    p: FunctionalParams,
    cfg: LangevinConfig,
    w_min: float,
    rng: np.random.Generator,
    n_states: int = 0,
) -> ChainResult:
    """Cadena Langevin; dt = dt_safety/L(n,p). Guarda observables cada `thin` pasos.

    `states` contiene hasta n_states matrices equiespaciadas en la parte post burn-in.
    """
    validate_weight_matrix(w0)
    _require_eta0(p)
    n = w0.shape[0]
    theta = temperature(cfg.theta_hat, p.beta)
    dt = cfg.dt_safety / lipschitz_bound(n, p)
    u = w0[np.triu_indices(n, 1)].copy()
    steps = list(range(cfg.thin, cfg.n_steps + 1, cfg.thin))
    n_burn = int(cfg.burn_in_fraction * cfg.n_steps)
    post = [i for i, s in enumerate(steps) if s > n_burn]
    keep_rel = pick_state_indices(len(post), n_states)
    keep = {post[i] for i in keep_rel}
    cols: dict[str, list[float]] = {k: [] for k in OBSERVABLE_KEYS}
    states: list[FloatArray] = []
    idx = 0
    for t in range(1, cfg.n_steps + 1):
        u = langevin_step(u, n, p, dt, theta, rng, cfg.boundary)
        if t % cfg.thin == 0:
            w = from_upper_triangle(u, n)
            for key, val in chain_observables(w, p, w_min).items():
                cols[key].append(val)
            if idx in keep:
                states.append(w)
            idx += 1
    return ChainResult(
        engine=Engine.LANGEVIN,
        theta=theta,
        samples={k: np.asarray(v, dtype=np.float64) for k, v in cols.items()},
        sample_steps=np.asarray(steps, dtype=np.int64),
        w_final=from_upper_triangle(u, n),
        states=tuple(states),
        acceptance=1.0,
        step_size=float(dt),
        n_steps=cfg.n_steps,
    )
