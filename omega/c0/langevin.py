"""Langevin con reflexion para S_C0 (prerregistro Fase 2, §2.1).

u <- B(u - dt*G + sqrt(2*Theta*dt)*xi) sobre una variable por par; estacionaria ~ exp(-S_C0/Theta).
Theta es temperatura estadistica, no tiempo fisico.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from omega.c0.functional import C0Params, _action_from_c, _grad_from_c
from omega.types import FloatArray

__all__ = ["reflect_unit", "background_fraction", "langevin_c0"]


def reflect_unit(x: FloatArray) -> FloatArray:
    y = np.mod(np.asarray(x, dtype=np.float64), 2.0)
    return np.asarray(np.clip(np.where(y <= 1.0, y, 2.0 - y), 0.0, 1.0), dtype=np.float64)


def background_fraction(w: FloatArray, rel: float = 0.1) -> float:
    """Fraccion de la fuerza total en pares con W < rel * max W."""
    total = float(w.sum())
    if total <= 0.0:
        return 0.0
    return float(w[w < rel * float(w.max())].sum() / total)


def langevin_c0(
    w0: FloatArray,
    p: C0Params,
    theta: float,
    *,
    dt: float,
    n_steps: int,
    rng: np.random.Generator,
    record_every: int = 1000,
) -> dict[str, Any]:
    if theta < 0.0 or not dt > 0.0:
        raise ValueError("theta >= 0 y dt > 0")
    n = int(w0.shape[0])
    iu = np.triu_indices(n, 1)
    w = np.array(w0, dtype=np.float64, copy=True)
    np.fill_diagonal(w, 0.0)
    amp = float(np.sqrt(2.0 * theta * dt))
    trace: list[dict[str, float]] = []
    for t in range(1, n_steps + 1):
        g = _grad_from_c(w, w @ w, p)
        u = w[iu] - dt * g[iu]
        if amp > 0.0:
            u = u + amp * rng.standard_normal(u.shape[0])
        u = reflect_unit(u)
        w = np.zeros((n, n), dtype=np.float64)
        w[iu] = u
        w = w + w.T
        if t % record_every == 0 or t == n_steps:
            trace.append({"step": float(t), "S": _action_from_c(w, w @ w, p), "f_bg": background_fraction(w)})
    return {"w": w, "trace": trace, "s_final": trace[-1]["S"] if trace else _action_from_c(w, w @ w, p)}
