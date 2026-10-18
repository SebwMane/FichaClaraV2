"""Descenso de gradiente proyectado en la caja [0,1] con paso adaptativo monotono (prerregistro §0)."""

from __future__ import annotations

from typing import Any

import numpy as np

from omega.c0.functional import C0Params, _action_from_c, _grad_from_c
from omega.types import FloatArray

__all__ = ["evolve_c0"]

DT_MIN = 1e-14


def _sym0(w: FloatArray) -> FloatArray:
    out = 0.5 * (w + w.T)
    np.fill_diagonal(out, 0.0)
    return np.asarray(out, dtype=np.float64)


def evolve_c0(
    w0: FloatArray,
    p: C0Params,
    *,
    max_steps: int = 20000,
    tol: float = 1e-10,
    patience: int = 50,
    dt0: float | None = None,
) -> dict[str, Any]:
    """Minimiza S_C0 en la caja; S no sube (si sube, dt/2 y se reintenta; si se acepta, dt*1.1)."""
    n = int(w0.shape[0])
    w = np.clip(_sym0(np.asarray(w0, dtype=np.float64)), 0.0, 1.0)
    dt = dt0 if dt0 is not None else 1.0 / (4.0 * p.kappa * n + 3.0 * n * (p.b + p.lam) + 1.0)
    c = w @ w
    s = _action_from_c(w, c, p)
    g = _grad_from_c(w, c, p)
    calm = 0
    rejections = 0
    status = "max_steps"
    steps = 0
    while steps < max_steps:
        accepted = False
        while True:
            wn = _sym0(np.clip(w - dt * g, 0.0, 1.0))
            cn = wn @ wn
            sn = _action_from_c(wn, cn, p)
            if sn <= s:
                accepted = True
                break
            rejections += 1
            dt *= 0.5
            if dt < DT_MIN:
                break
        if not accepted:
            status = "stalled"
            break
        steps += 1
        delta = float(np.max(np.abs(wn - w)))
        w, c, s = wn, cn, sn
        g = _grad_from_c(w, c, p)
        dt *= 1.1
        calm = calm + 1 if delta < tol else 0
        if calm >= patience:
            status = "converged"
            break
    resid = float(np.max(np.abs(w - np.clip(w - g, 0.0, 1.0))))  # gradiente proyectado (paso unitario)
    return {
        "w": w, "status": status, "steps": steps, "rejections": rejections,
        "s_final": s, "dt_final": dt, "kkt_residual": resid,
    }
