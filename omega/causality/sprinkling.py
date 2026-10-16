"""Sprinkling uniforme en un diamante causal de Minkowski (bloque C0, Rev. 3 §4).

Las coordenadas solo existen para construir la relacion de orden y como referencia en C0-2.
Ningun estimador de `omega.causality.estimators` las recibe.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

__all__ = ["DIAMOND_HALF_HEIGHT", "sprinkle_diamond", "interval_volume", "diamond_volume", "causal_matrix"]

DIAMOND_HALF_HEIGHT = 0.5  # T: el diamante va de t = -T a t = +T, altura propia tau = 2T = 1


def interval_volume(tau: float, d: int) -> float:
    """Volumen de un intervalo causal de altura propia tau en Minkowski d-dimensional: c_d tau^d."""
    if d < 1:
        raise ValueError("d debe ser >= 1")
    if tau < 0:
        raise ValueError("tau debe ser >= 0")
    c_d = math.pi ** ((d - 1) / 2) / (2 ** (d - 1) * d * math.gamma((d + 1) / 2))
    return float(c_d * tau**d)


def diamond_volume(d: int) -> float:
    """Volumen del diamante usado por `sprinkle_diamond` (tau = 2T)."""
    return interval_volume(2 * DIAMOND_HALF_HEIGHT, d)


def sprinkle_diamond(n: int, d: int, rng: np.random.Generator) -> npt.NDArray[np.float64]:
    """n puntos uniformes en {|x| <= T - |t|}; columna 0 = t. Rechazo en la caja [-T, T]^d."""
    if n < 0:
        raise ValueError("n debe ser >= 0")
    if d < 1:
        raise ValueError("d debe ser >= 1")
    t_half = DIAMOND_HALF_HEIGHT
    out = np.empty((n, d), dtype=np.float64)
    filled = 0
    while filled < n:
        batch = max(64, 2 * (n - filled))
        cand = rng.uniform(-t_half, t_half, size=(batch, d))
        r = np.linalg.norm(cand[:, 1:], axis=1) if d > 1 else np.zeros(batch)
        ok = cand[r <= t_half - np.abs(cand[:, 0])]
        take = min(len(ok), n - filled)
        out[filled : filled + take] = ok[:take]
        filled += take
    return out


def causal_matrix(coords: npt.NDArray[np.float64]) -> npt.NDArray[np.bool_]:
    """C[i, j] = True si i precede estrictamente a j: dt > 0 y dt^2 >= |dx|^2."""
    if coords.ndim != 2:
        raise ValueError("coords debe ser (n, d)")
    t = coords[:, 0]
    dt = t[None, :] - t[:, None]
    if coords.shape[1] > 1:
        x = coords[:, 1:]
        dx2 = np.zeros_like(dt)
        for k in range(x.shape[1]):
            diff = x[None, :, k] - x[:, None, k]
            dx2 += diff * diff
    else:
        dx2 = np.zeros_like(dt)
    rel: npt.NDArray[np.bool_] = (dt > 0) & (dt * dt >= dx2)
    return rel
