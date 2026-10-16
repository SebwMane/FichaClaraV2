"""Promedios de fase e^{i lam S} para la rama lorentziana de L-2c (Rev. 3 §10.4)."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

__all__ = ["phase_average", "bipartite_phase_exact"]


def phase_average(actions: Sequence[float], lam: float) -> tuple[float, float]:
    """(|media de e^{i lam S}|, SE aproximado) con SE = sqrt(var cos + var sin) / sqrt(M)."""
    s = np.asarray(actions, dtype=np.float64)
    m = s.size
    if m < 2:
        raise ValueError("se necesitan al menos 2 muestras")
    c, sn = np.cos(lam * s), np.sin(lam * s)
    mag = float(math.hypot(float(c.mean()), float(sn.mean())))
    se = float(math.sqrt(float(c.var(ddof=1) + sn.var(ddof=1)) / m))
    return mag, se


def bipartite_phase_exact(n: int, lam: float) -> float:
    """|cos lam|^{(n//2)(n - n//2)}: valor exacto para el bipartito aleatorio con p = 1/2."""
    return float(abs(math.cos(lam)) ** ((n // 2) * (n - n // 2)))
