"""Estadistica de resumen para replicas (D-5)."""

from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm

from omega.types import FloatArray

__all__ = ["bootstrap_ci", "proportion_ci"]


def _check_level(level: float) -> None:
    if not 0.0 < level < 1.0:
        raise ValueError(f"level debe estar en (0,1), recibido {level}")


def bootstrap_ci(
    values: FloatArray, rng: np.random.Generator, n_boot: int = 2000, level: float = 0.95
) -> tuple[float, float]:
    """IC percentil bootstrap de la media; values 1D finito y no vacio. Devuelve (inf, sup)."""
    v = np.asarray(values, dtype=np.float64)
    if v.ndim != 1 or v.size == 0:
        raise ValueError("values debe ser vector 1D no vacio")
    if not np.all(np.isfinite(v)):
        raise ValueError("values contiene valores no finitos")
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 1:
        raise ValueError("n_boot debe ser entero >= 1")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng debe ser numpy.random.Generator")
    _check_level(level)
    idx = rng.integers(0, v.size, size=(n_boot, v.size))
    means = v[idx].mean(axis=1)
    lo, hi = np.quantile(means, [(1.0 - level) / 2.0, (1.0 + level) / 2.0])
    return float(lo), float(hi)


def proportion_ci(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    """IC de Wilson para una proporcion k/n; 0<=k<=n, n>=1."""
    if isinstance(k, bool) or isinstance(n, bool) or not isinstance(k, int) or not isinstance(n, int):
        raise TypeError("k y n deben ser enteros")
    if n < 1 or not 0 <= k <= n:
        raise ValueError(f"se exige n>=1 y 0<=k<=n, recibido k={k}, n={n}")
    _check_level(level)
    z = float(norm.ppf(0.5 + level / 2.0))
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)
