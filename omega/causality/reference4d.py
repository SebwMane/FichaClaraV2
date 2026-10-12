"""Referencia independiente de E[S4_eps] en el diamante 4D por muestreo de pares de puntos (Rev. 3 §19.1).

No construye ninguna matriz de orden: usa coordenadas, volumen del intervalo y la forma cerrada de E_K f4(K, eps).
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from omega.causality.sprinkling import diamond_volume, interval_volume, sprinkle_diamond

__all__ = ["pair_volume_fractions", "expected_action_4d"]


def pair_volume_fractions(pairs: int, rng: np.random.Generator, chunk: int = 1_000_000) -> npt.NDArray[np.float64]:
    """v = V(I(p,q))/V(diamante) para pares relacionados, 0 si no; se genera en tramos (AUDITORIA A-21: memoria)."""
    if pairs > chunk:
        parts = [pair_volume_fractions(min(chunk, pairs - k), rng, chunk) for k in range(0, pairs, chunk)]
        return np.concatenate(parts)
    p = sprinkle_diamond(pairs, 4, rng)
    q = sprinkle_diamond(pairs, 4, rng)
    dt = np.abs(q[:, 0] - p[:, 0])
    dx2 = np.sum((q[:, 1:] - p[:, 1:]) ** 2, axis=1)
    rel = dt * dt >= dx2
    tau = np.sqrt(np.maximum(dt * dt - dx2, 0.0))
    vol: npt.NDArray[np.float64] = np.where(rel, interval_volume(1.0, 4) * tau**4 / diamond_volume(4), 0.0)
    return vol


def _ef4(v: npt.NDArray[np.float64], m: int, eps: float) -> npt.NDArray[np.float64]:
    y = 1 - eps * v
    out: npt.NDArray[np.float64] = (
        y**m
        - 9 * eps * m * v * y ** (m - 1)
        + 8 * eps**2 * m * (m - 1) * v**2 * y ** (m - 2)
        - (4 / 3) * eps**3 * m * (m - 1) * (m - 2) * v**3 * y ** (m - 3)
    )
    return out


def expected_action_4d(n: int, eps: float, v: npt.NDArray[np.float64]) -> tuple[float, float]:
    """(E[S4_eps], SE) a partir de las fracciones v de pares (un par no relacionado aporta 0)."""
    m = n - 2
    terms = np.where(v > 0, _ef4(v, m, eps), 0.0)
    pairs = math.comb(n, 2)
    mean = float(terms.mean())
    se = float(terms.std(ddof=1) / math.sqrt(terms.size))
    # un par ordenado al azar: P(relacionado) cuenta ambos sentidos, y C(n,2) pares no ordenados -> factor 1
    return n - eps * pairs * mean, eps * pairs * se
