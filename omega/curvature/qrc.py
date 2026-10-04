"""Curvatura de Ricci cuantica de Klitgaard-Loll, solo informe (OMEGA_1_1_DESIGN §1.6).

Para x, y a distancia delta: d_bar(S_delta(x), S_delta(y)) / delta, con S_delta(x) la esfera de
radio delta y d_bar la distancia media entre todos los pares (a, b) de ambas esferas. En una
variedad plana el cociente tiende a una constante; en curvatura positiva es menor.
No es utilizable con delta<=3 (artefactos de reticula: en T2, d_bar/delta = 1.88/1.68/1.65,
hallazgo 6): NO entra en el certificado. Referencia: Klitgaard & Loll (2018),
*Introducing quantum Ricci curvature*, Phys. Rev. D 97, 046008.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from omega.contracts import QRCProfile
from omega.types import FloatArray

__all__ = ["average_sphere_distance", "quantum_ricci_profile"]


def average_sphere_distance(d: FloatArray, x: int, y: int, delta: int) -> float:
    """Distancia media entre S_delta(x) y S_delta(y) (esferas por distancia exacta == delta).

    NaN si alguna esfera es vacia o si la media no es finita.
    """
    if isinstance(delta, bool) or not isinstance(delta, (int, np.integer)) or delta < 1:
        raise ValueError("delta debe ser entero >= 1")
    sx = np.flatnonzero(d[x] == delta)
    sy = np.flatnonzero(d[y] == delta)
    if sx.size == 0 or sy.size == 0:
        return float("nan")
    m = float(np.mean(d[np.ix_(sx, sy)]))
    return m if np.isfinite(m) else float("nan")


def quantum_ricci_profile(
    d: FloatArray, deltas: Sequence[int], max_pairs: int, rng: np.random.Generator
) -> QRCProfile:
    """Perfil d_bar/delta para cada delta, sobre hasta `max_pairs` pares (x<y, d(x,y)=delta).

    Los pares se muestrean sin reemplazo con `rng`. Se informan media, error estandar (ddof=1;
    0 con <2 pares) y numero de pares validos; sin pares la media es NaN.
    """
    if not isinstance(d, np.ndarray) or d.ndim != 2 or d.shape[0] != d.shape[1]:
        raise ValueError("d debe ser matriz cuadrada")
    if isinstance(max_pairs, bool) or not isinstance(max_pairs, int) or max_pairs < 1:
        raise ValueError("max_pairs debe ser entero >= 1")
    means: list[float] = []
    ses: list[float] = []
    counts: list[int] = []
    for delta in deltas:
        xs, ys = np.nonzero(np.triu(d == delta, k=1))
        if xs.size > max_pairs:
            pick = rng.choice(xs.size, size=max_pairs, replace=False)
            xs, ys = xs[pick], ys[pick]
        vals = np.asarray([average_sphere_distance(d, int(x), int(y), int(delta)) / delta
                           for x, y in zip(xs, ys, strict=True)], dtype=np.float64)
        vals = vals[np.isfinite(vals)]
        counts.append(int(vals.size))
        means.append(float(vals.mean()) if vals.size else float("nan"))
        ses.append(float(vals.std(ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0)
    return QRCProfile(
        deltas=np.asarray(list(deltas), dtype=np.int64),
        ratio_mean=np.asarray(means, dtype=np.float64),
        ratio_se=np.asarray(ses, dtype=np.float64),
        n_pairs=np.asarray(counts, dtype=np.int64),
    )
