"""Estimadores que reciben SOLO la matriz de orden C (C[i, j] = i precede a j)."""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt
from scipy.optimize import brentq

__all__ = [
    "mm_fraction",
    "ordering_fraction",
    "mm_dimension",
    "between_counts",
    "interval_abundances",
    "bd_action_2d",
    "bd_action_2d_smeared",
    "bd_action_4d",
    "bd_action_4d_smeared",
    "smearing_weight_4d",
    "midpoint_dimension",
]

BoolMatrix = npt.NDArray[np.bool_]
D_MIN, D_MAX = 1.0, 12.0


def mm_fraction(d: float) -> float:
    """r(d) = Gamma(d+1) Gamma(d/2) / (2 Gamma(3d/2)): fraccion esperada de pares relacionados."""
    return float(math.exp(math.lgamma(d + 1) + math.lgamma(d / 2) - math.lgamma(1.5 * d)) / 2)


def ordering_fraction(c: BoolMatrix) -> float:
    n = int(c.shape[0])
    if n < 2:
        raise ValueError("se necesitan al menos 2 elementos")
    return float(np.count_nonzero(c)) / (n * (n - 1) / 2)


def mm_dimension(c: BoolMatrix) -> float:
    """Invierte r(d) en [D_MIN, D_MAX]; nan si la fraccion cae fuera del rango."""
    r = ordering_fraction(c)
    lo, hi = mm_fraction(D_MAX), mm_fraction(D_MIN)
    if not (lo < r <= hi):
        return float("nan")
    if r == hi:
        return D_MIN
    return float(brentq(lambda d: mm_fraction(d) - r, D_MIN, D_MAX, xtol=1e-10))


def between_counts(c: BoolMatrix) -> npt.NDArray[np.int64]:
    """B[i, j] = elementos estrictamente entre i y j si i precede a j; -1 si no estan relacionados."""
    cf = c.astype(np.float32)
    b: npt.NDArray[np.int64] = np.rint(cf @ cf).astype(np.int64)
    b[~c] = -1
    return b


def interval_abundances(c: BoolMatrix, kmax: int) -> npt.NDArray[np.int64]:
    """N_k para k = 0..kmax: pares relacionados con exactamente k elementos entre ellos."""
    b = between_counts(c)
    vals = b[b >= 0]
    return np.bincount(vals, minlength=kmax + 1)[: kmax + 1].astype(np.int64)


def bd_action_2d(c: BoolMatrix) -> int:
    """Accion de Benincasa-Dowker 2D sin suavizar: S = N - 2 N0 + 4 N1 - 2 N2."""
    n0, n1, n2 = (int(x) for x in interval_abundances(c, 2))
    return int(c.shape[0]) - 2 * n0 + 4 * n1 - 2 * n2


def midpoint_dimension(c: BoolMatrix) -> float:
    """En el intervalo mas grande I(a, b), m maximiza min(|I(a,m)|, |I(m,b)|); d = log2(|I(a,b)| / min)."""
    b = between_counts(c)
    if np.all(b < 0):
        return float("nan")
    a, z = np.unravel_index(int(np.argmax(b)), b.shape)
    total = int(b[a, z])
    inside = np.flatnonzero(c[a] & c[:, z])
    if total < 2 or inside.size == 0:
        return float("nan")
    half = np.minimum(b[a, inside], b[inside, z])
    best = int(half.max())
    if best <= 0:
        return float("nan")
    return math.log2(total / best)


def smearing_weight(n: npt.NDArray[np.float64], eps: float) -> npt.NDArray[np.float64]:
    """f(n, eps) = (1-eps)^n [1 - 2 eps n/(1-eps) + eps^2 n(n-1) / (2(1-eps)^2)]  (Rev. 3 §11.2)."""
    z = 1.0 - eps
    out: npt.NDArray[np.float64] = z**n * (1 - 2 * eps * n / z + eps * eps * n * (n - 1) / (2 * z * z))
    return out


def bd_action_2d_smeared(c: BoolMatrix, eps: float) -> float:
    """S_eps = eps (N - 2 eps sum_n N_n f(n, eps)); eps -> 1 recupera bd_action_2d."""
    if not 0 < eps < 1:
        raise ValueError("eps debe estar en (0, 1)")
    b = between_counts(c)
    vals = b[b >= 0].astype(np.float64)
    return float(eps * (c.shape[0] - 2 * eps * float(np.sum(smearing_weight(vals, eps)))))


def bd_action_4d(c: BoolMatrix) -> int:
    """Accion BD 4D sin suavizar (prefactor absorbido): S4 = N - N0 + 9 N1 - 16 N2 + 8 N3."""
    n0, n1, n2, n3 = (int(x) for x in interval_abundances(c, 3))
    return int(c.shape[0]) - n0 + 9 * n1 - 16 * n2 + 8 * n3


def smearing_weight_4d(n: npt.NDArray[np.float64], eps: float) -> npt.NDArray[np.float64]:
    """f4(n, eps) de Benincasa-Dowker (Rev. 3 §19)."""
    z = 1.0 - eps
    out: npt.NDArray[np.float64] = z**n * (
        1 - 9 * eps * n / z + 8 * eps**2 * n * (n - 1) / z**2 - 4 * eps**3 * n * (n - 1) * (n - 2) / (3 * z**3)
    )
    return out


def bd_action_4d_smeared(c: BoolMatrix, eps: float) -> float:
    """S4_eps = sqrt(eps) (N - eps sum_n N_n f4(n, eps)) (BD 2010 ec. 8; Dowker-Glaser 2013 ec. 42, d = 4).

    Prefactor (4/sqrt(6)) (l/l_p)^2 absorbido como en bd_action_4d; eps = (l/l_k)^4; eps = 1 da bd_action_4d.
    Corregido en la auditoria de sqrt(eps) (docs/omega_c/AUDITORIA_SQRT_EPS.md); antes devolvia S4_eps / sqrt(eps).
    """
    if eps == 1.0:
        return float(bd_action_4d(c))
    if not 0 < eps < 1:
        raise ValueError("eps debe estar en (0, 1]")
    b = between_counts(c)
    vals = b[b >= 0].astype(np.float64)
    return float(math.sqrt(eps) * (c.shape[0] - eps * float(np.sum(smearing_weight_4d(vals, eps)))))
