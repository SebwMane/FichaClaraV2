"""Dimension efectiva D_eff (M§9-M§10; D-6, D-19, D-22; R3, R6, R7).

Estimador "shell" (D = 1 + dlnS/dlnr, S = numero de nodos en la cascara) por defecto y "ball"
(D = dlnN/dlnr) como validacion cruzada. Si no hay ventana de escala valida el estado es
`no_window` y el valor es NaN: nunca un numero espurio.
"""

from __future__ import annotations

import numpy as np

from omega.config.settings import DimensionConfig
from omega.geometry.distances import ball_counts, giant_nodes_from_distances, validate_distance_matrix
from omega.types import DimensionEstimate, FloatArray, IntArray

__all__ = [
    "radius_grid",
    "mean_profile",
    "scaling_window",
    "local_slopes",
    "fit_loglog",
    "shell_profile",
    "profile_dimension",
    "effective_dimension",
]

_NAN = float("nan")


def _finite_offdiag_values(d: FloatArray) -> FloatArray:
    n = d.shape[0]
    v = d[~np.eye(n, dtype=np.bool_)]
    return np.asarray(v[np.isfinite(v)], dtype=np.float64)


def _is_integer_valued(d: FloatArray, tol: float) -> bool:
    v = _finite_offdiag_values(d)
    return bool(v.size > 0 and np.all(np.abs(v - np.rint(v)) <= tol))


def radius_grid(d: FloatArray, cfg: DimensionConfig) -> FloatArray:
    """Rejilla de radios (D-6): enteros 1..d_max si las distancias son enteras (tol), si no
    `n_radii` valores log-espaciados en [mediana del vecino mas cercano, d_max]."""
    validate_distance_matrix(d)
    v = _finite_offdiag_values(d)
    if v.size == 0:
        return np.zeros(0, dtype=np.float64)
    dmax = float(v.max())
    if _is_integer_valued(d, cfg.integer_tol):
        top = int(np.rint(dmax))
        return np.arange(1, top + 1, dtype=np.float64)
    dd = np.where(np.isfinite(d), d, np.inf)
    np.fill_diagonal(dd, np.inf)
    nn = dd.min(axis=1)
    r_lo = float(np.median(nn[np.isfinite(nn)]))
    if not r_lo < dmax:
        return np.array([dmax], dtype=np.float64)
    return np.asarray(np.geomspace(r_lo, dmax, cfg.n_radii), dtype=np.float64)


def mean_profile(counts: FloatArray | IntArray) -> FloatArray:
    """Perfil medio N_bar(r) = media sobre nodos de N_i(r) (matriz (N,K) -> (K,))."""
    c = np.asarray(counts, dtype=np.float64)
    if c.ndim != 2 or c.size == 0:
        raise ValueError("counts debe ser matriz (N,K) no vacia")
    return np.asarray(c.mean(axis=0), dtype=np.float64)


def _best_run(mask: np.ndarray, scales: FloatArray, min_ratio: float, min_points: int) -> tuple[int, int] | None:
    """Mejor tramo contiguo de `mask` (mayor razon de escala) con razon >= min_ratio y >= min_points."""
    best: tuple[int, int] | None = None
    best_span = -1.0
    n = mask.size
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and mask[j + 1]:
            j += 1
        span = float(np.log(scales[j] / scales[i]))
        if span > best_span:
            best_span, best = span, (i, j)
        i = j + 1
    if best is None:
        return None
    lo, hi = best
    if hi - lo + 1 < min_points or scales[hi] / scales[lo] < min_ratio:
        return None
    return best


def scaling_window(
    scales: FloatArray, profile: FloatArray, upper_limit: float, min_ratio: float, min_points: int
) -> tuple[int, int] | None:
    """Ventana (lo, hi) de indices: perfil >1 y <= upper_limit (saturacion), razon de escala
    >= min_ratio y >= min_points; None si no existe (D-6)."""
    s = np.asarray(scales, dtype=np.float64)
    p = np.asarray(profile, dtype=np.float64)
    if s.ndim != 1 or s.shape != p.shape:
        raise ValueError("scales y profile deben ser vectores 1D del mismo tamano")
    if s.size == 0:
        return None
    if np.any(s <= 0.0) or np.any(np.diff(s) <= 0.0):
        raise ValueError("scales debe ser estrictamente creciente y positiva")
    mask = np.isfinite(p) & (p > 1.0 + 1e-9) & (p <= upper_limit)
    return _best_run(mask, s, min_ratio, min_points)


def local_slopes(x: FloatArray, y: FloatArray) -> FloatArray:
    """Pendiente local d ln y / d ln x (diferencias centrales en la malla log); y>0, x>0."""
    xa = np.asarray(x, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64)
    if xa.shape != ya.shape or xa.ndim != 1 or xa.size < 2:
        raise ValueError("x e y deben ser vectores 1D del mismo tamano (>=2)")
    if np.any(xa <= 0.0) or np.any(ya <= 0.0):
        raise ValueError("x e y deben ser positivos")
    return np.asarray(np.gradient(np.log(ya), np.log(xa)), dtype=np.float64)


def fit_loglog(x: FloatArray, y: FloatArray) -> tuple[float, float, float]:
    """MCO de ln y sobre ln x: (pendiente, error estandar de la pendiente, R^2). stderr=NaN con 2 puntos."""
    xa = np.asarray(x, dtype=np.float64)
    ya = np.asarray(y, dtype=np.float64)
    if xa.shape != ya.shape or xa.ndim != 1 or xa.size < 2:
        raise ValueError("x e y deben ser vectores 1D del mismo tamano (>=2)")
    if np.any(xa <= 0.0) or np.any(ya <= 0.0):
        raise ValueError("x e y deben ser positivos")
    lx, ly = np.log(xa), np.log(ya)
    sxx = float(np.sum((lx - lx.mean()) ** 2))
    if sxx <= 0.0:
        raise ValueError("x sin variacion")
    slope = float(np.sum((lx - lx.mean()) * (ly - ly.mean())) / sxx)
    icpt = float(ly.mean() - slope * lx.mean())
    res = ly - (icpt + slope * lx)
    ss_res = float(np.sum(res**2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 if ss_tot <= 1e-300 else max(0.0, 1.0 - ss_res / ss_tot)
    n = xa.size
    stderr = float(np.sqrt(ss_res / (n - 2) / sxx)) if n > 2 else _NAN
    return slope, stderr, r2


def shell_profile(scales: FloatArray, profile: FloatArray, integer: bool) -> tuple[FloatArray, FloatArray]:
    """Perfil de cascara a partir de N_bar(r) (N_bar(0)=1). Entero: S(r_k)=N_k-N_{k-1} en r_k.
    Continuo: dN/dr ~ DeltaN/Deltar en el punto medio geometrico (K-1 puntos)."""
    s = np.asarray(scales, dtype=np.float64)
    p = np.asarray(profile, dtype=np.float64)
    if s.ndim != 1 or s.shape != p.shape or s.size == 0:
        raise ValueError("scales y profile deben ser vectores 1D no vacios del mismo tamano")
    if integer:
        r_prev = np.concatenate(([0.0], s[:-1]))
        n_prev = np.concatenate(([1.0], p[:-1]))
        return s.copy(), np.asarray((p - n_prev) / (s - r_prev), dtype=np.float64)
    if s.size < 2:
        return np.zeros(0), np.zeros(0)
    mid = np.sqrt(s[1:] * s[:-1])
    return np.asarray(mid, dtype=np.float64), np.asarray(np.diff(p) / np.diff(s), dtype=np.float64)


def _empty(status: str, method: str) -> DimensionEstimate:
    z = np.zeros(0, dtype=np.float64)
    st = "insufficient_component" if status == "insufficient_component" else "no_window"
    return DimensionEstimate(_NAN, _NAN, st, None, z, z, z, False, method)  # type: ignore[arg-type]


def profile_dimension(
    scales: FloatArray,
    profile: FloatArray,
    n_nodes: int,
    cfg: DimensionConfig,
    integer: bool,
    *,
    trim_first: int = 0,
) -> DimensionEstimate:
    """D_eff desde un perfil medio N_bar(scales) con N_gc=n_nodes (D-6, D-19, D-22).
    Ventana por saturacion N_bar<=cfg.saturation*n_nodes, razon>=min_scale_ratio, >=min_points.
    `trim_first` (solo shell con radios enteros): descarta las primeras cascaras de la ventana del
    ajuste (la cascara r=1 es el numero de coordinacion, no escala continua); la ventana devuelta
    es la de los puntos ajustados (>=2) y, con 2 puntos, stderr=NaN."""
    if trim_first < 0:
        raise ValueError("trim_first debe ser >= 0")
    s = np.asarray(scales, dtype=np.float64)
    p = np.asarray(profile, dtype=np.float64)
    limit = cfg.saturation * n_nodes
    if s.size == 0:
        return _empty("no_window", cfg.estimator)
    ball_ok = np.isfinite(p) & (p > 1.0 + 1e-9) & (p <= limit)
    if cfg.estimator == "ball":
        x, y, offset = s, p, 0.0
        mask = ball_ok
    else:
        x, y = shell_profile(s, p, integer)
        offset = 1.0
        if x.size == 0:
            return _empty("no_window", "shell")
        right = ball_ok if integer else ball_ok[1:]
        mask = right & np.isfinite(y) & (y > 0.0)
    win = _best_run(mask, x, cfg.min_scale_ratio, cfg.min_points) if x.size else None
    if win is None:
        return DimensionEstimate(_NAN, _NAN, "no_window", None, x, y, np.full(x.shape, _NAN), False, cfg.estimator)
    lo, hi = win
    if cfg.estimator == "shell" and integer:
        lo += trim_first
        if hi - lo + 1 < 2:
            return DimensionEstimate(_NAN, _NAN, "no_window", None, x, y, np.full(x.shape, _NAN), False, "shell")
    y_pos = np.where(y > 0.0, y, 1e-300)
    loc = offset + local_slopes(x, y_pos) if x.size >= 2 else np.zeros_like(x)
    slope, err, _ = fit_loglog(x[lo : hi + 1], y[lo : hi + 1])
    plateau = bool(np.std(loc[lo : hi + 1]) <= cfg.plateau_tol)
    return DimensionEstimate(offset + slope, err, "ok", (lo, hi), x, y, loc, plateau, cfg.estimator)


def effective_dimension(
    d: FloatArray, cfg: DimensionConfig, *, min_nodes: int = 10, integer_trim: int = 1
) -> DimensionEstimate:
    """D_eff sobre la componente gigante de la matriz de distancias `d` (M§9, D-6, D-22).
    `integer_trim`: cascaras iniciales excluidas del ajuste shell con distancias enteras (0 = literal).
    `insufficient_component` si la gigante tiene < min_nodes nodos; `no_window` si no hay escala."""
    validate_distance_matrix(d)
    nodes = giant_nodes_from_distances(d)
    if nodes.size < min_nodes:
        return _empty("insufficient_component", cfg.estimator)
    dg = d[np.ix_(nodes, nodes)]
    integer = _is_integer_valued(dg, cfg.integer_tol)
    radii = radius_grid(dg, cfg)
    if radii.size == 0:
        return _empty("no_window", cfg.estimator)
    probe = radii + cfg.integer_tol if integer else radii
    prof = mean_profile(ball_counts(dg, probe))
    return profile_dimension(radii, prof, int(nodes.size), cfg, integer, trim_first=integer_trim)
