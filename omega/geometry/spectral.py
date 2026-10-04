"""Dimension espectral D_s por paseo aleatorio lazy (M§11; D-7, D-12; R4).

P = (1-q) I + q D^-1 W sobre la componente gigante. Retorno medio exacto
P_bar(t) = (1/N) sum_k mu_k^t (eigh de la forma simetrica). D_s = -2 dlnP_bar/dlnt.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components

from omega.config.settings import SpectralConfig
from omega.geometry.dimension import _best_run, fit_loglog, local_slopes
from omega.network.weights import validate_weight_matrix
from omega.types import DimensionEstimate, FloatArray

__all__ = [
    "walk_operator_spectrum",
    "heat_spectrum",
    "time_grid",
    "mean_return_probability",
    "spectral_dimension",
]

_NAN = float("nan")


def _sym_normalized(w: FloatArray) -> FloatArray:
    """D^-1/2 W D^-1/2 (exige fuerza > 0 en todos los nodos)."""
    validate_weight_matrix(w)
    k = w.sum(axis=1)
    if np.any(k <= 0.0):
        raise ValueError("hay nodos aislados (fuerza 0); restringir a la componente gigante")
    inv = 1.0 / np.sqrt(k)
    return np.asarray(w * inv[:, None] * inv[None, :], dtype=np.float64)


def walk_operator_spectrum(w: FloatArray, laziness: float) -> FloatArray:
    """Autovalores (ascendentes) de (1-q)I + q D^-1/2 W D^-1/2, en [1-2q, 1] (D-12)."""
    if not 0.0 < laziness <= 1.0:
        raise ValueError(f"laziness debe estar en (0,1], recibido {laziness}")
    s = _sym_normalized(w)
    m = (1.0 - laziness) * np.eye(w.shape[0]) + laziness * s
    return np.asarray(np.linalg.eigvalsh(m), dtype=np.float64)


def heat_spectrum(w: FloatArray) -> FloatArray:
    """Autovalores (ascendentes) del Laplaciano normalizado I - D^-1/2 W D^-1/2, en [0, 2] (D-12)."""
    s = _sym_normalized(w)
    return np.asarray(np.linalg.eigvalsh(np.eye(w.shape[0]) - s), dtype=np.float64)


def time_grid(n: int, cfg: SpectralConfig) -> FloatArray:
    """Tiempos enteros log-espaciados en [1, 10 n] (a lo sumo n_times, sin repetidos) (D-7)."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError(f"n debe ser entero >= 1, recibido {n!r}")
    t = np.unique(np.rint(np.geomspace(1.0, 10.0 * n, cfg.n_times)))
    return np.asarray(t, dtype=np.float64)


def mean_return_probability(spectrum: FloatArray, times: FloatArray, method: str) -> FloatArray:
    """P_bar(t) = (1/N) sum_k f(spectrum_k, t); lazy_walk: mu^t (t entero); heat_normalized: exp(-lambda t)."""
    mu = np.asarray(spectrum, dtype=np.float64)
    t = np.asarray(times, dtype=np.float64)
    if mu.ndim != 1 or mu.size == 0 or t.ndim != 1 or np.any(t < 0.0):
        raise ValueError("spectrum y times deben ser vectores 1D (t >= 0)")
    out = np.empty(t.size, dtype=np.float64)
    if method == "lazy_walk":
        if np.any(t != np.rint(t)):
            raise ValueError("lazy_walk exige tiempos enteros")
        for i, ti in enumerate(t):
            out[i] = float(np.sum(np.power(mu, ti))) / mu.size
    elif method == "heat_normalized":
        for i, ti in enumerate(t):
            out[i] = float(np.sum(np.exp(-mu * ti))) / mu.size
    else:
        raise ValueError(f"method desconocido: {method!r}")
    return out


def _insufficient(method: str) -> DimensionEstimate:
    z = np.zeros(0, dtype=np.float64)
    return DimensionEstimate(_NAN, _NAN, "insufficient_component", None, z, z, z, False, method)


def spectral_dimension(
    w: FloatArray,
    cfg: SpectralConfig,
    *,
    plateau_tol: float = 0.3,
    min_points: int = 3,
    min_nodes: int = 10,
) -> DimensionEstimate:
    """D_s sobre la componente gigante del soporte de `w` (D-7, D-12). Ventana: t>=t_min,
    P_bar>=saturation_factor/N_gc, razon>=min_scale_ratio, >=min_points; si no, `no_window`.
    `w` ya debe estar umbralizado/binarizado segun cfg.graph (lo hace observables)."""
    validate_weight_matrix(w)
    n_all = w.shape[0]
    _, labels = connected_components(csr_array(w > 0.0), directed=False)
    giant = int(np.argmax(np.bincount(labels)))
    nodes = np.flatnonzero(labels == giant)
    n = int(nodes.size)
    if n < min_nodes:
        return _insufficient(cfg.method)
    sub = w[np.ix_(nodes, nodes)] if n < n_all else w
    if cfg.method == "lazy_walk":
        spec = walk_operator_spectrum(sub, cfg.laziness)
        spec = np.maximum(spec, 0.0) if cfg.laziness <= 0.5 else spec
    else:
        spec = heat_spectrum(sub)
    t = time_grid(n, cfg)
    p = mean_return_probability(spec, t, cfg.method)
    pos = p > 0.0
    mask = pos & (t >= cfg.t_min) & (p >= cfg.saturation_factor / n)
    win = _best_run(mask, t, cfg.min_scale_ratio, min_points)
    if win is None:
        return DimensionEstimate(_NAN, _NAN, "no_window", None, t, p, np.full(t.shape, _NAN), False, cfg.method)
    lo, hi = win
    loc = -2.0 * local_slopes(t[pos], p[pos]) if pos.all() else np.full(t.shape, _NAN)
    slope, err, _ = fit_loglog(t[lo : hi + 1], p[lo : hi + 1])
    plateau = bool(np.std(loc[lo : hi + 1]) <= plateau_tol)
    return DimensionEstimate(-2.0 * slope, 2.0 * err, "ok", (lo, hi), t, p, loc, plateau, cfg.method)
