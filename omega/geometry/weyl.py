"""Dimensión de Weyl D_W (DESIGN §1.3): N(λ) ~ λ^{D/2} para el Laplaciano del grafo.

Ley de Weyl para el conteo de autovalores (Weyl, 1911); versión discreta: Laplaciano combinatorio
L = D − A del grafo binario sobre la componente gigante de A=(W>w_min) (umbral estricto). Se
exporta además la longitud de Fiedler ξ = λ₂^{-1/2} para el análisis de tamaño finito.
Funciones puras, sin estado global.
"""

from __future__ import annotations

import math
from typing import Literal

import numpy as np

from omega.config.settings11 import WeylConfig
from omega.geometry.dimension import fit_loglog
from omega.geometry.distances import threshold_adjacency
from omega.network.topology import component_labels, giant_component_nodes, submatrix
from omega.network.weights import validate_weight_matrix
from omega.types import DimensionEstimate, FloatArray

__all__ = ["laplacian_spectrum", "weyl_staircase", "weyl_dimension", "fiedler_length", "weyl_and_fiedler"]

_NAN = float("nan")


def laplacian_spectrum(w: FloatArray, kind: Literal["combinatorial", "normalized"]) -> FloatArray:
    """Autovalores ascendentes del Laplaciano de la matriz de pesos/adyacencia `w` (simétrica).

    combinatorial: L = D − W. normalized: I − D^{-1/2} W D^{-1/2} (nodos de grado 0: fila nula).
    """
    validate_weight_matrix(w)
    if kind not in ("combinatorial", "normalized"):
        raise ValueError(f"kind desconocido: {kind!r}")
    deg = w.sum(axis=1)
    if kind == "combinatorial":
        lap = np.diag(deg) - w
    else:
        inv = np.zeros_like(deg)
        pos = deg > 0.0
        inv[pos] = 1.0 / np.sqrt(deg[pos])
        lap = np.diag(pos.astype(np.float64)) - inv[:, None] * w * inv[None, :]
    return np.asarray(np.linalg.eigvalsh(0.5 * (lap + lap.T)), dtype=np.float64)


def weyl_staircase(eigs: FloatArray, rtol: float) -> tuple[FloatArray, FloatArray]:
    """Escalera "upper" sin modo cero: (niveles λ_u, N(λ_u)=#{λ_j>0 : λ_j ≤ λ_u}).

    Se descarta exactamente el menor autovalor (el modo cero); los restantes se agrupan en niveles
    con tolerancia relativa `rtol` (cada nivel es el máximo de su grupo; la cuenta incluye todo el
    grupo). Niveles estrictamente crecientes y cuentas estrictamente crecientes.
    """
    e = np.sort(np.asarray(eigs, dtype=np.float64))
    if e.ndim != 1 or e.size < 1 or not np.all(np.isfinite(e)):
        raise ValueError("eigs debe ser vector 1D no vacio y finito")
    if not rtol >= 0.0:
        raise ValueError("rtol debe ser >= 0")
    pos = e[1:]
    if pos.size == 0:
        return np.zeros(0, dtype=np.float64), np.zeros(0, dtype=np.float64)
    breaks = np.flatnonzero(np.diff(pos) > rtol * np.maximum(np.abs(pos[1:]), 1e-300)) + 1
    ends = np.concatenate((breaks, [pos.size])) - 1  # indice del ultimo elemento de cada grupo
    return np.asarray(pos[ends], dtype=np.float64), np.asarray(ends + 1, dtype=np.float64)


def _giant_adjacency(w: FloatArray, w_min: float) -> tuple[FloatArray, FloatArray]:
    """(w, A) restringidos a la gigante: pesos y adyacencia binaria float."""
    validate_weight_matrix(w)
    a = threshold_adjacency(w, w_min)
    _, labels = component_labels(a)
    nodes = giant_component_nodes(labels)
    ws = submatrix(w, nodes)
    asub = a[np.ix_(nodes, nodes)]
    return ws, asub.astype(np.float64)


def _weyl_core(ws: FloatArray, a: FloatArray, cfg: WeylConfig) -> tuple[DimensionEstimate, FloatArray | None]:
    """(estimacion, espectro combinatorio binario si se calculo con la configuracion binaria/combinatoria)."""
    method = f"weyl_{cfg.laplacian}_{cfg.graph}"
    z = np.zeros(0, dtype=np.float64)
    n = ws.shape[0]
    if n < cfg.min_nodes:
        return DimensionEstimate(_NAN, _NAN, "insufficient_component", None, z, z, z, False, method), None
    m = a if cfg.graph == "binary" else a * ws
    eigs = laplacian_spectrum(m, cfg.laplacian)
    reusable = eigs if (cfg.graph == "binary" and cfg.laplacian == "combinatorial") else None
    levels, counts = weyl_staircase(eigs, cfg.degeneracy_rtol)
    sel = np.flatnonzero((counts >= cfg.count_min) & (counts <= cfg.count_max_frac * n) & (levels > 0.0))
    if sel.size < cfg.min_levels:
        est = DimensionEstimate(_NAN, _NAN, "no_window", None, levels, counts, np.full(levels.shape, _NAN), False, method)
        return est, reusable
    lo, hi = int(sel[0]), int(sel[-1])
    slope, err, r2 = fit_loglog(levels[lo : hi + 1], counts[lo : hi + 1])
    loc = 2.0 * np.gradient(np.log(counts), np.log(levels)) if levels.size >= 2 else np.zeros_like(levels)
    est = DimensionEstimate(
        2.0 * slope, 2.0 * err, "ok", (lo, hi), levels, counts, np.asarray(loc, dtype=np.float64), bool(r2 >= cfg.r2_min), method
    )
    return est, reusable


def _fiedler_from(a: FloatArray, eigs: FloatArray | None) -> float:
    if a.shape[0] < 2:
        return _NAN
    spec = eigs if eigs is not None else laplacian_spectrum(a, "combinatorial")
    lam2 = float(spec[1])
    if not lam2 > 1e-12:
        return _NAN
    return float(1.0 / math.sqrt(lam2))


def weyl_dimension(w: FloatArray, w_min: float, cfg: WeylConfig) -> DimensionEstimate:
    """D_W = 2·pendiente MCO de ln N frente a ln λ en la ventana de cuenta [count_min, count_max_frac·n].

    Gigante con n < cfg.min_nodes: `insufficient_component`. Ventana con menos de `cfg.min_levels`
    niveles: `no_window`. `plateau := R² ≥ r2_min`; stderr = 2·se. method "weyl_<lap>_<graph>".
    """
    ws, a = _giant_adjacency(w, w_min)
    return _weyl_core(ws, a, cfg)[0]


def fiedler_length(w: FloatArray, w_min: float) -> float:
    """ξ = λ₂^{-1/2} (Laplaciano combinatorio binario de la gigante); NaN si n<2 o λ₂ ≤ 0."""
    _, a = _giant_adjacency(w, w_min)
    return _fiedler_from(a, None)


def weyl_and_fiedler(w: FloatArray, w_min: float, cfg: WeylConfig) -> tuple[DimensionEstimate, float]:
    """(`weyl_dimension`, `fiedler_length`) con un unico calculo de la gigante y, si la configuracion de Weyl es
    binaria/combinatorial, un unico espectro (B18). Resultados identicos a las dos funciones por separado."""
    ws, a = _giant_adjacency(w, w_min)
    est, eigs = _weyl_core(ws, a, cfg)
    return est, _fiedler_from(a, eigs)
