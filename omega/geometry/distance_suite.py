"""Suite de distancias de Ω-1.1 (DESIGN §1.1, §3.2): HOP, WEIGHTED_INVERSE, WEIGHTED_LOG y RESISTANCE.

Todas las funciones son puras. La suite opera sobre la componente gigante de A=(W>w_min) (umbral
estricto), igual que `omega.geometry.observables.geometry_observables` (Ω-1.0). Referencias:
Klein & Randić (1993) para la distancia de resistencia; Doyle & Snell (1984) para la analogía
con redes eléctricas; OMEGA_1_1_DESIGN §0.2 (RESISTANCE no vota: se informa ζ_R).
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components, dijkstra

from omega.config.settings import DimensionConfig, GraphConfig
from omega.config.settings11 import GEODESIC_MODES, DistanceMode, DistanceSuiteConfig
from omega.contracts import DistanceSuiteResult
from omega.geometry.dimension import effective_dimension
from omega.geometry.distances import (
    distance_matrix,
    hop_distance_matrix,
    threshold_adjacency,
    validate_distance_matrix,
)
from omega.network.topology import component_labels, giant_component_nodes, submatrix
from omega.network.weights import validate_weight_matrix
from omega.types import BoolArray, DimensionEstimate, FloatArray

__all__ = [
    "edge_length_matrix_mode",
    "resistance_matrix",
    "mode_distance_matrix",
    "mode_dimension",
    "resistance_exponent",
    "metric_spread",
    "distance_suite",
]

_NAN = float("nan")


def _check_adjacency_for(w: FloatArray, a: BoolArray) -> None:
    validate_weight_matrix(w)
    if not isinstance(a, np.ndarray) or a.dtype != np.bool_ or a.shape != w.shape:
        raise ValueError("a debe ser ndarray bool con la misma forma que w")
    if not np.array_equal(a, a.T) or np.any(np.diag(a)):
        raise ValueError("a debe ser simetrica con diagonal False")


def edge_length_matrix_mode(
    w: FloatArray, a: BoolArray, mode: DistanceMode, epsilon: float, log_floor: float
) -> csr_array:
    """Longitudes de arista sobre A segun el modo geodesico (DESIGN §1.1).

    HOP: 1; WEIGHTED_INVERSE: 1/(W+ε); WEIGHTED_LOG: 1 − ln(max(W, log_floor)) (W=1 da 1, W→0 da
    ∞, monótona). RESISTANCE no es geodésica y se rechaza (usar `resistance_matrix`).
    """
    _check_adjacency_for(w, a)
    if not epsilon > 0.0:
        raise ValueError(f"epsilon debe ser > 0, recibido {epsilon}")
    if not 0.0 < log_floor < 1.0:
        raise ValueError(f"log_floor debe estar en (0,1), recibido {log_floor}")
    if mode is DistanceMode.HOP:
        lengths = np.where(a, 1.0, 0.0)
    elif mode is DistanceMode.WEIGHTED_INVERSE:
        lengths = np.where(a, 1.0 / (w + epsilon), 0.0)
    elif mode is DistanceMode.WEIGHTED_LOG:
        lengths = np.where(a, 1.0 - np.log(np.maximum(w, log_floor)), 0.0)
    else:
        raise ValueError(f"modo no geodesico: {mode}")
    return csr_array(np.asarray(lengths, dtype=np.float64))


def resistance_matrix(w: FloatArray, a: BoolArray, max_nodes: int = 2000) -> FloatArray:
    """Resistencia efectiva R_ij = L⁺_ii + L⁺_jj − 2L⁺_ij con L el Laplaciano de W⊙A (DESIGN §1.1).

    Se calcula por componente conexa (pinv hermítica, O(n³)); inf entre componentes distintas y
    0 en la diagonal. Se rechaza (ValueError) si alguna componente tiene más de `max_nodes` nodos.
    """
    _check_adjacency_for(w, a)
    if isinstance(max_nodes, bool) or not isinstance(max_nodes, int) or max_nodes < 2:
        raise ValueError("max_nodes debe ser entero >= 2")
    n = w.shape[0]
    n_comp, labels = connected_components(csr_array(a.astype(np.float64)), directed=False)
    sizes = np.bincount(labels, minlength=n_comp)
    if int(sizes.max()) > max_nodes:
        raise ValueError(f"componente de {int(sizes.max())} nodos excede max_nodes={max_nodes}")
    out = np.full((n, n), np.inf, dtype=np.float64)
    for c in range(n_comp):
        idx = np.flatnonzero(labels == c)
        if idx.size == 1:
            out[idx[0], idx[0]] = 0.0
            continue
        cond = np.where(a[np.ix_(idx, idx)], w[np.ix_(idx, idx)], 0.0)
        lap = np.diag(cond.sum(axis=1)) - cond
        lp = np.linalg.pinv(lap, hermitian=True)
        dg = np.diag(lp)
        r = dg[:, None] + dg[None, :] - 2.0 * lp
        r = np.maximum(0.5 * (r + r.T), 0.0)
        np.fill_diagonal(r, 0.0)
        out[np.ix_(idx, idx)] = r
    return np.asarray(out, dtype=np.float64)


def mode_distance_matrix(
    w: FloatArray, w_min: float, mode: DistanceMode, graph: GraphConfig, suite: DistanceSuiteConfig
) -> FloatArray:
    """Matriz de distancias del modo dado sobre A=(W>w_min); inf entre componentes.

    HOP delega en `hop_distance_matrix` y WEIGHTED_INVERSE en `distance_matrix` (Ω-1.0, bit a bit).
    `graph.epsilon` es ε y `suite.log_floor`/`suite.resistance_max_nodes` los parámetros de LOG/RES.
    """
    a = threshold_adjacency(w, w_min)
    if mode is DistanceMode.HOP:
        return hop_distance_matrix(a)
    if mode is DistanceMode.WEIGHTED_INVERSE:
        return distance_matrix(w, w_min, graph.epsilon)
    if mode is DistanceMode.WEIGHTED_LOG:
        g = edge_length_matrix_mode(w, a, mode, graph.epsilon, suite.log_floor)
        d = np.asarray(dijkstra(g, directed=False), dtype=np.float64)
        return np.asarray(0.5 * (d + d.T), dtype=np.float64)
    if mode is DistanceMode.RESISTANCE:
        return resistance_matrix(w, a, suite.resistance_max_nodes)
    raise ValueError(f"modo desconocido: {mode}")


def _is_integer_valued(d: FloatArray, tol: float) -> bool:
    n = d.shape[0]
    v = d[~np.eye(n, dtype=np.bool_)]
    v = v[np.isfinite(v)]
    return bool(v.size > 0 and np.all(np.abs(v - np.rint(v)) <= tol))


def mode_dimension(
    d: FloatArray, dim: DimensionConfig, suite: DistanceSuiteConfig, integer_metric: bool
) -> tuple[DimensionEstimate, bool]:
    """D_eff de una matriz de distancias con reintento preregistrado (DESIGN §1.1).

    Si la métrica no es entera (`integer_metric` False y distancias no enteras) y la estimación da
    `no_window`, se reintenta UNA vez con n_radii=`suite.fallback_n_radii` y se marca
    `fallback_used`. Devuelve (estimación, fallback_used).
    """
    validate_distance_matrix(d)
    est = effective_dimension(d, dim)
    if est.status != "no_window" or integer_metric or _is_integer_valued(d, dim.integer_tol):
        return est, False
    retry = effective_dimension(d, dataclasses.replace(dim, n_radii=suite.fallback_n_radii))
    return retry, True


def resistance_exponent(r: FloatArray, hop: FloatArray, r_hi: float) -> float:
    """ζ_R = ln(R̄(r_hi)/R̄(1))/ln r_hi, con R̄(s) la resistencia media de los pares a s saltos.

    NaN si r_hi < 2 o si alguna cáscara está vacía o tiene resistencia no finita/no positiva.
    Valores medidos (DESIGN §0): 1D≈0.92–0.97, 2D 0.37–0.48, 3D 0.21–0.30.
    """
    if r.shape != hop.shape or r.ndim != 2:
        raise ValueError("r y hop deben tener la misma forma (N,N)")
    if not math.isfinite(r_hi) or r_hi < 2.0:
        return _NAN
    hi = int(round(r_hi))
    m1 = hop == 1.0
    mh = hop == float(hi)
    if not m1.any() or not mh.any():
        return _NAN
    r1 = float(r[m1].mean())
    rh = float(r[mh].mean())
    if not (math.isfinite(r1) and math.isfinite(rh)) or r1 <= 0.0 or rh <= 0.0:
        return _NAN
    return float(math.log(rh / r1) / math.log(float(hi)))


def metric_spread(estimates: Mapping[DistanceMode, DimensionEstimate]) -> float:
    """Dispersión máxima por pares (max−min) de D entre HOP, INV y LOG; inf si alguno no está "ok"."""
    vals: list[float] = []
    for m in GEODESIC_MODES:
        e = estimates.get(m)
        if e is None or e.status != "ok" or not math.isfinite(e.value):
            return float("inf")
        vals.append(e.value)
    return float(max(vals) - min(vals))




def _empty_estimate(method: str) -> DimensionEstimate:
    z = np.zeros(0, dtype=np.float64)
    return DimensionEstimate(_NAN, _NAN, "insufficient_component", None, z, z, z, False, method)


def distance_suite(
    w: FloatArray, graph: GraphConfig, dim: DimensionConfig, suite: DistanceSuiteConfig
) -> DistanceSuiteResult:
    """Suite completa sobre la componente gigante de A=(W>w_min) (DESIGN §1.1).

    `distance_sensitive` es True salvo que HOP, INV y LOG estén "ok" y su dispersión máxima sea
    ≤ `suite.metric_tol`. Invariantes: estimates[INV] == d_eff de `geometry_observables` y
    estimates[HOP] == d_eff_hops (Ω-1.0).
    """
    validate_weight_matrix(w)
    a = threshold_adjacency(w, graph.w_min)
    _, labels = component_labels(a)
    nodes = giant_component_nodes(labels)
    estimates: dict[DistanceMode, DimensionEstimate] = {}
    fallback: dict[DistanceMode, bool] = {}
    zeta = _NAN
    notes: list[str] = []
    if nodes.size < 2:
        for m in suite.modes:
            estimates[m] = _empty_estimate(dim.estimator)
            fallback[m] = False
        return DistanceSuiteResult(estimates, fallback, _NAN, float("inf"), True, "componente gigante trivial")
    ws = submatrix(w, nodes)
    hop: FloatArray | None = None
    rmat: FloatArray | None = None
    for m in suite.modes:
        if m is DistanceMode.RESISTANCE and nodes.size > suite.resistance_max_nodes:
            z = np.zeros(0, dtype=np.float64)
            estimates[m] = DimensionEstimate(_NAN, _NAN, "no_window", None, z, z, z, False, dim.estimator)
            fallback[m] = False
            notes.append(f"RESISTANCE omitida: n={nodes.size} > {suite.resistance_max_nodes}")
            continue
        d = mode_distance_matrix(ws, graph.w_min, m, graph, suite)
        estimates[m], fallback[m] = mode_dimension(d, dim, suite, m is DistanceMode.HOP)
        if m is DistanceMode.HOP:
            hop = d
        elif m is DistanceMode.RESISTANCE:
            rmat = d
    if rmat is not None:
        hop_est = estimates.get(DistanceMode.HOP)
        if hop is None:
            hop = mode_distance_matrix(ws, graph.w_min, DistanceMode.HOP, graph, suite)
            hop_est = effective_dimension(hop, dim)
        if hop_est is not None and hop_est.status == "ok" and hop_est.window is not None:
            r_hi = float(hop_est.scales[hop_est.window[1]])
            zeta = resistance_exponent(rmat, hop, r_hi)
    spread = metric_spread(estimates)
    sensitive = not (spread <= suite.metric_tol)
    reason = (
        f"dispersion geodesica {spread:.3g} "
        f"{'<=' if not sensitive else '>'} metric_tol {suite.metric_tol:g}"
        if math.isfinite(spread)
        else "algun modo geodesico sin ventana"
    )
    if notes:
        reason += "; " + "; ".join(notes)
    return DistanceSuiteResult(estimates, fallback, zeta, spread, sensitive, reason)
