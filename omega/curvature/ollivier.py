"""Curvatura de Ollivier-Ricci sobre aristas (OMEGA_1_1_DESIGN §1.6).

kappa(x, y) = 1 - W1(m_x, m_y) / d(x, y), con la medida perezosa
m_x = idleness * delta_x + (1 - idleness) * W_xy / sum_y W_xy sobre los vecinos de x en A, y W1 el
coste de transporte optimo (LP, HiGHS) con coste de suelo d del modo de distancia elegido
(HOP por defecto) sobre la componente gigante. Validaciones exactas: toros hipercubicos 2D/3D
kappa = 0; ciclo C_n (n>=6) kappa = 0; K_n: kappa = 1 - |a - (1-a)/(n-1)|.
Referencias: Ollivier (2009), *Ricci curvature of Markov chains on metric spaces*;
Lin-Lu-Yau (2011) para la version perezosa.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linprog
from scipy.sparse.csgraph import connected_components

from omega.config.settings import GraphConfig
from omega.config.settings11 import CurvatureConfig, DistanceMode, DistanceSuiteConfig
from omega.contracts import CurvatureSummary
from omega.geometry.distances import distance_matrix, hop_distance_matrix, threshold_adjacency
from omega.network.weights import validate_weight_matrix
from omega.types import BoolArray, FloatArray, IntArray

__all__ = ["wasserstein1", "neighbor_measure", "ollivier_edge", "ollivier_curvature"]


def wasserstein1(mu: FloatArray, nu: FloatArray, cost: FloatArray) -> float:
    """W1 entre mu (n) y nu (m) con matriz de costes `cost` (n, m), por LP de transporte (HiGHS).

    mu y nu deben ser no negativas con igual masa total (tol 1e-9). Costes finitos.
    """
    mu = np.asarray(mu, dtype=np.float64)
    nu = np.asarray(nu, dtype=np.float64)
    cost = np.asarray(cost, dtype=np.float64)
    n, m = mu.size, nu.size
    if cost.shape != (n, m):
        raise ValueError(f"cost debe tener forma ({n}, {m}), recibido {cost.shape}")
    if np.any(mu < 0) or np.any(nu < 0) or abs(mu.sum() - nu.sum()) > 1e-9:
        raise ValueError("mu y nu deben ser no negativas con la misma masa")
    if not np.all(np.isfinite(cost)):
        raise ValueError("cost debe ser finito")
    a_eq = np.zeros((n + m, n * m), dtype=np.float64)
    for i in range(n):
        a_eq[i, i * m : (i + 1) * m] = 1.0
    for j in range(m):
        a_eq[n + j, j::m] = 1.0
    res = linprog(cost.ravel(), A_eq=a_eq, b_eq=np.concatenate([mu, nu]), bounds=(0, None), method="highs")
    if not res.success:
        raise RuntimeError(f"linprog fallo: {res.message}")
    return float(res.fun)


def neighbor_measure(w: FloatArray, a: BoolArray, x: int, idleness: float) -> tuple[IntArray, FloatArray]:
    """Soporte y masas de m_x: idleness en x y (1-idleness) repartido proporcional a W sobre vecinos.

    Un nodo aislado en A concentra toda su masa en x.
    """
    if not 0.0 <= idleness < 1.0:
        raise ValueError("idleness debe estar en [0,1)")
    nb = np.flatnonzero(a[x])
    if nb.size == 0:
        return np.asarray([x], dtype=np.int64), np.asarray([1.0], dtype=np.float64)
    ww = w[x, nb]
    tot = float(ww.sum())
    support = np.concatenate([[x], nb]).astype(np.int64)
    mass = np.concatenate([[idleness], (1.0 - idleness) * ww / tot]).astype(np.float64)
    return support, mass


def ollivier_edge(w: FloatArray, a: BoolArray, d: FloatArray, x: int, y: int, idleness: float) -> float:
    """kappa(x, y) = 1 - W1(m_x, m_y)/d(x, y) con coste de suelo d."""
    dxy = float(d[x, y])
    if not (np.isfinite(dxy) and dxy > 0.0):
        raise ValueError("d(x, y) debe ser finito y > 0")
    sx, mx = neighbor_measure(w, a, x, idleness)
    sy, my = neighbor_measure(w, a, y, idleness)
    return 1.0 - wasserstein1(mx, my, d[np.ix_(sx, sy)]) / dxy


def _mode_distances(
    w: FloatArray, a: BoolArray, graph: GraphConfig, suite: DistanceSuiteConfig, mode: DistanceMode
) -> FloatArray:
    if mode is DistanceMode.HOP:
        return hop_distance_matrix(a)
    if mode is DistanceMode.WEIGHTED_INVERSE:
        return distance_matrix(w, graph.w_min, graph.epsilon)
    from omega.geometry.distance_suite import mode_distance_matrix  # import diferido (WP-B)

    return mode_distance_matrix(w, graph.w_min, mode, graph, suite)


def ollivier_curvature(
    w: FloatArray,
    w_min: float,
    graph: GraphConfig,
    suite: DistanceSuiteConfig,
    cfg: CurvatureConfig,
    rng: np.random.Generator,
) -> CurvatureSummary:
    """Curvatura de Ollivier sobre las aristas de la componente gigante de {W > w_min}.

    Si hay mas de cfg.max_edges aristas se muestrean sin reemplazo con `rng` (sampled=True).
    Salida: media, desviacion (ddof=1), SE = std/sqrt(n), tail_fraction = frac(kappa <
    cfg.tail_cut) y ok = (media >= cfg.mean_min) y (tail_fraction <= cfg.tail_frac_max).
    Sin aristas: n_edges=0, media 0 y ok=False. `graph.w_min` se sustituye por `w_min`.
    """
    validate_weight_matrix(w)
    a_full = threshold_adjacency(w, w_min)
    n_comp, labels = connected_components(a_full.astype(np.int8), directed=False)
    giant = np.flatnonzero(labels == np.argmax(np.bincount(labels)))
    ws = np.ascontiguousarray(w[np.ix_(giant, giant)])
    a = np.ascontiguousarray(a_full[np.ix_(giant, giant)])
    xs, ys = np.nonzero(np.triu(a, k=1))
    m = int(xs.size)
    if m == 0:
        return CurvatureSummary(mean=0.0, std=0.0, se=0.0, tail_fraction=0.0, n_edges=0, sampled=False,
                                edge_values=np.zeros(0, dtype=np.float64), ok=False)
    sampled = m > cfg.max_edges
    if sampled:
        pick = np.sort(rng.choice(m, size=cfg.max_edges, replace=False))
        xs, ys = xs[pick], ys[pick]
    g = GraphConfig(w_min=w_min, epsilon=graph.epsilon, w_min_sensitivity=graph.w_min_sensitivity)
    d = _mode_distances(ws, a, g, suite, cfg.distance_mode)
    vals = np.asarray([ollivier_edge(ws, a, d, int(x), int(y), cfg.idleness)
                       for x, y in zip(xs, ys, strict=True)], dtype=np.float64)
    n = int(vals.size)
    std = float(vals.std(ddof=1)) if n > 1 else 0.0
    mean = float(vals.mean())
    tail = float(np.mean(vals < cfg.tail_cut))
    return CurvatureSummary(
        mean=mean, std=std, se=std / float(np.sqrt(n)), tail_fraction=tail, n_edges=n, sampled=sampled,
        edge_values=vals, ok=bool(mean >= cfg.mean_min and tail <= cfg.tail_frac_max),
    )
