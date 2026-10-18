"""Solapamiento de difusion alpha_tau y escalado difusivo por arista (prerregistro Omega-D, secciones 1 y 4.2)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from omega.types import BoolArray, FloatArray, IntArray

DEFAULT_TAUS: tuple[int, ...] = (1, 2, 4, 8, 16, 32)
DEGENERATE_EPS = 1e-12


def giant_dense(adj: FloatArray | BoolArray) -> BoolArray:
    """Adyacencia binaria simetrica sin diagonal restringida a la componente conexa mayor."""
    a = np.asarray(adj, dtype=np.float64) > 0.0
    a = a | a.T
    np.fill_diagonal(a, False)
    _, labels = connected_components(csr_matrix(a), directed=False)
    big = int(np.argmax(np.bincount(labels)))
    idx = np.flatnonzero(labels == big)
    return np.asarray(a[np.ix_(idx, idx)], dtype=np.bool_)


def lazy_sym_operator(a_bool: BoolArray) -> FloatArray:
    """M = 1/2 (I + K^-1/2 A K^-1/2) con K = diag(grados). Requiere grado > 0 en todo nodo."""
    a = np.asarray(a_bool, dtype=np.float64)
    k = a.sum(axis=1)
    if np.any(k <= 0.0):
        raise ValueError("nodo aislado: usar giant_dense primero")
    s = 1.0 / np.sqrt(k)
    m = 0.5 * (a * s[:, None] * s[None, :])
    m[np.diag_indices_from(m)] += 0.5
    return np.asarray(m, dtype=np.float64)


def overlap_powers(m: FloatArray, taus: Sequence[int] = DEFAULT_TAUS) -> dict[int, FloatArray]:
    """alpha_tau = M^{2 tau}_ij / sqrt(M^{2 tau}_ii M^{2 tau}_jj), con M^{2 tau} por cuadrados sucesivos.

    Los tau deben ser potencias de 2 (M^2, M^4, ... se obtienen elevando al cuadrado).
    """
    want = sorted(set(int(t) for t in taus))
    for t in want:
        if t < 1 or (t & (t - 1)) != 0:
            raise ValueError("tau debe ser potencia de 2")
    out: dict[int, FloatArray] = {}
    p = np.asarray(m @ m, dtype=np.float64)  # M^2 = M^{2 tau} con tau = 1
    tau = 1
    while tau <= want[-1]:
        if tau in want:
            d = np.sqrt(np.diag(p))
            al = p / (d[:, None] * d[None, :])
            out[tau] = np.asarray(np.clip(al, 0.0, 1.0), dtype=np.float64)
        p = np.asarray(p @ p, dtype=np.float64)
        tau *= 2
    return out


def edge_scaling(
    alphas: dict[int, FloatArray],
    edges_i: IntArray,
    edges_j: IntArray,
    fit_taus: Sequence[int] = (2, 4, 8),
) -> dict[str, Any]:
    """Q_tau = -tau ln alpha por arista, pendiente s_e de ln Q vs ln tau en fit_taus y resumen por grafo."""
    taus = sorted(alphas)
    ae = {t: np.minimum(alphas[t][edges_i, edges_j], 1.0) for t in taus}
    q: dict[int, FloatArray] = {}
    degen_t: dict[int, BoolArray] = {}
    for t in taus:
        a = ae[t]
        dg = np.asarray(a >= 1.0 - DEGENERATE_EPS, dtype=np.bool_)
        with np.errstate(divide="ignore", invalid="ignore"):
            qt = -t * np.log(np.where(dg, 0.5, np.maximum(a, 1e-300)))
        qt = np.where(dg | (qt <= DEGENERATE_EPS), 0.0, qt)
        dg = dg | (qt <= DEGENERATE_EPS)
        q[t] = np.asarray(qt, dtype=np.float64)
        degen_t[t] = dg
    degenerate = np.zeros(edges_i.shape[0], dtype=np.bool_)
    for t in fit_taus:
        degenerate |= degen_t[t]
    lt = np.log(np.asarray(fit_taus, dtype=np.float64))
    lt_c = lt - lt.mean()
    lq = np.stack([np.log(np.where(degenerate, 1.0, np.maximum(q[t], 1e-300))) for t in fit_taus], axis=1)
    slope = (lq * lt_c[None, :]).sum(axis=1) / float((lt_c**2).sum())
    slope = np.where(degenerate, np.nan, slope)
    nd = ~degenerate
    n_e = int(edges_i.shape[0])
    summary: dict[str, Any] = {
        "n_edges": n_e,
        "median_s": float(np.median(slope[nd])) if nd.any() else float("nan"),
        "f_nd": float(np.mean(degenerate | (np.nan_to_num(slope, nan=0.0) < -0.5))),
        "frac_degenerate": float(degenerate.mean()),
        "median_Q": {str(t): float(np.median(q[t])) for t in taus},
        "median_alpha": {str(t): float(np.median(ae[t])) for t in taus},
        "p10_alpha": {str(t): float(np.percentile(ae[t], 10)) for t in taus},
    }
    return {"Q": q, "slope": slope, "degenerate": degenerate, "summary": summary}


def distance2_contrast(
    a_bool: BoolArray,
    alphas: dict[int, FloatArray],
    rng: np.random.Generator,
    n_pairs: int = 5000,
) -> dict[str, Any]:
    """Mediana de alpha_tau sobre pares (i<j) a distancia de saltos exactamente 2 (muestra de hasta n_pairs)."""
    a = np.asarray(a_bool, dtype=np.float64)
    two = (a @ a) > 0.0
    two &= ~np.asarray(a_bool, dtype=np.bool_)
    np.fill_diagonal(two, False)
    iu, ju = np.nonzero(np.triu(two, 1))
    total = int(iu.size)
    if total == 0:
        return {"n_pairs": 0, "n_available": 0, "median_alpha": {str(t): float("nan") for t in alphas}}
    if total > n_pairs:
        pick = np.sort(rng.choice(total, size=n_pairs, replace=False))
        iu, ju = iu[pick], ju[pick]
    return {
        "n_pairs": int(iu.size),
        "n_available": total,
        "median_alpha": {str(t): float(np.median(alphas[t][iu, ju])) for t in sorted(alphas)},
    }
