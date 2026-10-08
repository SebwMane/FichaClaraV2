"""Coherencia multiescala de vecindades (docs/OMEGA_COH_L0.md §1 y §3.1). Diagnostico, no energia ni certificado.

Para una arista e = (u, v): I_r(e) = 1 - |B_r(u) & B_r(v)| / |B_r(u) | B_r(v)|; gamma_e es la pendiente de
-ln I_r frente a ln r en r in [2, r_w]. Sin coordenadas, sin dimension objetivo, sin valores objetivo de I.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.sparse import coo_array, csr_array
from scipy.sparse.csgraph import connected_components, shortest_path

from omega.types import FloatArray, IntArray

__all__ = ["edge_coherence_profile", "coherence_status", "R_CAP_MAX", "GAMMA_LOW"]

R_CAP_MAX = 200
GAMMA_LOW = 0.5
MIN_WINDOW = 4
CHUNK = 20  # aristas por llamada a shortest_path (2 * CHUNK filas de N)


def _giant_edges(adj: Any) -> tuple[IntArray, IntArray, int, csr_array]:
    """Aristas (u < v) de la componente gigante, reetiquetadas 0..Ng-1, Ng y el subgrafo."""
    a = csr_array(adj)
    _, lab = connected_components(a, directed=False)
    big = int(np.bincount(lab).argmax())
    idx = np.flatnonzero(lab == big)
    g = csr_array(a[idx][:, idx])
    sub = coo_array(g)
    keep = sub.row < sub.col
    u = np.asarray(sub.row[keep], dtype=np.int64)
    v = np.asarray(sub.col[keep], dtype=np.int64)
    return u, v, int(idx.size), g


def _cum(counts: IntArray, r_max: int) -> IntArray:
    c = np.zeros(r_max + 1, dtype=np.int64)
    m = min(counts.size, r_max + 1)
    c[:m] = counts[:m]
    return np.asarray(np.cumsum(c), dtype=np.int64)


def edge_coherence_profile(adj_sparse: Any, rng: np.random.Generator, n_edges: int = 400) -> dict[str, Any]:
    """Perfil de coherencia por arista (§3.1). Devuelve el resumen por grafo (sin estado; ver `coherence_status`)."""
    u_all, v_all, n, g = _giant_edges(adj_sparse)
    m = int(u_all.size)
    if m == 0 or n < 3:
        return _summary(n, 0, 0, 0, None, 1.0, 0.0, None, np.empty(0))
    k = min(int(n_edges), m)
    sel = rng.choice(m, size=k, replace=False)
    flip = rng.random(k) < 0.5
    eu = np.where(flip, v_all[sel], u_all[sel])
    ev = np.where(flip, u_all[sel], v_all[sel])
    # perfiles acumulados por arista en r = 0..R_CAP_MAX (se recorta luego)
    r_max = R_CAP_MAX
    size_u = np.zeros((k, r_max + 1), dtype=np.int64)
    size_v = np.zeros((k, r_max + 1), dtype=np.int64)
    inter = np.zeros((k, r_max + 1), dtype=np.int64)
    union = np.zeros((k, r_max + 1), dtype=np.int64)
    for s in range(0, k, CHUNK):
        e = slice(s, min(s + CHUNK, k))
        cu, cv = eu[e], ev[e]
        src = np.concatenate([cu, cv])
        d = shortest_path(g, method="D", directed=False, unweighted=True, indices=src)
        du, dv = d[: cu.size], d[cu.size :]
        mx = np.maximum(du, dv)
        mn = np.minimum(du, dv)
        for j in range(cu.size):
            row = s + j
            size_u[row] = _cum(np.bincount(du[j].astype(np.int64)), r_max)
            size_v[row] = _cum(np.bincount(dv[j].astype(np.int64)), r_max)
            inter[row] = _cum(np.bincount(mx[j].astype(np.int64)), r_max)
            union[row] = _cum(np.bincount(mn[j].astype(np.int64)), r_max)
    med = np.median(size_u, axis=0)
    ge_half = np.flatnonzero(med >= n / 2.0)
    r_cap = int(min(R_CAP_MAX, ge_half[0] if ge_half.size else R_CAP_MAX))
    r_cap = max(r_cap, 1)
    le_q = np.flatnonzero(med[1 : r_cap + 1] <= n / 4.0)
    r_w = int(le_q[-1] + 1) if le_q.size else 0
    ii = 1.0 - inter[:, : r_cap + 1] / np.maximum(union[:, : r_cap + 1], 1)
    median_i1 = float(np.median(ii[:, 1]))
    if r_w < MIN_WINDOW:
        return _summary(n, k, r_cap, r_w, None, 0.0, 0.0, median_i1, np.empty(0))
    win = ii[:, 1 : r_w + 1]
    degenerate = np.any(win <= 0.0, axis=1)
    f_deg = float(np.mean(degenerate))
    gam = _slopes(ii[~degenerate][:, 2 : r_w + 1], r_w)
    med_gamma = float(np.median(gam)) if gam.size else None
    f_low = float(np.mean(gam < GAMMA_LOW)) if gam.size else 0.0
    return _summary(n, k, r_cap, r_w, med_gamma, f_deg, f_low, median_i1, gam)


def _slopes(i_win: FloatArray, r_w: int) -> FloatArray:
    """Pendiente por minimos cuadrados de -ln I frente a ln r, r = 2..r_w, por fila."""
    if i_win.shape[0] == 0:
        return np.empty(0, dtype=np.float64)
    x = np.log(np.arange(2, r_w + 1, dtype=np.float64))
    y = -np.log(i_win)
    xc = x - x.mean()
    return np.asarray((y - y.mean(axis=1, keepdims=True)) @ xc / float(xc @ xc), dtype=np.float64)


def _summary(
    n: int, k: int, r_cap: int, r_w: int, med_gamma: float | None, f_deg: float, f_low: float,
    median_i1: float | None, gam: FloatArray,
) -> dict[str, Any]:
    q = [float(x) for x in np.quantile(gam, [0.1, 0.5, 0.9])] if gam.size else None
    return {
        "n_giant": n, "n_edges": k, "r_cap": r_cap, "r_w": r_w, "median_gamma": med_gamma, "f_deg": f_deg,
        "f_low": f_low, "median_I1": median_i1, "gamma_q10_50_90": q,
    }


def coherence_status(summary: dict[str, Any]) -> str:
    """Estado por grafo segun la tabla de §3.1 (umbrales fijados a priori)."""
    if int(summary["r_w"]) < MIN_WINDOW:
        return "SIN_VENTANA"
    f_deg, f_low = float(summary["f_deg"]), float(summary["f_low"])
    if f_deg > 0.5:
        return "DEGENERADO"
    mg = summary["median_gamma"]
    if mg is None:
        return "DEGENERADO"
    if mg >= 0.7 and f_low <= 0.10 and f_deg <= 0.10:
        return "COHERENTE"
    if mg < 0.5 or f_low > 0.25:
        return "INCOHERENTE"
    return "INTERMEDIO"
