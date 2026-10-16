"""Redundancia ciclica multiescala como anillos conexos (docs/OMEGA_CIC_L0.md §1 y §3). Diagnostico, no energia.

Para un nodo u y una escala r, el anillo es A_r(u) = {x : r <= d(u, x) <= 2r}; se estudian las componentes conexas del
subgrafo inducido. Sin coordenadas, sin dimension, sin escala ni densidad objetivo.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components, shortest_path

from omega.diagnostics.sampled_growth import giant_component

__all__ = ["annulus_profile", "annulus_status", "MIN_ANNULUS", "F_CONEXO"]

MIN_ANNULUS = 10
F_CONEXO = 0.95
_CHUNK = 50


def _empty(n: int, n_centers: int, r_w: int, scales: list[int]) -> dict[str, Any]:
    return {
        "n_giant": n, "n_centers": n_centers, "r_w": r_w, "scales": scales, "n_pairs": 0, "f_med": None,
        "c_med": None, "f2_med": None, "by_scale": {},
    }


def annulus_profile(
    adj_sparse: Any, rng: np.random.Generator, n_centers: int = 200, max_scales: int = 12
) -> dict[str, Any]:
    """Resumen por grafo de §3: f~ (mediana de f_r(u)), c~ (mediana de c_r(u)), f2~ (fraccion en las 2 mayores)."""
    g: csr_array = giant_component(adj_sparse)
    n = int(g.shape[0])
    if n < 3:
        return _empty(n, 0, 0, [])
    k = min(int(n_centers), n)
    centers = np.sort(rng.choice(n, size=k, replace=False)).astype(np.int64)
    dist = np.empty((k, n), dtype=np.int32)
    for lo in range(0, k, _CHUNK):
        d = shortest_path(g, method="D", directed=False, unweighted=True, indices=centers[lo : lo + _CHUNK])
        dist[lo : lo + d.shape[0]] = np.minimum(d, 2**30).astype(np.int32)
    ecc = int(dist.max())
    # r_w = mayor r con mediana |B_r| <= N/4 (|B_r| es no decreciente en r, y la mediana tambien)
    r_w = 0
    for r in range(1, ecc + 1):
        if float(np.median((dist <= r).sum(axis=1))) <= n / 4.0:
            r_w = r
        else:
            break
    r_hi = r_w // 2
    if r_hi < 2:
        return _empty(n, k, r_w, [])
    allr = np.arange(2, r_hi + 1)
    if allr.size > max_scales:
        pos = np.unique(np.round(np.linspace(0, allr.size - 1, int(max_scales))).astype(np.int64))
        allr = allr[pos]
    scales = [int(r) for r in allr]
    fs: list[float] = []
    cs: list[int] = []
    f2s: list[float] = []
    by_scale: dict[str, list[float]] = {}
    for r in scales:
        fr: list[float] = []
        for i in range(k):
            idx = np.flatnonzero((dist[i] >= r) & (dist[i] <= 2 * r))
            if idx.size < MIN_ANNULUS:
                continue
            sub = csr_array(g[idx][:, idx])
            nc, lab = connected_components(sub, directed=False)
            sz = np.sort(np.bincount(lab, minlength=nc))[::-1]
            f = float(sz[0]) / idx.size
            fs.append(f)
            fr.append(f)
            cs.append(int(nc))
            f2s.append(float(sz[:2].sum()) / idx.size)
        if fr:
            by_scale[str(r)] = [float(np.median(fr)), float(len(fr))]
    if not fs:
        return _empty(n, k, r_w, scales)
    return {
        "n_giant": n, "n_centers": k, "r_w": r_w, "scales": scales, "n_pairs": len(fs), "f_med": float(np.median(fs)),
        "c_med": float(np.median(cs)), "f2_med": float(np.median(f2s)), "by_scale": by_scale,
    }


def annulus_status(summary: dict[str, Any]) -> str:
    """Estado por grafo (§3): SIN_VENTANA, CONEXO, DOS_EXTREMOS o RAMIFICADO."""
    if summary["f_med"] is None:
        return "SIN_VENTANA"
    if float(summary["f_med"]) >= F_CONEXO:
        return "CONEXO"
    if float(summary["c_med"]) == 2.0 and float(summary["f2_med"]) >= F_CONEXO:
        return "DOS_EXTREMOS"
    return "RAMIFICADO"
