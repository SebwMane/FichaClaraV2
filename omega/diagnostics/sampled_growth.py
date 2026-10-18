"""Crecimiento de bolas por BFS muestreado + dimension espectral (prerregistro P1-D.3 §2). Diagnostico, no energia ni certificado.

Version dispersa y muestreada de `ball_growth`: apta para N ~ 2e4-5e4. Solo opera sobre la componente gigante de una
adyacencia binaria simetrica dispersa. Nivel I (coherencia multiescala), Nivel II (estabilidad dimensional con D_B2) y
categoria de informe (sin Nivel III: la mejora a VARIEDAD la hace la herramienta del panel).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components, shortest_path

__all__ = [
    "giant_component",
    "sampled_ball_profile",
    "level1_status",
    "dimension_status",
    "spectral_return",
    "report_category",
    "variety_category",
    "LEVEL2_THRESHOLDS",
]

R_MAX = 200
CV_ZERO = 1e-9
LEVEL2_THRESHOLDS = (0.05, 0.03, 0.08)
NON_GROWTH_TOL = 0.02
CLASS_TOL = 0.25
_CHUNK = 50


def _binary_symmetric(adj: Any) -> csr_array:
    """Binaria, simetrica, sin diagonal."""
    coo = sparse.coo_array(adj)
    rows = np.asarray(coo.row, dtype=np.int64)
    cols = np.asarray(coo.col, dtype=np.int64)
    keep = (rows != cols) & (np.asarray(coo.data) != 0)
    r, c = rows[keep], cols[keep]
    both = sparse.coo_array(
        (np.ones(2 * r.size, dtype=np.float64), (np.concatenate([r, c]), np.concatenate([c, r]))), shape=coo.shape
    )
    out = csr_array(both)  # suma duplicados
    out.data[:] = 1.0
    return out


def giant_component(adj: Any) -> csr_array:
    """Componente gigante de una adyacencia dispersa (simetrizada, binaria, sin diagonal)."""
    a = _binary_symmetric(adj)
    _, lab = connected_components(a, directed=False)
    big = int(np.bincount(lab).argmax())
    idx = np.flatnonzero(lab == big)
    sub = csr_array(a[idx][:, idx])
    return sub


def level1_status(cv: list[float], window: list[int]) -> tuple[str, float | None]:
    """Nivel I (P1-D sin cambios). `cv[r-1]` = CV(r)."""
    if cv[1] < CV_ZERO:
        return "CV_CERO", None
    if len(window) < 2:
        return "SIN_VENTANA", None
    rho = cv[window[-1] - 1] / cv[1]
    return ("HOMOGENEIZA" if rho < 1.0 else "NO_HOMOGENEIZA"), float(rho)


def _dim_class(d: float) -> int | str:
    k = int(round(d))
    return k if (abs(d - k) <= CLASS_TOL and k >= 1) else "NO_ENTERA"


def dimension_status(d_b2: list[float], threshold: float = 0.05) -> dict[str, Any]:
    """Nivel II (§2.3) sobre la sucesion D_B2 con umbral final `threshold` (la tolerancia de no-crecimiento 0.02 es fija)."""
    n = len(d_b2)
    if n < 4:
        return {"status": "SIN_VENTANA_D", "D_conv": None, "class": None, "threshold": threshold}
    delta = [d_b2[i + 1] - d_b2[i] for i in range(n - 1)]  # delta[i-1] = Δ_i
    last, prev, prev2 = delta[-1], delta[-2], delta[-3]
    if abs(last) <= threshold and abs(last) <= abs(prev) + NON_GROWTH_TOL and abs(prev) <= abs(prev2) + NON_GROWTH_TOL:
        status = "CONVERGE"
    elif last > threshold:
        status = "CRUCE"
    else:
        status = "NO_CONVERGE"
    d_conv = float(d_b2[-1])
    return {"status": status, "D_conv": d_conv, "class": _dim_class(d_conv), "threshold": threshold}


def sampled_ball_profile(adj_sparse: Any, rng: np.random.Generator, n_sources: int = 400) -> dict[str, Any]:
    """Perfil de bolas por BFS desde `n_sources` fuentes uniformes sin reemplazo de la componente gigante (§2.1-2.3)."""
    g = giant_component(adj_sparse)
    ng = int(g.shape[0])
    s = min(int(n_sources), ng)
    sources = np.sort(rng.choice(ng, size=s, replace=False)).astype(np.int64)
    counts = np.zeros((s, R_MAX + 2), dtype=np.int64)
    for lo in range(0, s, _CHUNK):
        part = sources[lo : lo + _CHUNK]
        dist = np.asarray(shortest_path(g, method="D", directed=False, unweighted=True, indices=part), dtype=np.float64)
        di = np.minimum(dist, R_MAX + 1).astype(np.int64)
        flat = (np.arange(part.size, dtype=np.int64)[:, None] * (R_MAX + 2) + di).ravel()
        counts[lo : lo + part.size] = np.bincount(flat, minlength=part.size * (R_MAX + 2)).reshape(part.size, R_MAX + 2)
    ball = np.cumsum(counts, axis=1)  # ball[:, r] = |B_r| (r = 0..R_MAX+1)
    r_cap = R_MAX
    for r in range(1, R_MAX + 1):
        if float(np.median(ball[:, r])) >= ng / 2.0:
            r_cap = r
            break
    mean: list[float] = []
    cv: list[float] = []
    med: list[float] = []
    for r in range(1, r_cap + 1):
        b = ball[:, r].astype(np.float64)
        m = float(b.mean())
        mean.append(m)
        cv.append(float(b.std() / m))
        med.append(float(np.median(b)))
    while len(cv) < 2:
        cv.append(0.0)
    window = [r for r in range(2, r_cap + 1) if med[r - 1] <= ng / 4.0]
    st, rho = level1_status(cv, window)
    max_w = max(window) if window else 0
    d_b2 = [
        float(math.log(mean[r + 1] / mean[r - 1]) / math.log((r + 2) / r)) for r in range(2, max_w - 1)
    ]  # r = 2 .. max_w - 2 ; mean[r-1] = m(r)
    level2 = {str(t): dimension_status(d_b2, t) for t in LEVEL2_THRESHOLDS}
    return {
        "n_giant": ng,
        "n_sources": s,
        "r_cap": r_cap,
        "mean": mean,
        "cv": cv[: len(mean)],
        "median": med,
        "window": window,
        "max_window": max_w,
        "level1": {"status": st, "rho": rho, "pass": st in ("HOMOGENEIZA", "CV_CERO")},
        "D_B2": d_b2,
        "level2": level2,
    }


def spectral_return(
    adj_sparse: Any, rng: np.random.Generator, n_sources: int = 64, t_max: int = 4000
) -> dict[str, Any]:
    """Dimension espectral descriptiva (§2.4): paseo perezoso P = (I + D^-1 A)/2; probabilidad de retorno media de
    `n_sources` fuentes; D_s(t) = -2 dln p/dln t en una rejilla logaritmica (factor 1.25). Mediana en la mitad superior."""
    g = giant_component(adj_sparse)
    ng = int(g.shape[0])
    s = min(int(n_sources), ng)
    src = np.sort(rng.choice(ng, size=s, replace=False)).astype(np.int64)
    deg = np.asarray(g.sum(axis=1), dtype=np.float64).ravel()
    p = csr_array(0.5 * (sparse.identity(ng, format="csr") + sparse.diags(1.0 / deg) @ g))
    grid: list[int] = [1]
    x = 1.0
    while True:
        x *= 1.25
        t = int(round(x))
        if t > t_max:
            break
        if t > grid[-1]:
            grid.append(t)
    wanted = set(grid)
    block = np.zeros((ng, s), dtype=np.float64)
    cols = np.arange(s)
    block[src, cols] = 1.0
    ret: dict[int, float] = {}
    for t in range(1, grid[-1] + 1):
        block = np.asarray(p @ block, dtype=np.float64)
        if t in wanted:
            ret[t] = float(block[src, cols].mean())
    tg = [t for t in grid if ret[t] > 0.0]
    ds: list[float] = []
    ts: list[float] = []
    for a, b in zip(tg[:-1], tg[1:]):
        ds.append(float(-2.0 * (math.log(ret[b]) - math.log(ret[a])) / (math.log(b) - math.log(a))))
        ts.append(float(math.sqrt(a * b)))
    upper = ds[len(ds) // 2 :]
    return {
        "t_max": grid[-1],
        "t_grid": tg,
        "p_return": [ret[t] for t in tg],
        "t_mid": ts,
        "D_s": ds,
        "D_s_median_upper": float(np.median(upper)) if upper else None,
    }


def report_category(level1_pass: bool, level2_status: str, dim_class: int | str | None) -> str:
    """Categoria de informe (§2.5) sin Nivel III."""
    if not (level1_pass and level2_status == "CONVERGE") or dim_class is None:
        return "NO_GEOMETRICO"
    if isinstance(dim_class, str):
        return "GEOMETRIA_GRUESA(no entera)"
    return f"GEOMETRIA_GRUESA({dim_class})"


def variety_category(coarse: str, dim_class: int | str | None, certificate_passes: bool | None) -> str | None:
    """Mejora de Nivel III (§2.5): None si la categoria no es gruesa con clase entera."""
    if not coarse.startswith("GEOMETRIA_GRUESA") or not isinstance(dim_class, int):
        return None
    if dim_class <= 2:
        return "VARIEDAD_NO_EVALUABLE"
    if certificate_passes is True:
        return f"GEOMETRIA_VARIEDAD({dim_class})"
    return None
