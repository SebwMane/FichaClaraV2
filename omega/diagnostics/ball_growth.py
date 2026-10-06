"""Observador de crecimiento de bolas (prerregistro P1-D §2). Diagnostico, no energia ni certificado.

Sobre la componente gigante de una adyacencia binaria: |B_i(r)|, CV entre nodos, D_B(r) y la firma
rho = CV(r_max)/CV(2) en la ventana intermedia W = {r >= 2 : mediana |B_r| <= N_g * frac}.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.sparse import csr_array
from scipy.sparse.csgraph import connected_components, dijkstra

from omega.types import FloatArray

__all__ = ["giant", "ball_profile"]

R_CAP = 15
CV_ZERO = 1e-9


def giant(adj: FloatArray) -> np.ndarray:
    a = np.asarray(adj) > 0
    np.fill_diagonal(a, False)
    a = a | a.T
    _, lab = connected_components(csr_array(a.astype(np.float64)), directed=False)
    big = np.bincount(lab).argmax()
    idx = np.flatnonzero(lab == big)
    return np.asarray(a[np.ix_(idx, idx)], dtype=np.bool_)


def _status(cv: list[float], window: list[int]) -> tuple[str, float | None]:
    if cv[1] < CV_ZERO:
        return "CV_CERO", None
    if len(window) < 2:
        return "SIN_VENTANA", None
    rho = cv[window[-1] - 1] / cv[1]
    return ("HOMOGENEIZA" if rho < 1.0 else "NO_HOMOGENEIZA"), float(rho)


def ball_profile(adj: FloatArray) -> dict[str, Any]:
    a = giant(adj)
    n = a.shape[0]
    d = dijkstra(csr_array(a.astype(np.float64)), directed=False, unweighted=True)
    sizes: list[np.ndarray] = []
    for r in range(1, R_CAP + 1):
        b = (d <= r).sum(axis=1).astype(np.float64)
        sizes.append(b)
        if np.median(b) >= n / 2 or r == R_CAP:
            break
    mean = [float(b.mean()) for b in sizes]
    cv = [float(b.std() / b.mean()) for b in sizes]
    med = [float(np.median(b)) for b in sizes]
    while len(cv) < 2:
        cv.append(0.0)
    d_b = [float(np.log(mean[i + 1] / mean[i]) / np.log((i + 2) / (i + 1))) for i in range(len(mean) - 1)]
    out: dict[str, Any] = {"n_giant": n, "mean": mean, "cv": cv[: len(mean)], "median": med, "D_B": d_b}
    for tag, frac in (("q4", 0.25), ("q2", 0.5)):
        w = [r for r in range(2, len(med) + 1) if med[r - 1] <= n * frac]
        st, rho = _status(cv, w)
        steps = [cv[r - 1] < cv[r - 2] for r in w[1:]]
        out[tag] = {"window": w, "status": st, "rho": rho,
                    "frac_steps_down": float(np.mean(steps)) if steps else None,
                    "D_B_range": [min(d_b[r - 1] for r in w if r - 1 < len(d_b)), max(d_b[r - 1] for r in w if r - 1 < len(d_b))]
                    if any(r - 1 < len(d_b) for r in w) else None}
    return out
