"""Entropia espaciotemporal de Sorkin (SSEE) sobre ordenes, con vacio de Sorkin-Johnston (Omega-E, E0/E1).

Convenciones (Sorkin y Yazdi 2018): C_xy = 1 si x precede a y; G_R = C/2; Delta = (C - C^T)/2; iDelta hermitica.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

__all__ = [
    "pauli_jordan",
    "sj_split",
    "sj_both",
    "ssee",
    "interval_closest",
    "interval_balanced",
    "causal_complement",
    "knee_cutoff",
    "knee_powerlaw",
    "sj_from_eig",
    "knee_from_spectrum",
]

CArr = npt.NDArray[np.complex128]
KERNEL_TOL = 1e-8


def pauli_jordan(c: npt.NDArray[np.bool_]) -> CArr:
    """iDelta = (i/2)(C - C^T)."""
    cf = c.astype(np.float64)
    out: CArr = (0.5j * (cf - cf.T)).astype(np.complex128)
    return out


def sj_split(i_delta: CArr, cutoff: float) -> tuple[CArr, CArr]:
    """(iDelta truncado, W_SJ truncado): elimina modos con |mu| < cutoff; W = parte positiva."""
    mu, u = np.linalg.eigh(i_delta)
    keep = np.abs(mu) >= cutoff
    pos = keep & (mu > 0)
    idt: CArr = (u[:, keep] * mu[keep]) @ u[:, keep].conj().T
    w: CArr = (u[:, pos] * mu[pos]) @ u[:, pos].conj().T
    return idt, w


def sj_both(i_delta: CArr, cutoff: float) -> tuple[CArr, CArr, CArr]:
    """Una sola diagonalizacion: (W sin truncar, iDelta truncado, W truncado)."""
    mu, u = np.linalg.eigh(i_delta)
    pos = mu > 0
    w_full: CArr = (u[:, pos] * mu[pos]) @ u[:, pos].conj().T
    keep = np.abs(mu) >= cutoff
    kp = keep & pos
    idt: CArr = (u[:, keep] * mu[keep]) @ u[:, keep].conj().T
    w_t: CArr = (u[:, kp] * mu[kp]) @ u[:, kp].conj().T
    return w_full, idt, w_t


def ssee(i_delta: CArr, w: CArr, region: npt.NDArray[np.int64], truncate: bool, cutoff: float | None = None) -> float:
    """S = sum lambda ln|lambda|, lambda autovalores de (iDelta_R)^-1 W_R en la imagen (truncada) de iDelta_R."""
    a = i_delta[np.ix_(region, region)]
    b = w[np.ix_(region, region)]
    nu, v = np.linalg.eigh(a)
    if truncate:
        thr = cutoff if cutoff is not None else math.sqrt(region.size) / (4 * math.pi)
        keep = np.abs(nu) >= thr
    else:
        keep = np.abs(nu) > KERNEL_TOL * max(float(np.abs(nu).max()), 1e-300)
    if not np.any(keep):
        return 0.0
    vk = v[:, keep]
    bk = vk.conj().T @ b @ vk
    lam = np.linalg.eigvals(bk / nu[keep][:, None]).real
    lam = lam[np.abs(lam) > 1e-12]
    return float(np.sum(lam * np.log(np.abs(lam))))


def interval_closest(c: npt.NDArray[np.bool_], target: int) -> npt.NDArray[np.int64]:
    """Intervalo I[a,b] con extremos cuyo tamano es el mas cercano a target; empate: menores indices (a, b)."""
    cf = c.astype(np.float32)
    between = np.rint(cf @ cf).astype(np.int64) + 2
    between[~c] = -(10**9)
    diff = np.abs(between - target)
    a, b = np.unravel_index(int(np.argmin(diff)), diff.shape)
    members = np.flatnonzero(c[a] & c[:, b])
    out: npt.NDArray[np.int64] = np.sort(np.concatenate([[a, b], members])).astype(np.int64)
    return out


def interval_balanced(c: npt.NDArray[np.bool_], target: int, tol: float = 0.05) -> npt.NDArray[np.int64]:
    """E1': entre los intervalos I[a,b] con |I| en target*(1 +- tol), el que maximiza min(|pasado(a)|, |futuro(b)|).

    Empate: menores indices (a, b). Si no hay candidatos, devuelve un arreglo vacio.
    """
    cf = c.astype(np.float32)
    size = np.rint(cf @ cf).astype(np.int64) + 2
    ok = c & (size >= math.floor(target * (1 - tol))) & (size <= math.ceil(target * (1 + tol)))
    if not np.any(ok):
        return np.zeros(0, dtype=np.int64)
    past = c.sum(axis=0)
    fut = c.sum(axis=1)
    score = np.where(ok, np.minimum(past[:, None], fut[None, :]), -1)
    a, b = np.unravel_index(int(np.argmax(score)), score.shape)
    members = np.flatnonzero(c[a] & c[:, b])
    out: npt.NDArray[np.int64] = np.sort(np.concatenate([[a, b], members])).astype(np.int64)
    return out


def causal_complement(c: npt.NDArray[np.bool_], region: npt.NDArray[np.int64]) -> npt.NDArray[np.int64]:
    """R^perp: elementos fuera de R sin relacion con ningun elemento de R."""
    related: npt.NDArray[np.bool_] = np.asarray(c[region, :].any(axis=0) | c[:, region].any(axis=1), dtype=bool)
    related[region] = True
    out: npt.NDArray[np.int64] = np.flatnonzero(~related).astype(np.int64)
    return out


def knee_cutoff(i_delta: CArr, min_dist: float = 0.05) -> float | None:
    """Rodilla del espectro positivo en log10 (rango) frente a log10 (autovalor): punto de maxima distancia
    a la recta entre el primer y el ultimo punto. Devuelve el autovalor de la rodilla o None si no hay rodilla."""
    mu = np.linalg.eigvalsh(i_delta)
    pos = np.sort(mu[mu > 1e-12])[::-1]
    if pos.size < 3:
        return None
    x = np.log10(np.arange(1, pos.size + 1))
    y = np.log10(pos)
    p0, p1 = np.array([x[0], y[0]]), np.array([x[-1], y[-1]])
    d = p1 - p0
    dist = np.abs(d[0] * (y - p0[1]) - d[1] * (x - p0[0])) / float(np.hypot(*d))
    k = int(np.argmax(dist))
    return float(pos[k]) if float(dist[k]) >= min_dist else None


def knee_powerlaw(i_delta: CArr, n_elements: int) -> tuple[float | None, float]:
    """R-rodilla-2 (Consejo Rev. 4 §8.1). Devuelve (umbral o None si no hay rodilla, pendiente beta ajustada)."""
    mu = np.linalg.eigvalsh(i_delta)
    pos = np.sort(mu[mu > 1e-12])[::-1]
    k_ref = max(6, math.ceil(math.sqrt(n_elements) / 4))
    if pos.size <= k_ref + 3:
        return None, float("nan")
    ranks = np.arange(1, pos.size + 1, dtype=np.float64)
    x = np.log(ranks[1:k_ref])
    y = np.log(pos[1:k_ref])
    beta, alpha = np.polyfit(x, y, 1)
    dev = np.abs(np.log(pos) - (alpha + beta * np.log(ranks))) > math.log(1.3)
    for k in range(k_ref, pos.size - 2):
        if dev[k] and dev[k + 1] and dev[k + 2]:
            return float(pos[k]), float(beta)
    return None, float(beta)


def sj_from_eig(mu: npt.NDArray[np.float64], u: CArr, cutoff: float) -> tuple[CArr, CArr]:
    """(iDelta truncado, W truncado) a partir de una diagonalizacion ya hecha de iDelta."""
    keep = np.abs(mu) >= cutoff
    kp = keep & (mu > 0)
    idt: CArr = (u[:, keep] * mu[keep]) @ u[:, keep].conj().T
    w: CArr = (u[:, kp] * mu[kp]) @ u[:, kp].conj().T
    return idt, w


def knee_from_spectrum(mu: npt.NDArray[np.float64], n_elements: int) -> tuple[float | None, float]:
    """R-rodilla-2 sobre autovalores ya calculados (mismo algoritmo que knee_powerlaw)."""
    pos = np.sort(mu[mu > 1e-12])[::-1]
    k_ref = max(6, math.ceil(math.sqrt(n_elements) / 4))
    if pos.size <= k_ref + 3:
        return None, float("nan")
    ranks = np.arange(1, pos.size + 1, dtype=np.float64)
    beta, alpha = np.polyfit(np.log(ranks[1:k_ref]), np.log(pos[1:k_ref]), 1)
    dev = np.abs(np.log(pos) - (alpha + beta * np.log(ranks))) > math.log(1.3)
    for k in range(k_ref, pos.size - 2):
        if dev[k] and dev[k + 1] and dev[k + 2]:
            return float(pos[k]), float(beta)
    return None, float(beta)
