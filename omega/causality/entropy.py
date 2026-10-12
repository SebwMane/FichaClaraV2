"""Entropia espaciotemporal de Sorkin (SSEE) sobre ordenes, con vacio de Sorkin-Johnston (Omega-E, E0/E1).

Convenciones (Sorkin y Yazdi 2018): C_xy = 1 si x precede a y; G_R = C/2; Delta = (C - C^T)/2; iDelta hermitica.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

__all__ = ["pauli_jordan", "sj_split", "sj_both", "ssee", "interval_closest"]

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


def ssee(i_delta: CArr, w: CArr, region: npt.NDArray[np.int64], truncate: bool) -> float:
    """S = sum lambda ln|lambda|, lambda autovalores de (iDelta_R)^-1 W_R en la imagen (truncada) de iDelta_R."""
    a = i_delta[np.ix_(region, region)]
    b = w[np.ix_(region, region)]
    nu, v = np.linalg.eigh(a)
    if truncate:
        keep = np.abs(nu) >= math.sqrt(region.size) / (4 * math.pi)
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
