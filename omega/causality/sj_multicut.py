"""A-E1' (docs/OMEGA_CONSEJO_REV17_ES0.md §7.2): el instrumento de E1'/E1'' con varios umbrales desde una sola
diagonalizacion, con diagnosticos (modos conservados, margenes al umbral).

Mismo algebra exacta que `sj_fast` (R17-3): M = K^T K real (K = C - C^T), |mu| = sqrt(eig(M))/2, y para cada
multiplicador c del umbral (global c sqrt(N)/(4 pi), region c sqrt(|R|)/(4 pi)):
    |iDelta_t| = U_k diag(mu) U_k^T,  iDelta_t = i A U_k U_k^T,  W_t = (|iDelta_t| + iDelta_t)/2,  A = K/2,
y solo los bloques de la region. Los umbrales anidados comparten la diagonalizacion de M y el producto K_R U.
Si algun autovalor queda a menos de `guard` de un umbral (global o de region), ese valor de S se recalcula por la
ruta de referencia (eigh complejo de pauli_jordan + sj_from_eig + ssee), como el respaldo de `sj_fast`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from omega.causality.entropy import (
    causal_complement,
    interval_balanced,
    pauli_jordan,
    sj_from_eig,
    ssee,
)
from omega.causality.sj_fast import GUARD, SUBSET_FRAC, SUBSET_MIN_N

try:
    import scipy.linalg as _sl
except ImportError:  # pragma: no cover
    _sl = None

__all__ = ["CutResult", "MultiCutResult", "multicut", "reference_s"]

Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]
F64 = npt.NDArray[np.float64]
C128 = npt.NDArray[np.complex128]


@dataclass
class CutResult:
    c_mult: float
    s: float
    n_kept_global: int  # cuenta +-mu por separado (como sj_fast.SJResult.n_kept)
    margin_global: float  # min | |mu| - corte global | (cota inferior si se usa el subconjunto)
    n_kept_region: int  # autovalores de iDelta_t restringido a R con |nu| >= corte de region
    margin_region: float  # min | |nu| - corte de region |
    fallback: bool


@dataclass
class MultiCutResult:
    region: Int
    cuts: list[CutResult] = field(default_factory=list)
    n_perp: int = 0
    s_perp: float | None = None  # S de R^perp con c = c_perp (None si R^perp vacio o no pedido)
    n_kept_perp: int = 0
    margin_perp: float = math.inf


def reference_s(c: Bool, region: Int, c_mult: float) -> float:
    """Ruta de referencia de E1'/E1'': eigh complejo de iDelta, sj_from_eig y ssee truncado."""
    n = c.shape[0]
    mu, u = np.linalg.eigh(pauli_jordan(c))
    idt, w = sj_from_eig(mu, u, c_mult * math.sqrt(n) / (4 * math.pi))
    return ssee(idt, w, region, True, c_mult * math.sqrt(region.size) / (4 * math.pi))


def _spectrum(m: F64, n: int, c_min: float) -> tuple[F64, F64, float]:
    """(mu, U, piso): autopares de M con mu >= SUBSET_FRAC * corte minimo (todos si N < SUBSET_MIN_N)."""
    cut_min = c_min * math.sqrt(n) / (4 * math.pi)
    if n >= SUBSET_MIN_N and _sl is not None:
        lo = (2 * SUBSET_FRAC * cut_min) ** 2
        lam, u = _sl.eigh(m, subset_by_value=(lo, np.inf), driver="evr", check_finite=False)
        floor = SUBSET_FRAC * cut_min  # los no devueltos tienen mu < floor
    else:
        lam, u = np.linalg.eigh(m)
        floor = -math.inf
    mu = np.sqrt(np.maximum(lam, 0.0)) / 2
    return mu, u, floor


def _block(k_idx_u: F64, u_idx: F64, mu_k: F64, keep: npt.NDArray[np.bool_], cut_r: float) -> tuple[float, int, float]:
    """S sobre el bloque idx dados (K U)[idx], U[idx] y los modos conservados; devuelve (S, modos de region, margen)."""
    ur = u_idx[:, keep]
    pr = (ur * mu_k) @ ur.T
    idt_r: C128 = (1j * ((0.5 * k_idx_u[:, keep]) @ ur.T)).astype(np.complex128)
    w_r = 0.5 * (pr + idt_r)
    nu, v = np.linalg.eigh(idt_r)
    margin = float(np.min(np.abs(np.abs(nu) - cut_r)))
    kr = np.abs(nu) >= cut_r
    if not np.any(kr):
        return 0.0, 0, margin
    vk = v[:, kr]
    bk = vk.conj().T @ w_r @ vk
    lam = np.linalg.eigvals(bk / nu[kr][:, None]).real
    lam = lam[np.abs(lam) > 1e-12]
    return float(np.sum(lam * np.log(np.abs(lam)))), int(kr.sum()), margin


def multicut(c: Bool, c_mults: tuple[float, ...] = (0.5, 1.0, 2.0), guard: float = GUARD,
             c_perp: float | None = None) -> MultiCutResult:
    """Region N/8 de E1' (interval_balanced) y S para cada c_mult desde una diagonalizacion de K^T K."""
    c = np.asarray(c, dtype=bool)
    n = c.shape[0]
    region = interval_balanced(c, n // 8)
    out = MultiCutResult(region=region)
    if region.size == 0:
        return out
    cf = c.astype(np.float64)
    k = cf - cf.T
    mults = tuple(c_mults) + ((c_perp,) if c_perp is not None else ())
    mu, u, floor = _spectrum(k.T @ k, n, min(mults))
    kur = k[region] @ u
    ur_all = u[region]
    for cm in c_mults:
        cut = cm * math.sqrt(n) / (4 * math.pi)
        keep = mu >= cut
        mg = float(np.min(np.abs(mu - cut))) if mu.size else math.inf
        mg = min(mg, cut - floor)
        cut_r = cm * math.sqrt(region.size) / (4 * math.pi)
        s, nkr, mr = _block(kur, ur_all, mu[keep], keep, cut_r)
        fb = mg < guard or mr < guard
        if fb:
            s = reference_s(c, region, cm)
        out.cuts.append(CutResult(float(cm), s, int(keep.sum()), mg, nkr, mr, bool(fb)))
    if c_perp is not None:
        rp = causal_complement(c, region)
        out.n_perp = int(rp.size)
        if rp.size:
            keep = mu >= c_perp * math.sqrt(n) / (4 * math.pi)
            s, nkr, mr = _block(k[rp] @ u, u[rp], mu[keep], keep, c_perp * math.sqrt(rp.size) / (4 * math.pi))
            out.s_perp, out.n_kept_perp, out.margin_perp = s, nkr, mr
    return out
