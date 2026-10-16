"""Funcional Omega-C0 (prerregistro §0): S = -sum_{i<j} W_ij psi(C_ij) + kappa sum_i k_i^2, C = W^2; puro."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from omega.types import FloatArray

__all__ = [
    "C0Params",
    "params_from_targets",
    "psi",
    "dpsi",
    "action_c0",
    "grad_c0",
    "action_and_grad",
    "kkt_box",
]

_UPPER_TOL = 1e-12


@dataclass(frozen=True)
class C0Params:
    """psi(c) = a + b log(1+c) - lam c; kappa pesa la competencia de grados. Gauge b = 1."""

    a: float
    lam: float
    kappa: float
    b: float = 1.0

    def __post_init__(self) -> None:
        if not (self.lam > 0.0 and self.kappa > 0.0 and self.b > 0.0):
            raise ValueError("lam, kappa y b deben ser > 0")

    @property
    def c_star(self) -> float:
        """argmax psi sobre c >= 0: max(0, b/lam - 1)."""
        return max(0.0, self.b / self.lam - 1.0)

    @property
    def psi_star(self) -> float:
        c = self.c_star
        return float(self.a + self.b * math.log1p(c) - self.lam * c)

    @property
    def k_star(self) -> float:
        return self.psi_star / (4.0 * self.kappa)

    def lower_bound(self, n: int) -> float:
        """LB = -n psi*^2 / (16 kappa) (cota inferior exacta, C0-T1)."""
        return -n * self.psi_star**2 / (16.0 * self.kappa)


def params_from_targets(c_star: float, k_star: float, a: float) -> C0Params:
    """Parametros con c*, k* dados (b = 1): lam = 2 si c* = 0, si no 1/(1+c*); kappa = psi*/(4 k*)."""
    lam = 2.0 if c_star == 0 else 1.0 / (1.0 + c_star)
    p0 = C0Params(a=a, lam=lam, kappa=1.0)
    psi_s = p0.psi_star
    if not psi_s > 0.0:
        raise ValueError(f"psi* = {psi_s} <= 0 para c*={c_star}, a={a}")
    if not k_star > 0.0:
        raise ValueError("k_star debe ser > 0")
    return C0Params(a=a, lam=lam, kappa=psi_s / (4.0 * k_star))


def psi(c: FloatArray, p: C0Params) -> FloatArray:
    return np.asarray(p.a + p.b * np.log1p(c) - p.lam * c, dtype=np.float64)


def dpsi(c: FloatArray, p: C0Params) -> FloatArray:
    return np.asarray(p.b / (1.0 + c) - p.lam, dtype=np.float64)


def _action_from_c(w: FloatArray, c: FloatArray, p: C0Params) -> float:
    k = w.sum(axis=1)
    return float(-0.5 * np.sum(w * psi(c, p)) + p.kappa * np.dot(k, k))


def action_c0(w: FloatArray, p: C0Params) -> float:
    """S_C0[W] para W simetrica con diagonal 0."""
    return _action_from_c(w, w @ w, p)


def _grad_from_c(w: FloatArray, c: FloatArray, p: C0Params) -> FloatArray:
    k = w.sum(axis=1)
    m = w * dpsi(c, p)
    g = -psi(c, p) - (w @ m + m @ w) + 2.0 * p.kappa * (k[:, None] + k[None, :])
    g = 0.5 * (g + g.T)
    np.fill_diagonal(g, 0.0)
    return np.asarray(g, dtype=np.float64)


def grad_c0(w: FloatArray, p: C0Params) -> FloatArray:
    """Derivada por arista (W_ij y W_ji se mueven juntas): matriz simetrica con diagonal 0."""
    return _grad_from_c(w, w @ w, p)


def action_and_grad(w: FloatArray, p: C0Params) -> tuple[float, FloatArray]:
    c = w @ w
    return _action_from_c(w, c, p), _grad_from_c(w, c, p)


def kkt_box(w: FloatArray, p: C0Params, tol_rel: float = 1e-9) -> dict[str, Any]:
    """Condiciones KKT en la caja [0,1] (minimizacion). Estricto: G < -tol en W=1 y G > tol en W=0."""
    g = grad_c0(w, p)
    n = w.shape[0]
    iu = np.triu_indices(n, k=1)
    wv = w[iu]
    gv = g[iu]
    tol = tol_rel * max(1.0, float(np.max(np.abs(gv))) if gv.size else 1.0)
    up = wv >= 1.0 - _UPPER_TOL
    lo = wv <= _UPPER_TOL
    mid = ~(up | lo)
    viol = np.zeros_like(gv)
    viol[up] = np.maximum(gv[up], 0.0)
    viol[lo] = np.maximum(-gv[lo], 0.0)
    viol[mid] = np.abs(gv[mid])
    max_viol = float(viol.max()) if viol.size else 0.0
    is_kkt = bool(max_viol <= tol)
    strict = bool(
        is_kkt
        and (not up.any() or bool(np.all(gv[up] < -tol)))
        and (not lo.any() or bool(np.all(gv[lo] > tol)))
    )
    return {
        "is_kkt": is_kkt,
        "strict": strict,
        "max_violation": max_viol,
        "tol": float(tol),
        "n_upper": int(up.sum()),
        "n_lower": int(lo.sum()),
        "n_interior": int(mid.sum()),
    }
