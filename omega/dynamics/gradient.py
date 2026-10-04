"""Gradiente exacto por arista a<b (ANALYSIS §1.2, §1.5, §1.6; M§17, D-3, D-13, D-25).

Todas las funciones devuelven matrices simetricas con diagonal 0: G_ab = dS/dw_ab.
"""

from __future__ import annotations

import numpy as np

from omega.config.settings import FunctionalParams
from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray

__all__ = [
    "grad_triangles",
    "grad_density",
    "grad_degree_irregularity",
    "grad_smoothness",
    "grad_mass",
    "grad_action",
    "lipschitz_bound",
    "theta_gradient",
]


def _zero_diag(g: FloatArray) -> FloatArray:
    out = np.array(g, dtype=np.float64, copy=True)
    np.fill_diagonal(out, 0.0)
    return out


def grad_triangles(w: FloatArray) -> FloatArray:
    """dT/dw_ab = (W^2)_ab."""
    validate_weight_matrix(w)
    return _zero_diag(w @ w)


def grad_density(w: FloatArray) -> FloatArray:
    """dS_dens/dw_ab = 2 W_ab."""
    validate_weight_matrix(w)
    return _zero_diag(2.0 * w)


def grad_degree_irregularity(w: FloatArray) -> FloatArray:
    """dS_deg/dw_ab = 2[(k_a - k_medio) + (k_b - k_medio)] (exacto, §1.2)."""
    validate_weight_matrix(w)
    k = w.sum(axis=1)
    u = k - k.mean()
    return _zero_diag(2.0 * (u[:, None] + u[None, :]))


def grad_smoothness(w: FloatArray) -> FloatArray:
    """dS_smooth/dw_ab = 2(s_a+s_b) + 4 W_ab (k_a+k_b) - 12 (W^2)_ab, s=diag(W^2) (M§15)."""
    validate_weight_matrix(w)
    w2 = w @ w
    s = np.diag(w2)
    k = w.sum(axis=1)
    g = 2.0 * (s[:, None] + s[None, :]) + 4.0 * w * (k[:, None] + k[None, :]) - 12.0 * w2
    return _zero_diag(g)


def grad_mass(w: FloatArray) -> FloatArray:
    """dS_mass/dw_ab = 1 fuera de la diagonal."""
    validate_weight_matrix(w)
    return _zero_diag(np.ones_like(w))


def grad_action(w: FloatArray, p: FunctionalParams) -> FloatArray:
    """G = -alpha G_T + beta G_dens + gamma G_deg [+ eta G_smooth] [- mu G_mass] (§1.5)."""
    validate_weight_matrix(w)
    if not isinstance(p, FunctionalParams):
        raise TypeError("p debe ser FunctionalParams")
    g = (
        -p.alpha * grad_triangles(w)
        + p.beta * grad_density(w)
        + p.gamma * grad_degree_irregularity(w)
    )
    if p.eta != 0.0:
        g = g + p.eta * grad_smoothness(w)
    if p.mu != 0.0:
        g = g - p.mu * grad_mass(w)
    return _zero_diag(g)


def lipschitz_bound(n: int, p: FunctionalParams) -> float:
    """L = 2 beta + 2 (gamma + |alpha|)(N-2) (§1.6, D-3); lanza ValueError si eta != 0."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError(f"n debe ser entero >= 2, recibido {n!r}")
    if not isinstance(p, FunctionalParams):
        raise TypeError("p debe ser FunctionalParams")
    if p.eta != 0.0:
        raise ValueError("cota de Lipschitz no definida con eta != 0: use dt_fixed (D-3)")
    return float(2.0 * p.beta + 2.0 * (p.gamma + abs(p.alpha)) * (n - 2))


def theta_gradient(g: FloatArray, w: FloatArray) -> FloatArray:
    """dS/dtheta_ab = G_ab W_ab (1 - W_ab) para W = sigmoide(theta) (D-13, M§16)."""
    validate_weight_matrix(w)
    g = np.asarray(g, dtype=np.float64)
    if g.shape != w.shape:
        raise ValueError(f"g y w deben tener la misma forma: {g.shape} vs {w.shape}")
    return _zero_diag(g * w * (1.0 - w))
