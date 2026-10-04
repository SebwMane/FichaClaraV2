"""Funcional S0 literal (M§13-M§15; ANALYSIS §1.1, D-14, D-15): sumas i<j e i<j<k."""

from __future__ import annotations

import numpy as np

from omega.config.settings import FunctionalParams
from omega.network.weights import validate_weight_matrix
from omega.types import FloatArray

__all__ = ["triangles", "density", "degree_irregularity", "smoothness", "mass", "action"]


def triangles(w: FloatArray) -> float:
    """T = sum_{i<j<k} W_ij W_jk W_ki = tr(W^3)/6 (M§13)."""
    validate_weight_matrix(w)
    return float(np.sum((w @ w) * w) / 6.0)


def density(w: FloatArray) -> float:
    """S_dens = sum_{i<j} W_ij^2 (M§13)."""
    validate_weight_matrix(w)
    return float(np.sum(w * w) / 2.0)


def degree_irregularity(w: FloatArray) -> float:
    """S_deg = sum_i (k_i - k_medio)^2 con k_i = sum_j W_ij (M§13, D-14)."""
    validate_weight_matrix(w)
    k = w.sum(axis=1)
    return float(np.sum((k - k.mean()) ** 2))


def smoothness(w: FloatArray) -> float:
    """S_smooth sin eta: sum_{i,j} W_ij C_ij = 2 sum_i k_i s_i - 2 tr(W^3), s=diag(W^2) (M§15)."""
    validate_weight_matrix(w)
    w2 = w @ w
    k = w.sum(axis=1)
    s = np.diag(w2)
    return float(2.0 * np.sum(k * s) - 2.0 * np.sum(w2 * w))


def mass(w: FloatArray) -> float:
    """S_mass sin mu: sum_{i<j} W_ij (M§15, D-30)."""
    validate_weight_matrix(w)
    return float(w.sum() / 2.0)


def action(w: FloatArray, p: FunctionalParams) -> float:
    """S0 = -alpha T + beta S_dens + gamma S_deg + eta S_smooth - mu S_mass (M§13-M§15)."""
    validate_weight_matrix(w)
    if not isinstance(p, FunctionalParams):
        raise TypeError("p debe ser FunctionalParams")
    total = -p.alpha * triangles(w) + p.beta * density(w) + p.gamma * degree_irregularity(w)
    if p.eta != 0.0:
        total += p.eta * smoothness(w)
    if p.mu != 0.0:
        total -= p.mu * mass(w)
    return float(total)
