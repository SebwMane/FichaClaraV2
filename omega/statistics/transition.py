"""Observables de transicion: susceptibilidad, Binder U4, bimodalidad (DESIGN §1.8, §1.14).

lambda es alpha_hat o Theta_hat (temperatura estadistica, no tiempo). En corridas de
gradiente la varianza entre semillas es `seed_variance`, no termica.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from omega.contracts import TransitionSummary
from omega.types import FloatArray

__all__ = [
    "fluctuation",
    "susceptibility",
    "binder_cumulant",
    "bimodality_coefficient",
    "order_parameter_histogram",
    "transition_summary",
]


def _as1d(x: FloatArray, minimum: int = 2) -> FloatArray:
    a = np.asarray(x, dtype=np.float64)
    if a.ndim != 1 or a.shape[0] < minimum:
        raise ValueError(f"se requiere serie 1-D con al menos {minimum} datos")
    return a


def fluctuation(x: FloatArray) -> float:
    """Varianza poblacional <dx^2> de la serie."""
    a = _as1d(x, 1)
    return float(np.mean((a - a.mean()) ** 2))


def susceptibility(m: FloatArray, n_edges: int) -> float:
    """chi_m = M * Var(m), con M = numero de aristas."""
    if n_edges < 1:
        raise ValueError("n_edges debe ser >= 1")
    return float(n_edges * fluctuation(m))


def binder_cumulant(m: FloatArray) -> float:
    """U4 = 1 - <dm^4> / (3 <dm^2>^2) con dm = m - <m>; nan si la varianza es 0."""
    a = _as1d(m, 1)
    d = a - a.mean()
    m2 = float(np.mean(d**2))
    if m2 <= 0.0:
        return float("nan")
    return float(1.0 - np.mean(d**4) / (3.0 * m2 * m2))


def bimodality_coefficient(m: FloatArray) -> float:
    """b = (g^2 + 1) / kappa (Sarle; kappa = curtosis no centrada); b > 5/9 sugiere bimodal. nan si var=0."""
    a = _as1d(m, 1)
    d = a - a.mean()
    m2 = float(np.mean(d**2))
    if m2 <= 0.0:
        return float("nan")
    g = float(np.mean(d**3)) / m2**1.5
    kappa = float(np.mean(d**4)) / (m2 * m2)
    return float((g * g + 1.0) / kappa)


def order_parameter_histogram(m: FloatArray, bins: int = 40) -> tuple[FloatArray, FloatArray]:
    """(densidad, bordes) del histograma de m."""
    a = _as1d(m, 1)
    if bins < 1:
        raise ValueError("bins debe ser >= 1")
    dens, edges = np.histogram(a, bins=bins, density=True)
    return np.asarray(dens, dtype=np.float64), np.asarray(edges, dtype=np.float64)


def transition_summary(samples: Mapping[float, FloatArray], n_edges: int) -> TransitionSummary:
    """Resumen sobre lambda ordenado: <m>, d<m>/dlambda (diferencias centrales), chi, U4, b, lambda*=argmax chi."""
    if len(samples) < 2:
        raise ValueError("se requieren al menos 2 valores de lambda")
    lams = np.array(sorted(float(k) for k in samples), dtype=np.float64)
    series = [_as1d(samples[k]) for k in sorted(samples, key=float)]
    mean = np.array([s.mean() for s in series], dtype=np.float64)
    deriv = np.asarray(np.gradient(mean, lams), dtype=np.float64)
    chi = np.array([susceptibility(s, n_edges) for s in series], dtype=np.float64)
    binder = np.array([binder_cumulant(s) for s in series], dtype=np.float64)
    bim = np.array([bimodality_coefficient(s) for s in series], dtype=np.float64)
    return TransitionSummary(
        lambdas=lams,
        mean=mean,
        derivative=deriv,
        susceptibility=chi,
        binder=binder,
        bimodality=bim,
        peak_lambda=float(lams[int(np.argmax(chi))]),
    )
