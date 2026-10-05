"""Residuos KKT de S0 y de Omega-B sobre la caja [0,1] (Consejo Omega Rev. 2, R2.2 y R2.5 L-3a).

Convencion. G = grad_action(W, p) es dS/dw_ij por arista i<j. Para el descenso proyectado en [0,1]
las condiciones de primer orden son (tol = umbral de "cero"/"uno"):

* S0:       W <= tol => G >= 0;  W >= 1-tol => G <= 0;  interior => G = 0.
* Omega-B:  igual con G + lam en lugar de G, donde lam = nu/dt es el multiplicador de la restriccion
            sum W = W0 (el mismo que `evolve_fixed_density`): G+lam >= 0 en ceros, G+lam <= 0 en
            unos, G+lam = 0 en el interior.

Solo se consideran los pares i<j. El residuo relativo divide la violacion maxima entre max(1, max|G|).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from omega.config.settings import FunctionalParams
from omega.dynamics.gradient import grad_action
from omega.network.weights import validate_weight_matrix
from omega.types import BoolArray, FloatArray, IntArray

__all__ = ["KKTResiduals", "codegree", "kkt_residuals_s0", "kkt_residuals_fixed_density", "open_p3_count"]


@dataclass(frozen=True, slots=True)
class KKTResiduals:
    """Violaciones KKT (absolutas y relativa) y conteos por estrato; `multiplier` solo en Omega-B."""

    max_violation_zero: float
    max_violation_one: float
    max_violation_interior: float
    max_violation: float
    relative_violation: float
    grad_scale: float
    n_zero: int
    n_one: int
    n_interior: int
    multiplier: float | None = None


def codegree(w: FloatArray) -> FloatArray:
    """c_ij = (W^2)_ij con diagonal 0 (codegrado ponderado)."""
    validate_weight_matrix(w)
    c = np.array(w @ w, dtype=np.float64)
    np.fill_diagonal(c, 0.0)
    return c


def _check_tol(tol: float) -> None:
    if not 0.0 <= tol < 0.5:
        raise ValueError(f"tol debe estar en [0, 0.5), recibido {tol}")


def _strata(w: FloatArray, tol: float) -> tuple[BoolArray, BoolArray, tuple[IntArray, IntArray]]:
    iu = np.triu_indices(w.shape[0], k=1)
    we = w[iu]
    return we <= tol, we >= 1.0 - tol, (iu[0], iu[1])


def _build(
    g: FloatArray, zero: BoolArray, one: BoolArray, shift: float, multiplier: float | None
) -> KKTResiduals:
    interior = ~(zero | one)
    h = g + shift
    v0 = float(max(0.0, np.max(-h[zero]))) if np.any(zero) else 0.0
    v1 = float(max(0.0, np.max(h[one & ~zero]))) if np.any(one & ~zero) else 0.0
    vi = float(np.max(np.abs(h[interior]))) if np.any(interior) else 0.0
    vmax = max(v0, v1, vi)
    scale = max(1.0, float(np.max(np.abs(g)))) if g.size else 1.0
    return KKTResiduals(
        max_violation_zero=v0,
        max_violation_one=v1,
        max_violation_interior=vi,
        max_violation=vmax,
        relative_violation=vmax / scale,
        grad_scale=scale,
        n_zero=int(np.sum(zero)),
        n_one=int(np.sum(one & ~zero)),
        n_interior=int(np.sum(interior)),
        multiplier=multiplier,
    )


def kkt_residuals_s0(w: FloatArray, p: FunctionalParams, tol: float = 1e-6) -> KKTResiduals:
    """Residuos KKT de S0 en [0,1] (descenso proyectado): ceros G>=0, unos G<=0, interior G=0."""
    validate_weight_matrix(w)
    _check_tol(tol)
    g_full = grad_action(w, p)
    zero, one, iu = _strata(w, tol)
    return _build(g_full[iu], zero, one, 0.0, None)


def kkt_residuals_fixed_density(w: FloatArray, p: FunctionalParams, tol: float = 1e-6) -> KKTResiduals:
    """Residuos KKT de Omega-B con multiplicador lam estimado (G+lam>=0 en ceros, <=0 en unos, =0 interior).

    Convencion de lam: si hay entradas interiores, lam = -media(G interior). Si no las hay, lam es el
    que minimiza la violacion maxima: con a = max_ceros(-G) y b = min_unos(-G) (cotas inferior y
    superior del intervalo factible [a, b]) se toma el punto medio (a+b)/2 si hay ceros y unos;
    a si solo hay ceros; b si solo hay unos; 0 si no hay ninguno. El valor se devuelve en `multiplier`.
    """
    validate_weight_matrix(w)
    _check_tol(tol)
    g_full = grad_action(w, p)
    zero, one, iu = _strata(w, tol)
    g = g_full[iu]
    one_only = one & ~zero
    interior = ~(zero | one)
    if np.any(interior):
        lam = -float(np.mean(g[interior]))
    else:
        has0, has1 = bool(np.any(zero)), bool(np.any(one_only))
        a = float(np.max(-g[zero])) if has0 else 0.0
        b = float(np.min(-g[one_only])) if has1 else 0.0
        lam = 0.5 * (a + b) if (has0 and has1) else (a if has0 else (b if has1 else 0.0))
    return _build(g, zero, one, lam, lam)


def open_p3_count(w: FloatArray, thr: float) -> int:
    """Numero de PARES i<j no adyacentes (W<=thr) con al menos un vecino comun en el soporte {W>thr}.

    Cuenta pares (no caminos): cada P3 inducido abierto i-k-j con ij fuera del soporte aporta una vez
    al par (i,j), por muchos que sean los vecinos comunes k. Vale 0 sii el soporte es union disjunta
    de cliques.
    """
    validate_weight_matrix(w)
    s = (w > thr).astype(np.float64)
    c = s @ s
    mask = np.triu(np.ones_like(s, dtype=bool), k=1) & (s == 0.0) & (c > 0.0)
    return int(np.sum(mask))
