"""Metropolis sobre ordenes con etiquetado natural y accion BD 4D suavizada S_eps (ES-C0, OMEGA_ESC_PREREGISTRO.md).

Peso pi(C) prop. a exp(-beta S_eps(C)), S_eps = sqrt(eps) (N - eps sum_{pares relacionados} f4(n, eps)), con n el
numero de elementos estrictamente entre el par (formula congelada de ``bd_action_4d_smeared``). Con eps = 1 es la
accion sin suavizar S = N + sum COEF4[min(n, 4)] (``mcmc_general.action4``).

S_eps = const + sum_{pares relacionados} g[n], con const = sqrt(eps) N y g[n] = -sqrt(eps) eps f4(n, eps)
(eps = 1: g[n] = COEF4[min(n, 4)]). El movimiento, el consumo del generador y el mantenimiento de
B[a, b] = #{m : a < m < b} son los de ``mcmc_general.sweep``; con eps = 1 la cadena es identica bit a bit.

Delta S de quitar el link (i, j): sale el termino g[0] del par (i, j) y n baja en 1 para los pares (i, s), s sucesor
de j, y (p, j), p predecesor de i. Agregar la relacion (i, j) es el movimiento inverso.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from omega.causality.estimators import smearing_weight_4d
from omega.causality.mcmc_general import COEF4

__all__ = ["EpsState", "action_eps", "coef_table", "init_state", "sweep"]

FloatArray = npt.NDArray[np.float64]


def coef_table(n: int, eps: float) -> tuple[FloatArray, float]:
    """(g, const): g[k] = aporte de un par relacionado con k elementos entre (k = 0..n), const = sqrt(eps) N."""
    if eps == 1.0:
        g = np.zeros(n + 1, dtype=np.float64)
        m = min(n + 1, COEF4.size)
        g[:m] = COEF4[:m]
        return g, float(n)
    if not 0 < eps < 1:
        raise ValueError("eps debe estar en (0, 1]")
    k = np.arange(n + 1, dtype=np.float64)
    out: FloatArray = -math.sqrt(eps) * eps * smearing_weight_4d(k, eps)
    return out, math.sqrt(eps) * n


@dataclass
class EpsState:
    c: npt.NDArray[np.bool_]
    b: npt.NDArray[np.int64]
    s: float
    eps: float
    g: FloatArray
    const: float


def action_eps(c: npt.NDArray[np.bool_], b: npt.NDArray[np.int64], g: FloatArray, const: float) -> float:
    """S_eps a partir de C y B (B solo se lee en los pares relacionados)."""
    return float(const + float(np.sum(g[b[c]])))


def init_state(c: npt.NDArray[np.bool_], eps: float) -> EpsState:
    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError("c debe ser cuadrada")
    if np.any(np.tril(c)):
        raise ValueError("el orden debe estar etiquetado naturalmente (triangular superior estricta)")
    n = int(c.shape[0])
    g, const = coef_table(n, eps)
    cf = c.astype(np.float32)
    b = np.rint(cf @ cf).astype(np.int64)
    return EpsState(c.copy(), b, action_eps(c, b, g, const), float(eps), g, const)


def sweep(st: EpsState, beta: float, rng: np.random.Generator) -> tuple[int, int]:
    """N(N-1)/2 propuestas (mismo flujo del generador que mcmc_general.sweep). Devuelve (aceptadas, validas)."""
    c, b, g = st.c, st.b, st.g
    n = c.shape[0]
    m = n * (n - 1) // 2
    i_arr = rng.integers(0, n, size=m)
    j_arr = rng.integers(0, n - 1, size=m)
    j_arr = j_arr + (j_arr >= i_arr)
    u = rng.random(size=m)
    g0 = float(g[0])
    acc = val = 0
    for t in range(m):
        i, j = int(i_arr[t]), int(j_arr[t])
        if i > j:
            i, j = j, i
        succ = np.flatnonzero(c[j])
        pred = np.flatnonzero(c[:, i])
        if c[i, j]:
            if b[i, j] != 0:
                continue  # relacion que no es link: no se puede quitar
            sign = -1
            ds = -g0  # sale el par (i, j) con n = 0
        else:
            if np.any(~c[i, succ]) or np.any(~c[pred, j]):
                continue  # agregarla romperia la transitividad
            sign = 1
            ds = g0
        val += 1
        k1 = b[i, succ]
        k2 = b[pred, j]
        ds += float(np.sum(g[k1 + sign] - g[k1]) + np.sum(g[k2 + sign] - g[k2]))
        # valido tambien con beta < 0 (auditoria A-22): no hay atajo 'ds <= 0 acepta'
        if -beta * ds >= 0 or u[t] < math.exp(-beta * ds):
            c[i, j] = sign > 0
            b[i, succ] += sign
            b[pred, j] += sign
            st.s += ds
            acc += 1
    return acc, val
