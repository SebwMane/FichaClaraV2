"""Metropolis sobre ordenes con etiquetado natural y accion BD 4D sin suavizar (piloto P4, Rev. 3 §21).

Movimiento simetrico: par i < j al azar; quitar si es link, agregar si deja el conjunto transitivo.
Se mantiene B[a, b] = #{m : a < m < b en el orden} para todos los pares, actualizado en O(N).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

__all__ = ["COEF4", "GeneralState", "init_state", "action4", "sweep", "height"]

COEF4 = np.array([-1.0, 9.0, -16.0, 8.0, 0.0])  # aporte de un par con k = 0..3 elementos entre (k >= 4: 0)


def _c(k: npt.NDArray[np.int64]) -> npt.NDArray[np.float64]:
    out: npt.NDArray[np.float64] = COEF4[np.minimum(k, 4)]
    return out


@dataclass
class GeneralState:
    c: npt.NDArray[np.bool_]
    b: npt.NDArray[np.int64]
    s: float


def action4(c: npt.NDArray[np.bool_], b: npt.NDArray[np.int64]) -> float:
    return float(c.shape[0] + float(np.sum(_c(b[c]))))


def init_state(c: npt.NDArray[np.bool_]) -> GeneralState:
    if np.any(np.tril(c)):
        raise ValueError("el orden debe estar etiquetado naturalmente (triangular superior estricta)")
    cf = c.astype(np.float32)
    b = np.rint(cf @ cf).astype(np.int64)
    return GeneralState(c.copy(), b, action4(c, b))


def sweep(st: GeneralState, beta: float, rng: np.random.Generator) -> tuple[int, int]:
    """N(N-1)/2 propuestas. Devuelve (aceptadas, validas)."""
    c, b = st.c, st.b
    n = c.shape[0]
    m = n * (n - 1) // 2
    i_arr = rng.integers(0, n, size=m)
    j_arr = rng.integers(0, n - 1, size=m)
    j_arr = j_arr + (j_arr >= i_arr)
    u = rng.random(size=m)
    acc = val = 0
    for t in range(m):
        i, j = int(i_arr[t]), int(j_arr[t])
        if i > j:
            i, j = j, i
        succ = np.flatnonzero(c[j])
        pred = np.flatnonzero(c[:, i])
        if c[i, j]:
            if b[i, j] != 0:
                continue  # relacion no-link: no se puede quitar
            sign = -1
            ds = 1.0  # sale un link: -(-1)
        else:
            if np.any(~c[i, succ]) or np.any(~c[pred, j]):
                continue  # agregarla romperia la transitividad
            sign = 1
            ds = -1.0
        val += 1
        k1 = b[i, succ]
        k2 = b[pred, j]
        ds += float(np.sum(_c(k1 + sign) - _c(k1)) + np.sum(_c(k2 + sign) - _c(k2)))
        # AUDITORIA A-22: el atajo 'ds <= 0 acepta' solo vale para beta >= 0
        if -beta * ds >= 0 or u[t] < math.exp(-beta * ds):
            c[i, j] = sign > 0
            b[i, succ] += sign
            b[pred, j] += sign
            st.s += ds
            acc += 1
    return acc, val


def height(c: npt.NDArray[np.bool_]) -> int:
    """Cadena mas larga (numero de elementos) en un orden etiquetado naturalmente."""
    n = c.shape[0]
    h = np.ones(n, dtype=np.int64)
    for j in range(n):
        p = np.flatnonzero(c[:j, j])
        if p.size:
            h[j] = int(h[p].max()) + 1
    return int(h.max()) if n else 0
