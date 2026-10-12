"""Valores esperados exactos para 2-ordenes aleatorios (auditoria C0, enmienda C-A2).

Dos elementos de un 2-orden aleatorio estan relacionados con probabilidad 1/2; dado que lo estan,
las diferencias en los dos ordenes son independientes con densidad 2(1 - x) en [0, 1]. Por tanto

    E[N_k] = C(N, 2) * (1/2) * integral_0^1 integral_0^1 4(1-a)(1-b) Binom(k; N-2, a b) da db.

Con esto, E[S_BD] = N - 2E[N0] + 4E[N1] - 2E[N2] = 1 para todo N (verificado numericamente).
"""

from __future__ import annotations

import math

from scipy import integrate
from scipy.stats import binom

__all__ = ["expected_abundance_2order", "expected_bd_action_2order"]


def expected_abundance_2order(k: int, n: int) -> float:
    if n < 2 or k < 0 or k > n - 2:
        raise ValueError("se requiere n >= 2 y 0 <= k <= n-2")

    def f(b: float, a: float) -> float:
        return float(4 * (1 - a) * (1 - b) * binom.pmf(k, n - 2, a * b))

    val, _ = integrate.dblquad(f, 0, 1, 0, 1, epsabs=1e-11, epsrel=1e-10)
    return math.comb(n, 2) * 0.5 * float(val)


def expected_bd_action_2order(n: int) -> float:
    e = [expected_abundance_2order(k, n) for k in range(3)]
    return n - 2 * e[0] + 4 * e[1] - 2 * e[2]
