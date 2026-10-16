"""Metropolis sobre 2-ordenes con peso e^{-beta S} (L-3c-E, Rev. 3 §11.6).

Estado: permutacion v (u = identidad). i precede a j si i < j y v_i < v_j.
Movimiento: transponer v_a y v_b con a != b uniformes. Aceptacion min(1, e^{-beta dS}).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import numpy.typing as npt

from omega.causality.estimators import bd_action_2d

__all__ = ["order_from_v", "ChainResult", "run_chain"]


def order_from_v(v: npt.NDArray[np.int64]) -> npt.NDArray[np.bool_]:
    n = v.size
    idx = np.arange(n)
    rel: npt.NDArray[np.bool_] = (idx[:, None] < idx[None, :]) & (v[:, None] < v[None, :])
    return rel


@dataclass
class ChainResult:
    samples: list[npt.NDArray[np.bool_]] = field(default_factory=list)
    actions: list[int] = field(default_factory=list)
    acceptance: float = 0.0


def run_chain(
    v0: npt.NDArray[np.int64], beta: float, therm_sweeps: int, meas_sweeps: int, every: int, rng: np.random.Generator
) -> ChainResult:
    v = v0.copy()
    n = v.size
    s_cur = bd_action_2d(order_from_v(v))
    res = ChainResult()
    accepted = total = 0
    for sweep in range(therm_sweeps + meas_sweeps):
        for _ in range(n):
            a, b = rng.choice(n, size=2, replace=False)
            v[a], v[b] = v[b], v[a]
            s_new = bd_action_2d(order_from_v(v))
            ds = s_new - s_cur
            total += 1
            if ds <= 0 or rng.random() < math.exp(-beta * ds):
                s_cur = s_new
                accepted += 1
            else:
                v[a], v[b] = v[b], v[a]
        if sweep >= therm_sweeps and (sweep - therm_sweeps + 1) % every == 0:
            res.samples.append(order_from_v(v))
            res.actions.append(s_cur)
    res.acceptance = accepted / total if total else 0.0
    return res
