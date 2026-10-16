"""Wang-Landau 1/t y produccion multicanonica sobre 2-ordenes (L-3c-EQ, Rev. 3 §16).

Estado: permutacion v (u = identidad). Energia: accion BD 2D sin suavizar (entera, S = N mod 2).
Los numeros aleatorios se generan en lotes de BATCH pasos ligados al indice absoluto, de modo que
el resultado no depende de como se partan los tramos (siempre que sean multiplos de BATCH).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

from omega.causality.estimators import bd_action_2d
from omega.causality.mcmc import order_from_v

__all__ = ["BATCH", "FLAT_EVERY", "WLState", "new_state", "wl_run", "prod_run", "exact_dos", "normalized_lng"]

BATCH = 100_000
FLAT_EVERY = 10_000
FLATNESS = 0.8


@dataclass
class WLState:
    n: int
    smin: int
    smax: int
    v: list[int]
    s: int
    lng: list[float]
    hist: list[float]
    lnf: float
    t: int
    one_over_t: bool
    rng_state: dict[str, Any]
    prod_visits: list[int]
    prod_flags: list[int]
    prod_t: int

    def nbins(self) -> int:
        return (self.smax - self.smin) // 2 + 1

    def to_json(self) -> dict[str, Any]:
        return self.__dict__.copy()

    @staticmethod
    def from_json(d: dict[str, Any]) -> WLState:
        return WLState(**d)


def _bin(st: WLState, s: int) -> int:
    return (s - st.smin) // 2


def new_state(n: int, smin: int, smax: int, rng: np.random.Generator) -> WLState:
    if (smin - n) % 2 or (smax - n) % 2:
        raise ValueError("smin y smax deben tener la paridad de N")
    v = rng.permutation(n).astype(np.int64)
    s = bd_action_2d(order_from_v(v))
    if not smin <= s <= smax:
        raise ValueError("estado inicial fuera de la ventana")
    nb = (smax - smin) // 2 + 1
    return WLState(n, smin, smax, [int(x) for x in v], int(s), [0.0] * nb, [0.0] * nb, 1.0, 0, False,
                   dict(rng.bit_generator.state), [0] * nb, [0] * nb, 0)


def _batches(st: WLState, steps: int) -> tuple[np.random.Generator, int]:
    if steps % BATCH:
        raise ValueError("los tramos deben ser multiplos de BATCH")
    rng = np.random.Generator(np.random.PCG64())
    rng.bit_generator.state = st.rng_state
    return rng, steps // BATCH


def wl_run(st: WLState, steps: int) -> WLState:
    rng, nbatch = _batches(st, steps)
    v = np.array(st.v, dtype=np.int64)
    lng = np.array(st.lng)
    hist = np.array(st.hist)
    n, s = st.n, st.s
    for _ in range(nbatch):
        a = rng.integers(0, n, size=BATCH)
        b = rng.integers(0, n - 1, size=BATCH)
        b = b + (b >= a)
        u = np.log(rng.random(size=BATCH) + 1e-300)
        for i in range(BATCH):
            ai, bi = a[i], b[i]
            v[ai], v[bi] = v[bi], v[ai]
            s2 = bd_action_2d(order_from_v(v))
            if st.smin <= s2 <= st.smax and u[i] <= lng[(s - st.smin) // 2] - lng[(s2 - st.smin) // 2]:
                s = s2
            else:
                v[ai], v[bi] = v[bi], v[ai]
            k = (s - st.smin) // 2
            lng[k] += st.lnf
            hist[k] += 1
            st.t += 1
            nvis = int(np.count_nonzero(lng))
            if st.one_over_t:
                st.lnf = nvis / st.t
            elif st.t % FLAT_EVERY == 0:
                vis = hist[hist > 0]
                if vis.min() > FLATNESS * vis.mean():
                    st.lnf /= 2
                    hist[:] = 0
                if st.lnf <= nvis / st.t:
                    st.one_over_t = True
                    st.lnf = nvis / st.t
    st.v, st.s, st.lng, st.hist = [int(x) for x in v], int(s), lng.tolist(), hist.tolist()
    st.rng_state = dict(rng.bit_generator.state)
    return st


def prod_run(st: WLState, steps: int, every: int, flag: Callable[[npt.NDArray[np.bool_]], bool]) -> WLState:
    """Multicanonica con ln g congelado; cada `every` pasos registra la bandera en el bin actual."""
    rng, nbatch = _batches(st, steps)
    v = np.array(st.v, dtype=np.int64)
    lng = np.array(st.lng)
    visits = np.array(st.prod_visits)
    flags = np.array(st.prod_flags)
    n, s = st.n, st.s
    for _ in range(nbatch):
        a = rng.integers(0, n, size=BATCH)
        b = rng.integers(0, n - 1, size=BATCH)
        b = b + (b >= a)
        u = np.log(rng.random(size=BATCH) + 1e-300)
        for i in range(BATCH):
            ai, bi = a[i], b[i]
            v[ai], v[bi] = v[bi], v[ai]
            s2 = bd_action_2d(order_from_v(v))
            if st.smin <= s2 <= st.smax and u[i] <= lng[(s - st.smin) // 2] - lng[(s2 - st.smin) // 2]:
                s = s2
            else:
                v[ai], v[bi] = v[bi], v[ai]
            st.prod_t += 1
            if st.prod_t % every == 0:
                k = (s - st.smin) // 2
                visits[k] += 1
                flags[k] += int(flag(order_from_v(v)))
    st.v, st.s = [int(x) for x in v], int(s)
    st.prod_visits, st.prod_flags = visits.tolist(), flags.tolist()
    st.rng_state = dict(rng.bit_generator.state)
    return st


def exact_dos(n: int) -> dict[int, int]:
    """Densidad de estados exacta por enumeracion de las n! permutaciones (solo n chico)."""
    out: dict[int, int] = {}
    for p in itertools.permutations(range(n)):
        s = bd_action_2d(order_from_v(np.array(p, dtype=np.int64)))
        out[s] = out.get(s, 0) + 1
    return out


def normalized_lng(st: WLState, total: float) -> dict[int, float]:
    lng = np.array(st.lng)
    vis = lng > 0
    m = lng[vis].max()
    lz = m + math.log(float(np.sum(np.exp(lng[vis] - m))))
    return {st.smin + 2 * int(k): float(lng[k] - lz + math.log(total)) for k in np.flatnonzero(vis)}
