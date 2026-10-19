"""Metropolis con peso generico e^{lam S} sobre ordenes (ES-0, docs/OMEGA_ES0_PILOTO.md).

E-a: 2-ordenes. Estado permutacion v (u = identidad); movimiento: transponer v_a y v_b con a != b uniformes.
E-b: ordenes con etiquetado natural fijo (triangular superior estricta). Par i < j al azar: quitar si es link,
agregar si deja el conjunto transitivo; si no se puede, la propuesta cuenta como rechazada (sin evaluar el puntaje).
Aceptacion min(1, e^{lam dS}); el puntaje se recalcula completo en cada propuesta valida. Barrido = N(N-1)/2 propuestas.
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, NamedTuple

import numpy as np
import numpy.typing as npt

from omega.causality.entropy import interval_balanced, pauli_jordan, sj_from_eig, ssee
from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_general import height

__all__ = [
    "GeneralOrderState",
    "Score",
    "SweepCounts",
    "TwoOrderState",
    "init_general",
    "init_two_order",
    "measure",
    "s_sj",
    "sweep_general",
    "sweep_two_order",
]

Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]
Score = Callable[[Bool], tuple[float, int]]


class SweepCounts(NamedTuple):
    proposals: int
    valid: int
    accepted: int


def s_sj(c: Bool) -> tuple[float, int]:
    """Instrumento de Omega-E (E1''): (S_SJ, |R|). Sin intervalo del tamano pedido: (0.0, 0) (convencion ES-0 §1)."""
    n = c.shape[0]
    r = interval_balanced(c, n // 8)
    if r.size == 0:
        return 0.0, 0  # el resultado no depende de la diagonalizacion: se omite
    mu, u = np.linalg.eigh(pauli_jordan(c))
    idt, w = sj_from_eig(mu, u, math.sqrt(n) / (4 * math.pi))
    return ssee(idt, w, r, True, math.sqrt(r.size) / (4 * math.pi)), int(r.size)


def _accept(lam: float, ds: float, u: float) -> bool:
    # vale para lam de cualquier signo (cf. AUDITORIA A-22)
    x = lam * ds
    return x >= 0 or u < math.exp(x)


@dataclass
class TwoOrderState:
    v: Int
    s: float
    n_r: int
    n_eval: int = 0
    t_eval: float = 0.0

    def order(self) -> Bool:
        return order_from_v(self.v)


@dataclass
class GeneralOrderState:
    c: Bool
    s: float
    n_r: int
    n_eval: int = 0
    t_eval: float = 0.0

    def order(self) -> Bool:
        return self.c


def _timed(st: TwoOrderState | GeneralOrderState, score: Score, c: Bool) -> tuple[float, int]:
    t0 = time.perf_counter()
    s, nr = score(c)
    st.t_eval += time.perf_counter() - t0
    st.n_eval += 1
    return float(s), int(nr)


def init_two_order(v: Int, score: Score = s_sj) -> TwoOrderState:
    v = np.asarray(v, dtype=np.int64).copy()
    if not np.array_equal(np.sort(v), np.arange(v.size)):
        raise ValueError("v debe ser una permutacion de 0..N-1")
    s, nr = score(order_from_v(v))
    return TwoOrderState(v, float(s), int(nr))


def init_general(c: Bool, score: Score = s_sj) -> GeneralOrderState:
    c = np.asarray(c, dtype=bool).copy()
    if np.any(np.tril(c)):
        raise ValueError("el orden debe estar etiquetado naturalmente (triangular superior estricta)")
    cf = c.astype(np.float32)
    if np.any(((cf @ cf) > 0) & ~c):
        raise ValueError("la relacion no es transitiva")
    s, nr = score(c)
    return GeneralOrderState(c, float(s), int(nr))


def _pairs(n: int, m: int, rng: np.random.Generator) -> tuple[Int, Int]:
    a = rng.integers(0, n, size=m)
    b = rng.integers(0, n - 1, size=m)
    b = b + (b >= a)
    return a, b


def sweep_two_order(st: TwoOrderState, lam: float, rng: np.random.Generator, score: Score = s_sj) -> SweepCounts:
    """N(N-1)/2 transposiciones propuestas (todas validas)."""
    v = st.v
    n = v.size
    m = n * (n - 1) // 2
    a_arr, b_arr = _pairs(n, m, rng)
    u = rng.random(size=m)
    acc = 0
    for t in range(m):
        a, b = int(a_arr[t]), int(b_arr[t])
        v[a], v[b] = v[b], v[a]
        s_new, nr_new = _timed(st, score, order_from_v(v))
        if _accept(lam, s_new - st.s, float(u[t])):
            st.s, st.n_r = s_new, nr_new
            acc += 1
        else:
            v[a], v[b] = v[b], v[a]
    return SweepCounts(m, m, acc)


def sweep_general(st: GeneralOrderState, lam: float, rng: np.random.Generator, score: Score = s_sj) -> SweepCounts:
    """N(N-1)/2 propuestas sobre pares i < j; las invalidas no evaluan el puntaje."""
    c = st.c
    n = c.shape[0]
    m = n * (n - 1) // 2
    i_arr, j_arr = _pairs(n, m, rng)
    u = rng.random(size=m)
    acc = val = 0
    for t in range(m):
        i, j = int(i_arr[t]), int(j_arr[t])
        if i > j:
            i, j = j, i
        if c[i, j]:
            if np.any(c[i] & c[:, j]):
                continue  # relacion no-link: no se puede quitar
        elif np.any(c[j] & ~c[i]) or np.any(c[:, i] & ~c[:, j]):
            continue  # agregarla romperia la transitividad
        val += 1
        new = not bool(c[i, j])
        c[i, j] = new
        s_new, nr_new = _timed(st, score, c)
        if _accept(lam, s_new - st.s, float(u[t])):
            st.s, st.n_r = s_new, nr_new
            acc += 1
        else:
            c[i, j] = not new
    return SweepCounts(m, val, acc)


def measure(st: TwoOrderState | GeneralOrderState) -> dict[str, Any]:
    """Traza por barrido: S, |R|, si hay region, fraccion de pares relacionados, links y altura."""
    c = st.order()
    n = c.shape[0]
    npairs = n * (n - 1) // 2
    cf = c.astype(np.float32)
    links = int(np.sum(c & ((cf @ cf) < 0.5)))
    return {
        "S": float(st.s),
        "nR": int(st.n_r),
        "has_region": bool(st.n_r > 0),
        "rel_frac": float(c.sum()) / npairs if npairs else 0.0,
        "links": links,
        "height": height(c),
    }
