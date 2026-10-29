"""Candidato S_SJ^max de R17-1 (docs/OMEGA_CONSEJO_REV17_ES0.md §5.1 y §6), separado del E1' firmado.

- Region: el intervalo I[a, b] (con extremos, c[a, b] verdadero) de mayor |I|; empate: se maximiza
  min(|pasado(a)|, |futuro(b)|) (sumas de columna y de fila, como `interval_balanced`); empate restante: menor (a, b)
  lexicografico.
- Soporte: todos los ordenes de N elementos menos la anticadena (unico orden sin relaciones). Evaluar la anticadena
  levanta `OutsideSupport`; los muestreadores de este modulo tratan una propuesta hacia ella como rechazada.
- Truncamiento global sqrt(N)/(4 pi) como `mcmc_sj.s_sj`; truncamiento de la region sqrt(|R|)/(4 pi) con un minimo de
  un modo (variante (ii) firmada en §6): si ningun modo pasa el umbral se conserva el de mayor |nu|.

Este modulo no cambia `entropy.ssee` ni `mcmc_sj`: los barridos replican la logica de `mcmc_sj` (mismo consumo del
generador) y solo agregan el predicado de soporte.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt

from omega.causality.entropy import CArr, pauli_jordan, sj_from_eig
from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_sj import (
    GeneralOrderState,
    Score,
    SweepCounts,
    TwoOrderState,
    _accept,
    _pairs,
    _timed,
)

__all__ = [
    "OutsideSupport",
    "has_relation",
    "interval_max",
    "region_shape",
    "s_sj_max",
    "s_sj_max_info",
    "ssee_min1",
    "sweep_general_max",
    "sweep_two_order_max",
]

Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]
Support = Callable[[Bool], bool]
NULL_TOL = 1e-12  # |nu| por debajo de esto es numericamente cero: el modo esta en el nucleo de iDelta_R


class OutsideSupport(ValueError):
    """El orden no pertenece al soporte de S_SJ^max (la anticadena)."""


def has_relation(c: Bool) -> bool:
    """Predicado de soporte: el orden tiene al menos una relacion (no es la anticadena)."""
    return bool(np.any(c))


def interval_max(c: Bool) -> Int:
    """Intervalo I[a, b] (con extremos) de mayor |I|; empates como E1' y luego menor (a, b). Vacio solo en la anticadena."""
    if not np.any(c):
        return np.zeros(0, dtype=np.int64)
    ci = c.astype(np.int64)
    size = ci @ ci + 2
    best = int(np.max(size[c]))
    ok = c & (size == best)
    past = c.sum(axis=0)
    fut = c.sum(axis=1)
    score = np.where(ok, np.minimum(past[:, None], fut[None, :]), -1)
    a, b = np.unravel_index(int(np.argmax(score)), score.shape)  # primer maximo en orden por filas = menor (a, b)
    members = np.flatnonzero(c[a] & c[:, b])
    out: Int = np.sort(np.concatenate([[a, b], members])).astype(np.int64)
    return out


def ssee_min1(i_delta: CArr, w: CArr, region: Int, cutoff: float | None = None) -> tuple[float, int, bool]:
    """Como `entropy.ssee(..., truncate=True, cutoff)` pero con un minimo de un modo (variante (ii)).

    Devuelve (S, modos conservados, si se aplico la regla del minimo). Si el modo de mayor |nu| es numericamente
    nulo (|nu| < NULL_TOL), iDelta_R no tiene imagen: la suma es vacia y S = 0 con 0 modos (lo que da la formula).
    """
    a = i_delta[np.ix_(region, region)]
    b = w[np.ix_(region, region)]
    nu, v = np.linalg.eigh(a)
    thr = cutoff if cutoff is not None else math.sqrt(region.size) / (4 * math.pi)
    keep = np.abs(nu) >= thr
    fired = False
    if not np.any(keep):
        fired = True
        k = int(np.argmax(np.abs(nu)))
        if abs(float(nu[k])) < NULL_TOL:
            return 0.0, 0, True
        keep = np.zeros(nu.size, dtype=bool)
        keep[k] = True
    vk = v[:, keep]
    bk = vk.conj().T @ b @ vk
    lam = np.linalg.eigvals(bk / nu[keep][:, None]).real
    lam = lam[np.abs(lam) > 1e-12]
    return float(np.sum(lam * np.log(np.abs(lam)))), int(keep.sum()), fired


def region_shape(c: Bool, region: Int) -> tuple[int, float]:
    """(altura = cadena mas larga en R, fraccion de pares de R relacionados)."""
    sub = c[np.ix_(region, region)]
    m = region.size
    order = np.argsort(sub.sum(axis=0), kind="stable")  # por numero de antecesores: extension lineal
    sub = sub[np.ix_(order, order)]
    depth = np.ones(m, dtype=np.int64)
    for j in range(m):
        pred = sub[:j, j]
        if np.any(pred):
            depth[j] = int(np.max(depth[:j][pred])) + 1
    npairs = m * (m - 1) // 2
    return int(np.max(depth)) if m else 0, float(sub.sum()) / npairs if npairs else 0.0


def _evaluate(c: Bool) -> tuple[float, Int, int, bool]:
    n = c.shape[0]
    r = interval_max(c)
    if r.size == 0:
        raise OutsideSupport("la anticadena esta fuera del soporte de S_SJ^max")
    mu, u = np.linalg.eigh(pauli_jordan(c))
    idt, w = sj_from_eig(mu, u, math.sqrt(n) / (4 * math.pi))
    s, kept, fired = ssee_min1(idt, w, r, math.sqrt(r.size) / (4 * math.pi))
    return s, r, kept, fired


def s_sj_max(c: Bool) -> tuple[float, int]:
    """Puntaje compatible con `mcmc_sj.Score`: (S_SJ^max, |R|). Levanta `OutsideSupport` en la anticadena."""
    s, r, _, _ = _evaluate(c)
    return s, int(r.size)


def s_sj_max_info(c: Bool) -> dict[str, Any]:
    """Registro de §5.1: S, |R|, modos conservados, regla del minimo, altura y fraccion de relaciones en R."""
    s, r, kept, fired = _evaluate(c)
    h, frac = region_shape(c, r)
    return {"S": s, "nR": int(r.size), "modes": kept, "min1": fired, "height": h, "rel_frac_R": frac}


def sweep_two_order_max(
    st: TwoOrderState,
    lam: float,
    rng: np.random.Generator,
    score: Score = s_sj_max,
    in_support: Support = has_relation,
) -> SweepCounts:
    """Como `mcmc_sj.sweep_two_order` (mismos sorteos); una transposicion que sale del soporte se rechaza sin
    evaluar el puntaje y no cuenta como valida."""
    v = st.v
    n = v.size
    m = n * (n - 1) // 2
    a_arr, b_arr = _pairs(n, m, rng)
    u = rng.random(size=m)
    acc = val = 0
    for t in range(m):
        a, b = int(a_arr[t]), int(b_arr[t])
        v[a], v[b] = v[b], v[a]
        c = order_from_v(v)
        if not in_support(c):
            v[a], v[b] = v[b], v[a]
            continue
        val += 1
        s_new, nr_new = _timed(st, score, c)
        if _accept(lam, s_new - st.s, float(u[t])):
            st.s, st.n_r = s_new, nr_new
            acc += 1
        else:
            v[a], v[b] = v[b], v[a]
    return SweepCounts(m, val, acc)


def sweep_general_max(
    st: GeneralOrderState,
    lam: float,
    rng: np.random.Generator,
    score: Score = s_sj_max,
    in_support: Support = has_relation,
) -> SweepCounts:
    """Como `mcmc_sj.sweep_general` (mismos sorteos); un movimiento valido que sale del soporte se rechaza sin
    evaluar el puntaje y no cuenta como valido."""
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
                continue
        elif np.any(c[j] & ~c[i]) or np.any(c[:, i] & ~c[:, j]):
            continue
        new = not bool(c[i, j])
        c[i, j] = new
        if not in_support(c):
            c[i, j] = not new
            continue
        val += 1
        s_new, nr_new = _timed(st, score, c)
        if _accept(lam, s_new - st.s, float(u[t])):
            st.s, st.n_r = s_new, nr_new
            acc += 1
        else:
            c[i, j] = not new
    return SweepCounts(m, val, acc)
