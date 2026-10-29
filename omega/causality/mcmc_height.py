"""Metropolis sesgado en la altura sobre ordenes con etiquetado natural (R18-1, OMEGA_R18_WL_ALTURA_PREREGISTRO.md).

Medida objetivo de la produccion: pi_W(C) prop. a exp(-W(h(C))) sobre todos los ordenes de N elementos con etiquetado
natural fijo (beta = 0), con h la altura (numero de elementos de la cadena mas larga). W se da en la ventana
h en [HMIN, hmax]; fuera de ella toma el valor del borde mas cercano, de modo que no se excluye ningun estado.

Movimiento: el de ``mcmc_general.sweep`` (par i < j al azar; quitar si es link, agregar si deja el conjunto
transitivo), con el mismo consumo del generador por barrido (i, j y u en lotes de N(N-1)/2). La aceptacion es
min(1, exp(-(W(h') - W(h)))); las propuestas invalidas son rechazos. Con W constante la cadena es identica bit a bit a
``mcmc_general.sweep`` con beta = 0.

Estado: C, B[a, b] = #{m : a < m < b}, depth[x] (cadena mas larga que termina en x, contando elementos), codepth[x]
(cadena mas larga que empieza en x), h = max depth, S4 (accion BD 4D sin suavizar, solo observable, como
``mcmc_general``) y el numero de relaciones.

Altura incremental exacta:
- agregar i < j: h' = max(h, depth[i] + codepth[j]) (la nueva cadena mas larga por el par (i, j));
- quitar el link (i, j): si depth[i] + codepth[j] < h, ninguna cadena de largo h usa el link y h' = h; si no, se
  recalcula depth en {j} U suc(j) sobre el orden propuesto (relajacion de Bellman-Ford desde abajo, exacta en un
  DAG) y h' = max depth'.
Tras aceptar, depth se actualiza en {j} U suc(j) y codepth en {i} U pred(i) (los unicos que pueden cambiar).

Sectores y viajes (§2): sector bajo h <= H_LOW, alto h >= H_HIGH. El sector se sigue en cada movimiento aceptado
(resolucion de propuesta). Un viaje completo es bajo -> alto -> bajo: se cuenta en cada paso alto -> bajo precedido de
algun paso bajo -> alto en la misma cadena.

Wang-Landau (fase A, una cadena): despues de cada propuesta (aceptada, rechazada o invalida) W[h] += ln f y H[h] += 1,
solo si h esta en la ventana (fuera de ella no se actualiza nada). W es la estimacion de ln g(h): el peso muestreado es
exp(-W(h)). Cada barrido se mira la planitud (toda barra >= FLATNESS x media); si es plano, ln f se parte por 2 y H
vuelve a cero. Se para con ln f < lnf_min o al agotar el presupuesto.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import numpy.typing as npt

from omega.causality.mcmc_general import COEF4

__all__ = [
    "FLATNESS",
    "HMIN",
    "H_HIGH",
    "H_LOW",
    "Bias",
    "HeightState",
    "Proposal",
    "accept_prob",
    "apply",
    "codepth_full",
    "depth_full",
    "init_state",
    "propose",
    "sector",
    "sweep",
    "wang_landau",
]

HMIN = 2
H_LOW, H_HIGH = 2, 5
FLATNESS = 0.8

Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]
Float = npt.NDArray[np.float64]


def _c(k: Int) -> Float:
    out: Float = COEF4[np.minimum(k, 4)]
    return out


@dataclass
class Bias:
    """W(h) en la ventana [HMIN, hmax] (w[k] = W(HMIN + k)); fuera de la ventana, el valor del borde."""

    w: Float
    hmax: int

    def __post_init__(self) -> None:
        if self.w.shape != (self.hmax - HMIN + 1,):
            raise ValueError("w debe tener hmax - HMIN + 1 valores")

    def at(self, h: int) -> float:
        return float(self.w[min(max(h, HMIN), self.hmax) - HMIN])

    def full(self, n: int) -> Float:
        """W(h) para h = 0..n (con el recorte en los bordes)."""
        idx = np.clip(np.arange(n + 1), HMIN, self.hmax) - HMIN
        out: Float = self.w[idx]
        return out


@dataclass
class HeightState:
    c: Bool
    b: Int
    depth: Int
    codepth: Int
    h: int
    s4: float
    nrel: int
    last_ext: int = -1  # ultimo sector extremo visitado: -1 ninguno, 0 bajo, 1 alto
    lh_seen: bool = False  # hubo algun paso bajo -> alto
    crossings: int = 0  # pasos entre sectores extremos (en cualquier sentido)
    trips: int = 0  # viajes completos bajo -> alto -> bajo


@dataclass
class Proposal:
    i: int
    j: int
    sign: int  # +1 agrega (i, j), -1 quita el link (i, j)
    ds: float  # Delta S4
    h_new: int
    succ: Int
    pred: Int
    depth_new: Int | None = field(default=None)  # depth del orden propuesto, si ya se calculo


def sector(h: int) -> int:
    """0 bajo, 1 alto, -1 intermedio."""
    return 0 if h <= H_LOW else (1 if h >= H_HIGH else -1)


def depth_full(c: Bool) -> Int:
    n = c.shape[0]
    d = np.ones(n, dtype=np.int64)
    for j in range(n):
        p = np.flatnonzero(c[:j, j])
        if p.size:
            d[j] = int(np.max(d[p])) + 1
    return d


def codepth_full(c: Bool) -> Int:
    n = c.shape[0]
    e = np.ones(n, dtype=np.int64)
    for i in range(n - 1, -1, -1):
        s = np.flatnonzero(c[i, i + 1:])
        if s.size:
            e[i] = int(np.max(e[s + i + 1])) + 1
    return e


def _relax_depth(c: Bool, depth: Int, a: Int, from_below: bool) -> Int:
    """depth exacto en los nodos a (el resto ya es exacto). Arranca de una cota inferior y relaja hasta el punto fijo."""
    d = depth.copy()
    if from_below:
        d[a] = 1
    sub = c[:, a]
    while True:
        new = 1 + (sub * d[:, None]).max(axis=0)
        if np.array_equal(new, d[a]):
            return d
        d[a] = new


def _relax_codepth(c: Bool, codepth: Int, a: Int, from_below: bool) -> Int:
    e = codepth.copy()
    if from_below:
        e[a] = 1
    sub = c[a, :]
    while True:
        new = 1 + (sub * e[None, :]).max(axis=1)
        if np.array_equal(new, e[a]):
            return e
        e[a] = new


def init_state(c: Bool) -> HeightState:
    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError("c debe ser cuadrada")
    if np.any(np.tril(c)):
        raise ValueError("el orden debe estar etiquetado naturalmente (triangular superior estricta)")
    cf = c.astype(np.float32)
    b = np.rint(cf @ cf).astype(np.int64)
    depth = depth_full(c)
    s4 = float(c.shape[0] + float(np.sum(_c(b[c]))))
    st = HeightState(c.copy(), b, depth, codepth_full(c), int(np.max(depth)) if c.shape[0] else 0, s4,
                     int(np.count_nonzero(c)))
    _track(st)
    return st


def _track(st: HeightState) -> None:
    sec = sector(st.h)
    if sec < 0:
        return
    if st.last_ext >= 0 and sec != st.last_ext:
        st.crossings += 1
        if sec == 1:
            st.lh_seen = True
        elif st.lh_seen:
            st.trips += 1
    st.last_ext = sec


def propose(st: HeightState, i: int, j: int) -> Proposal | None:
    """Propuesta sobre el par i < j (None si es invalida). No modifica el estado."""
    c, b = st.c, st.b
    succ = np.flatnonzero(c[j])
    pred = np.flatnonzero(c[:, i])
    if c[i, j]:
        if b[i, j] != 0:
            return None  # relacion que no es link: no se puede quitar
        sign = -1
        ds = 1.0
    else:
        if np.any(~c[i, succ]) or np.any(~c[pred, j]):
            return None  # agregarla romperia la transitividad
        sign = 1
        ds = -1.0
    k1 = b[i, succ]
    k2 = b[pred, j]
    ds += float(np.sum(_c(k1 + sign) - _c(k1)) + np.sum(_c(k2 + sign) - _c(k2)))
    through = int(st.depth[i] + st.codepth[j])  # cadena mas larga que usa el par (i, j)
    if sign > 0:
        return Proposal(i, j, 1, ds, max(st.h, through), succ, pred)
    if through < st.h:
        return Proposal(i, j, -1, ds, st.h, succ, pred)
    c[i, j] = False
    d = _relax_depth(c, st.depth, np.concatenate(([j], succ)), from_below=True)
    c[i, j] = True
    return Proposal(i, j, -1, ds, int(np.max(d)), succ, pred, d)


def accept_prob(bias: Bias, h: int, h_new: int) -> float:
    return min(1.0, math.exp(-(bias.at(h_new) - bias.at(h))))


def apply(st: HeightState, p: Proposal) -> None:
    c, b = st.c, st.b
    i, j, sign = p.i, p.j, p.sign
    c[i, j] = sign > 0
    b[i, p.succ] += sign
    b[p.pred, j] += sign
    st.s4 += p.ds
    st.nrel += sign
    a_dep = np.concatenate(([j], p.succ))
    a_cod = np.concatenate(([i], p.pred))
    if p.depth_new is not None:
        st.depth = p.depth_new
    else:
        st.depth = _relax_depth(c, st.depth, a_dep, from_below=sign < 0)
    st.codepth = _relax_codepth(c, st.codepth, a_cod, from_below=sign < 0)
    if p.h_new != st.h:
        st.h = p.h_new
        _track(st)


@dataclass
class WLAcc:
    lnf: float
    hist: Int


def sweep(st: HeightState, bias: Bias, rng: np.random.Generator, wl: WLAcc | None = None,
          check: Callable[[HeightState], None] | None = None) -> tuple[int, int]:
    """N(N-1)/2 propuestas (mismo flujo del generador que mcmc_general.sweep). Devuelve (aceptadas, validas).

    Con ``wl``, cada propuesta actualiza W (in situ en ``bias.w``) y el histograma si h esta en la ventana.
    ``check`` se llama despues de cada propuesta (solo para pruebas)."""
    n = st.c.shape[0]
    m = n * (n - 1) // 2
    i_arr = rng.integers(0, n, size=m)
    j_arr = rng.integers(0, n - 1, size=m)
    j_arr = j_arr + (j_arr >= i_arr)
    u = rng.random(size=m)
    w = bias.w
    nwin = w.size
    acc = val = 0
    for t in range(m):
        i, j = int(i_arr[t]), int(j_arr[t])
        if i > j:
            i, j = j, i
        p = propose(st, i, j)
        if p is not None:
            val += 1
            if p.h_new == st.h:
                dw = 0.0
            else:
                dw = bias.at(p.h_new) - bias.at(st.h)
            if dw <= 0 or u[t] < math.exp(-dw):
                apply(st, p)
                acc += 1
        if wl is not None:
            k = st.h - HMIN
            if 0 <= k < nwin:
                w[k] += wl.lnf
                wl.hist[k] += 1
        if check is not None:
            check(st)
    return acc, val


def wang_landau(st: HeightState, bias: Bias, rng: np.random.Generator, max_sweeps: int, lnf0: float = 1.0,
                lnf_min: float = 1e-4, flatness: float = FLATNESS, check_every: int = 1,
                on_sweep: Callable[[int, HeightState, float], None] | None = None) -> dict[str, Any]:
    """Fase de adaptacion de una sola cadena. Modifica ``bias.w`` in situ y devuelve el registro.

    La planitud se mira cada ``check_every`` barridos (R18-1: 1). ``on_sweep(barrido, estado, ln f)`` se llama al
    final de cada barrido (antes de mirar la planitud)."""
    wl = WLAcc(lnf0, np.zeros(bias.w.size, dtype=np.int64))
    schedule: list[dict[str, Any]] = []
    h_trace: list[int] = []
    sweeps = 0
    acc_tot = val_tot = 0
    while sweeps < max_sweeps and wl.lnf >= lnf_min:
        a, v = sweep(st, bias, rng, wl)
        acc_tot += a
        val_tot += v
        sweeps += 1
        h_trace.append(st.h)
        if on_sweep is not None:
            on_sweep(sweeps, st, wl.lnf)
        if sweeps % check_every:
            continue
        mean = float(wl.hist.mean())
        if mean > 0 and float(np.min(wl.hist)) >= flatness * mean:
            schedule.append({"sweep": sweeps, "lnf_antes": wl.lnf, "hist": wl.hist.tolist()})
            wl.lnf /= 2.0
            wl.hist[:] = 0
    n = st.c.shape[0]
    return {"sweeps_used": sweeps, "lnf_final": wl.lnf, "reached_lnf_min": wl.lnf < lnf_min,
            "halvings": schedule, "hist_final": wl.hist.tolist(), "h_trace": h_trace,
            "accept_rate": acc_tot / max(sweeps * n * (n - 1) // 2, 1),
            "valid_frac": val_tot / max(sweeps * n * (n - 1) // 2, 1)}
