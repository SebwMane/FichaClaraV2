"""R17-3: evaluacion exacta y mas rapida del S_SJ original de E1' (`mcmc_sj.s_sj`), docs/OMEGA_CONSEJO_REV17_ES0.md §5.2.

Mismo valor que `s_sj` (criterio: |dS| <= 1e-9 e igual |R|), con tres cambios de algebra exacta:

1. Estructura real. iDelta = i A con A = K/2, K = C - C^T entera antisimetrica. (iDelta)^2 = A^T A, asi que
   |mu| = sqrt(eig(K^T K))/2, con autovectores reales. M = K^T K es entera (exacta en float64). Con U_k los
   autovectores de M con |mu| >= corte (modos conservados, por pares +-mu):
       |iDelta_t| = U_k diag(mu) U_k^T,   iDelta_t = i A U_k U_k^T,   W_t = (|iDelta_t| + iDelta_t)/2,
   que es exactamente (iDelta truncado, W truncado) de `sj_from_eig` (W = parte positiva = (|H| + H)/2).
   Se reemplaza el eigh complejo N x N por un eigh real simetrico (3 a 7 veces mas rapido).
2. Solo los bloques de la region: iDelta_R y W_R (|R| x |R|) en vez de los productos N x N.
3. Evaluador con estado (muestreadores): M = K^T K y B = C C (conteos de intervalos) se actualizan en O(N^2)
   por movimiento (relacion agregada/quitada en E-b, transposicion en E-a), con aritmetica entera exacta;
   la region se elige desde B con la misma regla y el mismo desempate que `interval_balanced`.

Umbrales: el truncamiento global (|mu| >= sqrt(N)/(4 pi)) y el de la region (|nu| >= sqrt(|R|)/(4 pi)) hacen
discontinuo a S. Si algun autovalor queda a menos de `guard` (1e-8) de cualquiera de los dos umbrales, el valor
se calcula con `s_sj` original (respaldo), de modo que los cruces de umbral se deciden igual que en la referencia.
Fuera de esa banda ambos calculos (error ~1e-12) deciden lo mismo. El filtro |lambda| > 1e-12 de `ssee` no se
protege: un lambda en ese borde aporta menos de 3e-11.

Lo que sigue siendo O(N^3): la diagonalizacion de M (tridiagonalizacion). No se adopto la actualizacion de bajo
rango de la diagonalizacion (ver tools/r17_3_bench.py): con aceptacion 0.5-1 el producto U <- U U~ al aceptar ya
cuesta lo mismo que el eigh real para N <= 512.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import numpy.typing as npt

from omega.causality.entropy import interval_balanced
from omega.causality.mcmc_sj import SweepCounts, s_sj

try:  # el eigh por subconjunto (dsyevr) solo gana para N grande
    import scipy.linalg as _sl
except ImportError:  # pragma: no cover
    _sl = None

__all__ = [
    "GUARD",
    "SJEvaluator",
    "SJResult",
    "evaluate",
    "s_sj_fast",
    "sweep_general_fast",
    "sweep_two_order_fast",
]

Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]
F64 = npt.NDArray[np.float64]
C128 = npt.NDArray[np.complex128]

GUARD = 1e-8  # banda alrededor de los umbrales donde se usa el calculo original
SUBSET_MIN_N = 400  # desde aqui se piden solo los autovalores >= (0.9 corte) (dsyevr)
SUBSET_FRAC = 0.9
TOL = 0.05  # tolerancia de interval_balanced en s_sj


class SJResult(NamedTuple):
    s: float
    n_r: int
    margin_global: float  # min | |mu| - corte global | entre los autovalores calculados (inf sin region)
    margin_region: float  # min | |nu| - corte de la region | (inf sin region)
    n_kept: int  # modos globales conservados (cuenta +-mu por separado)
    fallback: bool  # True si el valor viene de s_sj original


_EMPTY = SJResult(0.0, 0, math.inf, math.inf, 0, False)


def _kept_modes(m: F64, n: int) -> tuple[F64, F64, float]:
    """Autopares de M = K^T K con mu = sqrt(lambda)/2 >= sqrt(N)/(4 pi); devuelve (mu_k, U_k, margen)."""
    cut = math.sqrt(n) / (4 * math.pi)
    if n >= SUBSET_MIN_N and _sl is not None:
        lo = (2 * SUBSET_FRAC * cut) ** 2
        lam, u = _sl.eigh(m, subset_by_value=(lo, np.inf), driver="evr", check_finite=False)
        floor = (1 - SUBSET_FRAC) * cut  # los no devueltos tienen mu < 0.9 corte
    else:
        lam, u = np.linalg.eigh(m)
        floor = math.inf
    mu = np.sqrt(np.maximum(lam, 0.0)) / 2
    margin = min(floor, float(np.min(np.abs(mu - cut)))) if mu.size else floor
    keep = mu >= cut
    return mu[keep], u[:, keep], margin


def _ssee_blocks(idt_r: C128, w_r: C128, cut_r: float) -> tuple[float, float]:
    """Mismo algebra que entropy.ssee(truncate=True) sobre bloques ya restringidos; devuelve (S, margen)."""
    nu, v = np.linalg.eigh(idt_r)
    margin = float(np.min(np.abs(np.abs(nu) - cut_r)))
    keep = np.abs(nu) >= cut_r
    if not np.any(keep):
        return 0.0, margin
    vk = v[:, keep]
    bk = vk.conj().T @ w_r @ vk
    lam = np.linalg.eigvals(bk / nu[keep][:, None]).real
    lam = lam[np.abs(lam) > 1e-12]
    return float(np.sum(lam * np.log(np.abs(lam)))), margin


def _s_region(c: Bool, k: F64, m: F64, region: Int, guard: float) -> SJResult:
    """S_SJ dado K = C - C^T, M = K^T K y la region (no vacia)."""
    n = c.shape[0]
    mu, uk, mg = _kept_modes(m, n)
    ur = uk[region]
    pr = (ur * mu) @ ur.T  # |iDelta_t|_R
    au = (k[region] @ uk) * 0.5  # (A U_k)_R
    idt_r: C128 = (1j * (au @ ur.T)).astype(np.complex128)  # (i A U_k U_k^T)_R
    w_r: C128 = (0.5 * (pr + idt_r)).astype(np.complex128)
    s, mr = _ssee_blocks(idt_r, w_r, math.sqrt(region.size) / (4 * math.pi))
    if mg < guard or mr < guard:
        s0, nr0 = s_sj(c)
        return SJResult(float(s0), int(nr0), mg, mr, int(mu.size), True)
    return SJResult(s, int(region.size), mg, mr, int(mu.size), False)


def evaluate(c: Bool, guard: float = GUARD) -> SJResult:
    """Evaluacion sin estado con diagnosticos (margenes, respaldo)."""
    c = np.asarray(c, dtype=bool)
    region = interval_balanced(c, c.shape[0] // 8)
    if region.size == 0:
        return _EMPTY
    k = c.astype(np.float64)
    k = k - k.T
    return _s_region(c, k, k.T @ k, region, guard)


def s_sj_fast(c: Bool) -> tuple[float, int]:
    """Reemplazo directo de mcmc_sj.s_sj (misma firma y mismo valor)."""
    r = evaluate(c)
    return r.s, r.n_r


def _select_region(c: Bool, size: Int, target: int) -> Int:
    """interval_balanced con los tamanos |I| = B + 2 ya calculados (misma regla y mismo desempate)."""
    ok = c & (size >= math.floor(target * (1 - TOL))) & (size <= math.ceil(target * (1 + TOL)))
    if not np.any(ok):
        return np.zeros(0, dtype=np.int64)
    past = c.sum(axis=0)
    fut = c.sum(axis=1)
    score = np.where(ok, np.minimum(past[:, None], fut[None, :]), -1)
    a, b = np.unravel_index(int(np.argmax(score)), score.shape)
    members = np.flatnonzero(c[a] & c[:, b])
    out: Int = np.sort(np.concatenate([[a, b], members])).astype(np.int64)
    return out


@dataclass
class _Undo:
    s_idx: Int
    rows: Bool
    cols: Bool
    m_old: F64
    b_old: F64
    swap: tuple[int, int] | None
    result: SJResult


class SJEvaluator:
    """Evaluador con estado: propose_toggle / propose_swap devuelven S del estado propuesto; luego commit o rollback.

    E-b: `propose_toggle(i, j)` invierte la relacion (i, j) (el llamador verifica que el resultado sea un orden).
    E-a: construir con `SJEvaluator.from_two_order(v)` y usar `propose_swap(a, b)` (transposicion de v_a y v_b).
    M = K^T K y B = C C se mantienen exactos (enteros en float64); `check_state()` los compara con el calculo completo.
    """

    def __init__(self, c: Bool, guard: float = GUARD, v: Int | None = None) -> None:
        self.c: Bool = np.array(c, dtype=bool)
        n = self.c.shape[0]
        if self.c.shape != (n, n):
            raise ValueError("c debe ser cuadrada")
        self.n = n
        self.guard = guard
        self.v: Int | None = None if v is None else np.array(v, dtype=np.int64)
        self.cf: F64 = self.c.astype(np.float64)
        self.k: F64 = self.cf - self.cf.T
        self.m: F64 = self.k.T @ self.k
        self.b: F64 = self.cf @ self.cf
        self._undo: _Undo | None = None
        self.n_eval = 0
        self.n_fallback = 0
        self.min_margin_global = math.inf
        self.min_margin_region = math.inf
        self.current = self._evaluate()

    @classmethod
    def from_two_order(cls, v: Int, guard: float = GUARD) -> SJEvaluator:
        v = np.asarray(v, dtype=np.int64)
        if not np.array_equal(np.sort(v), np.arange(v.size)):
            raise ValueError("v debe ser una permutacion de 0..N-1")
        idx = np.arange(v.size)
        c = (idx[:, None] < idx[None, :]) & (v[:, None] < v[None, :])
        return cls(c, guard, v)

    @property
    def s(self) -> float:
        return self.current.s

    @property
    def n_r(self) -> int:
        return self.current.n_r

    def _evaluate(self) -> SJResult:
        region = _select_region(self.c, np.rint(self.b).astype(np.int64) + 2, self.n // 8)
        if region.size == 0:
            r = _EMPTY
        else:
            r = _s_region(self.c, self.k, self.m, region, self.guard)
        self.n_eval += 1
        self.n_fallback += int(r.fallback)
        self.min_margin_global = min(self.min_margin_global, r.margin_global)
        self.min_margin_region = min(self.min_margin_region, r.margin_region)
        return r

    def _apply(self, s_idx: Int, new_rows: Bool, new_cols: Bool, swap: tuple[int, int] | None) -> SJResult:
        """Reemplaza las filas y columnas s_idx de C; actualiza K, M = K^T K y B = C C en O(|S| N^2), exacto."""
        if self._undo is not None:
            raise RuntimeError("hay una propuesta pendiente: commit() o rollback() primero")
        c, cf, k = self.c, self.cf, self.k
        old_rows, old_cols = c[s_idx, :].copy(), c[:, s_idx].copy()
        # pares (x, y) fuera de S: solo cambian los terminos con elemento intermedio m en S
        b_new = self.b - cf[:, s_idx] @ cf[s_idx, :]
        m_new = self.m - k[s_idx, :].T @ k[s_idx, :]
        c[s_idx, :] = new_rows
        c[:, s_idx] = new_cols
        cf[s_idx, :] = new_rows
        cf[:, s_idx] = new_cols
        k[s_idx, :] = cf[s_idx, :] - cf[:, s_idx].T
        k[:, s_idx] = cf[:, s_idx] - cf[s_idx, :].T
        b_new += cf[:, s_idx] @ cf[s_idx, :]
        m_new += k[s_idx, :].T @ k[s_idx, :]
        # filas y columnas S: recalculo directo
        b_new[s_idx, :] = cf[s_idx, :] @ cf
        b_new[:, s_idx] = cf @ cf[:, s_idx]
        ms = k.T @ k[:, s_idx]
        m_new[:, s_idx] = ms
        m_new[s_idx, :] = ms.T
        self._undo = _Undo(s_idx, old_rows, old_cols, self.m, self.b, swap, self.current)
        self.m, self.b = m_new, b_new
        if swap is not None and self.v is not None:
            a, b2 = swap
            self.v[a], self.v[b2] = self.v[b2], self.v[a]
        self.current = self._evaluate()
        return self.current

    def propose_toggle(self, i: int, j: int) -> SJResult:
        s_idx = np.array([i, j], dtype=np.int64)
        rows = self.c[s_idx, :].copy()
        rows[0, j] = not rows[0, j]
        cols = self.c[:, s_idx].copy()
        cols[i, 1] = rows[0, j]
        return self._apply(s_idx, rows, cols, None)

    def propose_swap(self, a: int, b: int) -> SJResult:
        if self.v is None:
            raise RuntimeError("propose_swap requiere un estado de 2-orden (from_two_order)")
        v = self.v.copy()
        v[a], v[b] = v[b], v[a]
        s_idx = np.array([a, b], dtype=np.int64)
        idx = np.arange(self.n)
        rows = (idx[None, :] > s_idx[:, None]) & (v[None, :] > v[s_idx][:, None])
        cols = (idx[:, None] < s_idx[None, :]) & (v[:, None] < v[s_idx][None, :])
        return self._apply(s_idx, rows, cols, (a, b))

    def commit(self) -> None:
        self._undo = None

    def rollback(self) -> None:
        u = self._undo
        if u is None:
            raise RuntimeError("no hay propuesta pendiente")
        self.c[u.s_idx, :] = u.rows
        self.c[:, u.s_idx] = u.cols
        self.cf[u.s_idx, :] = u.rows
        self.cf[:, u.s_idx] = u.cols
        self.k[u.s_idx, :] = self.cf[u.s_idx, :] - self.cf[:, u.s_idx].T
        self.k[:, u.s_idx] = self.cf[:, u.s_idx] - self.cf[u.s_idx, :].T
        self.m, self.b = u.m_old, u.b_old
        if u.swap is not None and self.v is not None:
            a, b = u.swap
            self.v[a], self.v[b] = self.v[b], self.v[a]
        self.current = u.result
        self._undo = None

    def check_state(self) -> bool:
        """M, B y K coinciden exactamente con el calculo completo desde C."""
        cf = self.c.astype(np.float64)
        k = cf - cf.T
        ok = np.array_equal(self.cf, cf) and np.array_equal(self.k, k)
        ok = ok and np.array_equal(self.m, k.T @ k) and np.array_equal(self.b, cf @ cf)
        if self.v is not None:
            idx = np.arange(self.n)
            ok = ok and np.array_equal(self.c, (idx[:, None] < idx[None, :]) & (self.v[:, None] < self.v[None, :]))
        return bool(ok)


# --- muestreadores con el mismo consumo de numeros aleatorios que mcmc_sj ------------------------------------


def _pairs(n: int, m: int, rng: np.random.Generator) -> tuple[Int, Int]:
    a = rng.integers(0, n, size=m)
    b = rng.integers(0, n - 1, size=m)
    b = b + (b >= a)
    return a, b


def _accept(lam: float, ds: float, u: float) -> bool:
    x = lam * ds
    return x >= 0 or u < math.exp(x)


def sweep_two_order_fast(ev: SJEvaluator, lam: float, rng: np.random.Generator) -> SweepCounts:
    """Igual que mcmc_sj.sweep_two_order (mismas propuestas y uniformes), con el evaluador incremental."""
    n = ev.n
    m = n * (n - 1) // 2
    a_arr, b_arr = _pairs(n, m, rng)
    u = rng.random(size=m)
    acc = 0
    for t in range(m):
        s_old = ev.s
        r = ev.propose_swap(int(a_arr[t]), int(b_arr[t]))
        if _accept(lam, r.s - s_old, float(u[t])):
            ev.commit()
            acc += 1
        else:
            ev.rollback()
    return SweepCounts(m, m, acc)


def sweep_general_fast(ev: SJEvaluator, lam: float, rng: np.random.Generator) -> SweepCounts:
    """Igual que mcmc_sj.sweep_general (mismas propuestas, validez y uniformes), con el evaluador incremental."""
    c = ev.c
    n = ev.n
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
        val += 1
        s_old = ev.s
        r = ev.propose_toggle(i, j)
        if _accept(lam, r.s - s_old, float(u[t])):
            ev.commit()
            acc += 1
        else:
            ev.rollback()
    return SweepCounts(m, val, acc)
