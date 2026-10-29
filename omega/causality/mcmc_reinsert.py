"""Movimiento colectivo R (reinsercion de un elemento) y nucleo mixto R + local (ES-C0b, OMEGA_ESC0B_PREREGISTRO.md §1).

Movimiento R sobre un orden C con etiquetado natural (triangular superior estricta, transitiva):

1. x uniforme en {0..N-1}.
2. C' = C sin las relaciones de x (fila y columna x en cero). Sigue siendo un orden.
3. Familia F(x, C') de pares (D, U):
   D in {vacio} U {{y} U pasado_{C'}(y) : y < x (etiqueta)}, U in {vacio} U {{z} U futuro_{C'}(z) : z > x};
   el par es valido si D x U esta contenido en C'. Como D es un conjunto inferior con maximo y y U uno superior
   con minimo z, eso equivale a C'[y, z] (o D o U vacio).
4. Si (pasado(x), futuro(x)) no esta en F, el movimiento se rechaza (R no aplicable).
5. Si esta, se propone un par de F uniforme y se acepta con min(1, exp(-beta Delta S_eps)).

F depende solo de C' (igual en la ida y la vuelta): la propuesta es simetrica. Los D de y distintos son distintos
(y es el unico maximo de D_y) y lo mismo los U, asi que |F| = 1 + x + (N - 1 - x) + #{(y, z) : y < x < z, C[y, z]};
``reinsert_family`` arma la familia explicita y la deduplica como conjuntos (prueba de que el conteo es correcto).

Delta S_eps es exacto e incremental: con P0 = pasado(x), U0 = futuro(x), (D, U) el par nuevo,
B_nuevo[a, b] = B[a, b] - P0[a] U0[b] + D[a] U[b] para a, b != x, B_nuevo[:, x] = C' D, B_nuevo[x, :] = U C'.
Solo cambian los terminos de los pares (a, b) en (P0 x U0) U (D x U) y los pares con x.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

from omega.causality.mcmc_bd_eps import EpsState

__all__ = ["family_size", "local_step", "pair_at", "reinsert_family", "reinsert_step", "sweep_mixed"]

Bool = npt.NDArray[np.bool_]
Int = npt.NDArray[np.int64]


def _current_applicable(c: Bool, x: int) -> bool:
    """(pasado(x), futuro(x)) es un par de F(x, C'): ambos vacios o principales en C'."""
    past = c[:, x]
    p = np.flatnonzero(past)
    if p.size:
        y0 = int(p[-1])  # la etiqueta mayor del pasado es un maximal; principal sii pasado = {y0} U pasado(y0)
        want = c[:, y0].copy()
        want[y0] = True
        if not np.array_equal(past, want):
            return False
    fut = c[x]
    f = np.flatnonzero(fut)
    if f.size:
        z0 = int(f[0])
        want = c[z0].copy()
        want[z0] = True
        if not np.array_equal(fut, want):
            return False
    return True


def family_size(c: Bool, x: int) -> int:
    """|F(x, C')| (C' = c sin las relaciones de x; las filas y < x y columnas z > x no tocan a x)."""
    n = c.shape[0]
    return int(1 + x + (n - 1 - x) + int(np.count_nonzero(c[:x, x + 1:])))


def pair_at(c: Bool, x: int, k: int) -> tuple[int | None, int | None]:
    """Par k-esimo de F como (y, z) (None = vacio). Orden: (0, 0), (y, 0) y = 0..x-1, (0, z) z = x+1..N-1,
    y luego los (y, z) con C[y, z] en orden de filas."""
    n = c.shape[0]
    if k == 0:
        return None, None
    k -= 1
    if k < x:
        return k, None
    k -= x
    if k < n - 1 - x:
        return None, x + 1 + k
    k -= n - 1 - x
    flat = np.flatnonzero(c[:x, x + 1:])
    idx = int(flat[k])
    w = n - 1 - x
    return idx // w, x + 1 + idx % w


def _down(c: Bool, y: int | None) -> Bool:
    n = c.shape[0]
    if y is None:
        return np.zeros(n, dtype=bool)
    d: Bool = c[:, y].copy()
    d[y] = True
    return d


def _up(c: Bool, z: int | None) -> Bool:
    n = c.shape[0]
    if z is None:
        return np.zeros(n, dtype=bool)
    u: Bool = c[z].copy()
    u[z] = True
    return u


def reinsert_family(c: Bool, x: int) -> list[tuple[Bool, Bool]]:
    """Familia explicita F(x, C') como pares de vectores (D, U), deduplicada por identidad de conjuntos.

    Se arma desde C' literal (sin usar la estructura de etiquetas para la validez) y en el mismo orden que
    ``pair_at``; sirve de referencia para las pruebas."""
    n = c.shape[0]
    cp = c.copy()
    cp[x, :] = False
    cp[:, x] = False
    downs: list[Bool] = [np.zeros(n, dtype=bool)]
    for y in range(x):
        d = cp[:, y].copy()
        d[y] = True
        downs.append(d)
    ups: list[Bool] = [np.zeros(n, dtype=bool)]
    for z in range(x + 1, n):
        u = cp[z].copy()
        u[z] = True
        ups.append(u)
    out: list[tuple[Bool, Bool]] = []
    seen: set[tuple[bytes, bytes]] = set()

    def add(d: Bool, u: Bool) -> None:
        if not np.all(cp[np.ix_(d, u)]):
            return
        key = (np.packbits(d).tobytes(), np.packbits(u).tobytes())
        if key not in seen:
            seen.add(key)
            out.append((d, u))

    add(downs[0], ups[0])
    for d in downs[1:]:
        add(d, ups[0])
    for u in ups[1:]:
        add(downs[0], u)
    for d in downs[1:]:
        for u in ups[1:]:
            add(d, u)
    return out


def _accept(beta: float, ds: float, u: float) -> bool:
    # valido tambien con beta < 0 (auditoria A-22): no hay atajo 'ds <= 0 acepta'
    return -beta * ds >= 0 or u < math.exp(-beta * ds)


def reinsert_step(st: EpsState, x: int, v: float, u: float, beta: float) -> tuple[bool, bool]:
    """Una propuesta R sobre x; v in [0, 1) elige el par floor(v |F|), u la decision. Devuelve (aplicable, aceptada)."""
    c, b, g = st.c, st.b, st.g
    if not _current_applicable(c, x):
        return False, False
    nf = family_size(c, x)
    k = min(int(v * nf), nf - 1)
    y, z = pair_at(c, x, k)
    p0 = c[:, x].copy()
    u0 = c[x].copy()
    d = _down(c, y)
    up = _up(c, z)
    if np.array_equal(d, p0) and np.array_equal(up, u0):
        return True, True  # se propuso el par actual: no cambia nada (Delta S = 0, aceptado)
    # C' (sin x): columnas/filas que no son x quedan iguales
    cp = c.copy()
    cp[x, :] = False
    cp[:, x] = False
    # pares (a, b) sin x cuyo n cambia
    rows = np.flatnonzero(p0 | d)
    cols = np.flatnonzero(u0 | up)
    b_sub_old = b[np.ix_(rows, cols)]
    b_sub_new = (b_sub_old - np.outer(p0[rows], u0[cols]).astype(np.int64)
                 + np.outer(d[rows], up[cols]).astype(np.int64))
    rel_sub = cp[np.ix_(rows, cols)]  # todos relacionados en C' (P0 x U0 y D x U lo estan)
    ds = float(np.sum(g[b_sub_new[rel_sub]]) - np.sum(g[b_sub_old[rel_sub]]))
    # pares con x
    col_new = np.count_nonzero(cp[:, d], axis=1).astype(np.int64)  # B_nuevo[a, x] = #{m in D : a < m}
    row_new = np.count_nonzero(cp[up, :], axis=0).astype(np.int64)  # B_nuevo[x, b] = #{m in U : m < b}
    ds += float(np.sum(g[col_new[d]]) + np.sum(g[row_new[up]]) - np.sum(g[b[p0, x]]) - np.sum(g[b[x, u0]]))
    if not _accept(beta, ds, u):
        return True, False
    b[np.ix_(rows, cols)] = b_sub_new
    b[:, x] = col_new
    b[x, :] = row_new
    c[:, x] = d
    c[x, :] = up
    st.s += ds
    return True, True


def local_step(st: EpsState, i: int, j: int, u: float, beta: float) -> tuple[bool, bool]:
    """Una propuesta del movimiento local de mcmc_bd_eps.sweep sobre el par {i, j} (i != j); misma logica exacta.

    Devuelve (valida, aceptada)."""
    c, b, g = st.c, st.b, st.g
    if i > j:
        i, j = j, i
    succ = np.flatnonzero(c[j])
    pred = np.flatnonzero(c[:, i])
    if c[i, j]:
        if b[i, j] != 0:
            return False, False  # relacion que no es link: no se puede quitar
        sign = -1
        ds = -float(g[0])
    else:
        if np.any(~c[i, succ]) or np.any(~c[pred, j]):
            return False, False  # agregarla romperia la transitividad
        sign = 1
        ds = float(g[0])
    k1 = b[i, succ]
    k2 = b[pred, j]
    ds += float(np.sum(g[k1 + sign] - g[k1]) + np.sum(g[k2 + sign] - g[k2]))
    if not _accept(beta, ds, u):
        return True, False
    c[i, j] = sign > 0
    b[i, succ] += sign
    b[pred, j] += sign
    st.s += ds
    return True, True


def sweep_mixed(st: EpsState, beta: float, rng: np.random.Generator, p_r: float = 0.1) -> dict[str, int]:
    """N(N-1)/2 propuestas; cada una R con prob. p_r, si no la local de mcmc_bd_eps.

    Devuelve conteos: local_proposed, local_valid, local_accepted, r_proposed, r_applicable, r_accepted
    (r_accepted incluye las propuestas del par actual, que no cambian nada) y r_changed (aceptadas que cambian C)."""
    n = st.c.shape[0]
    m = n * (n - 1) // 2
    is_r = rng.random(size=m) < p_r
    i_arr = rng.integers(0, n, size=m)
    j_arr = rng.integers(0, n - 1, size=m)
    j_arr = j_arr + (j_arr >= i_arr)
    v_arr = rng.random(size=m)
    u_arr = rng.random(size=m)
    cnt = dict.fromkeys(("local_proposed", "local_valid", "local_accepted", "r_proposed", "r_applicable",
                         "r_accepted", "r_changed"), 0)
    for t in range(m):
        u = float(u_arr[t])
        if is_r[t]:
            cnt["r_proposed"] += 1
            x = int(i_arr[t])  # x uniforme: se reutiliza el i de la propuesta
            before = st.c[:, x].copy(), st.c[x].copy()
            ap, ac = reinsert_step(st, x, float(v_arr[t]), u, beta)
            cnt["r_applicable"] += ap
            cnt["r_accepted"] += ac
            if ac and not (np.array_equal(before[0], st.c[:, x]) and np.array_equal(before[1], st.c[x])):
                cnt["r_changed"] += 1
        else:
            cnt["local_proposed"] += 1
            va, ac = local_step(st, int(i_arr[t]), int(j_arr[t]), u, beta)
            cnt["local_valid"] += va
            cnt["local_accepted"] += ac
    return cnt
