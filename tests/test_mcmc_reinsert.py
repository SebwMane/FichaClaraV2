"""Pruebas exactas del movimiento R y del nucleo mixto (ES-C0b, OMEGA_ESC0B_PREREGISTRO.md §1).

Enumeracion exacta con N <= 5: nucleo de una propuesta construido con el codigo real (reinsert_step y local_step),
filas que suman 1, balance detallado, pi P = pi, irreducibilidad, simetria de la propuesta R, Delta S exacto
frente a bd_action_4d(_smeared), transitividad y etiquetado natural tras cada R aceptado; y prueba de variacion
total con una cadena larga de sweep_mixed.
"""

from __future__ import annotations

import itertools
import math
from functools import cache

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality.estimators import bd_action_4d, bd_action_4d_smeared
from omega.causality.mcmc_bd_eps import action_eps, init_state, sweep
from omega.causality.mcmc_reinsert import (
    _current_applicable,
    family_size,
    local_step,
    pair_at,
    reinsert_family,
    reinsert_step,
    sweep_mixed,
)
from omega.causality.orders import is_partial_order, random_2order, transitive_closure
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond

BoolMat = npt.NDArray[np.bool_]
Mat = npt.NDArray[np.float64]


def natural(c: BoolMat) -> BoolMat:
    perm = np.argsort(c.sum(axis=0), kind="stable")
    out: BoolMat = c[np.ix_(perm, perm)]
    assert not np.any(np.tril(out))
    return out


def full_between(c: BoolMat) -> npt.NDArray[np.int64]:
    cf = c.astype(np.int64)
    out: npt.NDArray[np.int64] = cf @ cf
    return out


@cache
def _orders(n: int) -> tuple[list[BoolMat], dict[bytes, int]]:
    pairs = list(itertools.combinations(range(n), 2))
    out = []
    for bits in itertools.product((False, True), repeat=len(pairs)):
        c = np.zeros((n, n), dtype=bool)
        for (a, b), v in zip(pairs, bits, strict=True):
            c[a, b] = v
        if np.array_equal(transitive_closure(c), c):
            out.append(c)
    return out, {c.tobytes(): k for k, c in enumerate(out)}


def s_ref(c: BoolMat, eps: float) -> float:
    return float(bd_action_4d(c)) if eps == 1.0 else bd_action_4d_smeared(c, eps)


def _check_state(st_c: BoolMat, st_b: npt.NDArray[np.int64], st_s: float, eps: float) -> None:
    assert not np.any(np.tril(st_c))
    assert is_partial_order(st_c)
    assert np.array_equal(st_b, full_between(st_c))
    if eps == 1.0:
        assert st_s == bd_action_4d(st_c)
    else:
        assert abs(st_s - bd_action_4d_smeared(st_c, eps)) <= 1e-9


# ---------------------------------------------------------------- familia


def _key(d: BoolMat, u: BoolMat) -> tuple[bytes, bytes]:
    return d.tobytes(), u.tobytes()


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_family_matches_count_and_order(n: int) -> None:
    orders, _ = _orders(n)
    for c in orders:
        for x in range(n):
            fam = reinsert_family(c, x)
            keys = [_key(d, u) for d, u in fam]
            assert len(set(keys)) == len(keys)  # deduplicada como conjuntos
            assert len(fam) == family_size(c, x)
            for k, (d, u) in enumerate(fam):
                y, z = pair_at(c, x, k)
                assert (y is None) == (not d.any()) and (z is None) == (not u.any())
                if y is not None:
                    assert d[y] and not np.any(d[y + 1:]) and np.array_equal(d[:y], c[:y, y])
                if z is not None:
                    assert u[z] and not np.any(u[:z]) and np.array_equal(u[z + 1:], c[z, z + 1:])
            # aplicable sii el par actual esta en F (como conjuntos)
            cur = _key(c[:, x].copy(), c[x].copy())
            assert _current_applicable(c, x) == (cur in set(keys))


def test_family_dedup_is_by_set_identity() -> None:
    """Sin dedup, D_y distintos ya son conjuntos distintos (y es su unico maximo): la familia explicita sin
    deduplicar tiene el mismo largo que la deduplicada. Se verifica construyendo todos los D_y y U_z a mano."""
    rng = np.random.default_rng(5)
    for _ in range(30):
        c = natural(random_2order(9, rng))
        for x in range(9):
            downs = {(c[:, y] | (np.arange(9) == y)).tobytes() for y in range(x)}
            ups = {(c[z] | (np.arange(9) == z)).tobytes() for z in range(x + 1, 9)}
            assert len(downs) == x and len(ups) == 8 - x
            assert len(reinsert_family(c, x)) == family_size(c, x)


# ---------------------------------------------------------------- nucleo exacto


def _r_moves(c: BoolMat, eps: float, x: int) -> list[tuple[BoolMat, float]] | None:
    """Destinos (C nuevo, Delta S) de cada par de F con el codigo real; None si R no es aplicable."""
    n = c.shape[0]
    fam = reinsert_family(c, x)
    nf = len(fam)
    out = []
    for k in range(nf):
        st = init_state(c, eps)
        ap, ac = reinsert_step(st, x, (k + 0.5) / nf, 0.0, 0.0)
        if not ap:
            assert not _current_applicable(c, x)
            assert np.array_equal(st.c, c)
            return None
        assert ac
        d, u = fam[k]
        exp = c.copy()
        exp[:, x] = d
        exp[x, :] = u
        assert np.array_equal(st.c, exp)
        _check_state(st.c, st.b, st.s, eps)
        ds = st.s - s_ref(c, eps)
        full = s_ref(exp, eps) - s_ref(c, eps)
        if eps == 1.0:
            assert ds == full
        else:
            assert abs(ds - full) <= 1e-12
        out.append((st.c.copy(), ds))
    assert n >= 1
    return out


def _check_threshold_r(c: BoolMat, eps: float, x: int, k: int, nf: int, beta: float, ds: float) -> None:
    a = min(1.0, math.exp(-beta * ds))
    v = (k + 0.5) / nf
    if a < 1.0:
        st = init_state(c, eps)
        assert reinsert_step(st, x, v, a * (1 - 1e-9), beta)[1]
        st = init_state(c, eps)
        acc = reinsert_step(st, x, v, a * (1 + 1e-9), beta)[1]
        assert not acc and np.array_equal(st.c, c) and st.s == init_state(c, eps).s
    else:
        st = init_state(c, eps)
        assert reinsert_step(st, x, v, 1 - 1e-12, beta)[1]


def _kernels(n: int, eps: float, beta: float) -> tuple[Mat, Mat, Mat]:
    """(P_local, P_R, Q_R): nucleos de una propuesta local y R, y la matriz de propuesta R (sin aceptacion)."""
    orders, index = _orders(n)
    m = len(orders)
    pl = np.zeros((m, m))
    pr = np.zeros((m, m))
    qr = np.zeros((m, m))
    w_pair = 2 / (n * (n - 1))
    for a, c in enumerate(orders):
        s0 = s_ref(c, eps)
        for i, j in itertools.combinations(range(n), 2):
            st = init_state(c, eps)
            va, ac = local_step(st, i, j, 0.0, 0.0)
            if not va:
                pl[a, a] += w_pair
                continue
            assert ac
            _check_state(st.c, st.b, st.s, eps)
            ds = st.s - init_state(c, eps).s
            assert abs(ds - (s_ref(st.c, eps) - s0)) <= 1e-12
            acc = min(1.0, math.exp(-beta * ds))
            if acc < 1.0:
                assert local_step(init_state(c, eps), j, i, acc * (1 - 1e-9), beta)[1]
                assert not local_step(init_state(c, eps), i, j, acc * (1 + 1e-9), beta)[1]
            b = index[st.c.tobytes()]
            pl[a, b] += w_pair * acc
            pl[a, a] += w_pair * (1 - acc)
        for x in range(n):
            moves = _r_moves(c, eps, x)
            if moves is None:
                pr[a, a] += 1 / n
                qr[a, a] += 1 / n
                continue
            nf = len(moves)
            for k, (cn, ds) in enumerate(moves):
                _check_threshold_r(c, eps, x, k, nf, beta, ds)
                acc = min(1.0, math.exp(-beta * ds))
                b = index[cn.tobytes()]
                qr[a, b] += 1 / (n * nf)
                pr[a, b] += acc / (n * nf)
                pr[a, a] += (1 - acc) / (n * nf)
    return pl, pr, qr


def _pi(n: int, eps: float, beta: float) -> Mat:
    orders, _ = _orders(n)
    s = np.array([s_ref(c, eps) for c in orders])
    w = np.exp(-beta * (s - s.min()))
    out: Mat = w / w.sum()
    return out


@cache
def _kernels_cached(n: int, eps: float, beta: float) -> tuple[Mat, Mat, Mat]:
    return _kernels(n, eps, beta)


CASES = [(4, 1.0, -0.5), (4, 0.25, 2.0), (5, 1.0, 0.7), (5, 1.0, -0.4), (5, 0.25, -3.0), (5, 0.25, 1.5)]


@pytest.mark.parametrize(("n", "eps", "beta"), CASES)
def test_exact_mixed_kernel(n: int, eps: float, beta: float) -> None:
    pl, pr, _ = _kernels_cached(n, eps, beta)
    pi = _pi(n, eps, beta)
    m = pi.size
    for p_r in (0.1, 0.5, 1.0):
        p = (1 - p_r) * pl + p_r * pr
        np.testing.assert_allclose(p.sum(axis=1), 1.0, atol=1e-12)
        flow = pi[:, None] * p
        np.testing.assert_allclose(flow, flow.T, atol=1e-12, rtol=0)
        np.testing.assert_allclose(pi @ p, pi, atol=1e-12, rtol=0)
        if p_r < 1.0:  # el nucleo mixto es irreducible
            reach = (p > 0) | np.eye(m, dtype=bool)
            for _ in range(m.bit_length() + 1):  # cierre por cuadrados booleanos
                reach = (reach.astype(np.int64) @ reach.astype(np.int64)) > 0
            assert np.all(reach)


@pytest.mark.parametrize("n", [3, 4, 5])
def test_r_proposal_symmetric(n: int) -> None:
    _, _, qr = _kernels_cached(n, 1.0, 0.7)
    np.testing.assert_allclose(qr.sum(axis=1), 1.0, atol=1e-14)
    np.testing.assert_allclose(qr, qr.T, atol=1e-15, rtol=0)
    off = qr - np.diag(np.diag(qr))
    assert np.count_nonzero(off) > 0  # R mueve algo


def test_r_moves_layers_n5() -> None:
    """R conecta en un paso ordenes que el movimiento local no conecta (salto colectivo)."""
    pl, pr, _ = _kernels_cached(5, 1.0, 0.7)
    jump = (pr > 0) & (pl == 0) & ~np.eye(pr.shape[0], dtype=bool)
    assert np.any(jump)


# ---------------------------------------------------------------- cadenas


def _starts(rng: np.random.Generator) -> list[BoolMat]:
    return [
        natural(causal_matrix(sprinkle_diamond(16, 4, rng))),
        natural(random_2order(12, rng)),
        np.zeros((10, 10), dtype=bool),
        np.triu(np.ones((9, 9), dtype=bool), k=1),
    ]


@pytest.mark.parametrize("beta", [-0.5, 0.0, 0.3])
@pytest.mark.parametrize("eps", [1.0, 0.25])
def test_mixed_chains_keep_action_and_between(eps: float, beta: float) -> None:
    rng = np.random.default_rng(int(1000 * (beta + 2)) + int(100 * eps))
    r_changed = 0
    for c in _starts(rng):
        st = init_state(c, eps)
        for _ in range(6):
            cnt = sweep_mixed(st, beta, rng, p_r=0.5)
            assert cnt["r_proposed"] + cnt["local_proposed"] == c.shape[0] * (c.shape[0] - 1) // 2
            assert cnt["r_changed"] <= cnt["r_accepted"] <= cnt["r_applicable"] <= cnt["r_proposed"]
            assert cnt["local_accepted"] <= cnt["local_valid"] <= cnt["local_proposed"]
            r_changed += cnt["r_changed"]
            _check_state(st.c, st.b, st.s, eps)
            assert st.s == pytest.approx(action_eps(st.c, st.b, st.g, st.const), abs=1e-9)
    assert r_changed > 0


def test_local_step_equals_bd_eps_sweep() -> None:
    """Con p_r = 0, cada paso local reproduce la logica de mcmc_bd_eps.sweep (mismas propuestas y u)."""
    rng0 = np.random.default_rng(3)
    for c in _starts(rng0):
        for eps, beta in ((1.0, 0.3), (0.25, -2.0)):
            a, e = init_state(c, eps), init_state(c, eps)
            ra = np.random.default_rng(11)
            for _ in range(3):
                n = c.shape[0]
                m = n * (n - 1) // 2
                state = ra.bit_generator.state
                acc, val = sweep(a, beta, ra)
                rb = np.random.default_rng()
                rb.bit_generator.state = state
                i_arr = rb.integers(0, n, size=m)
                j_arr = rb.integers(0, n - 1, size=m)
                j_arr = j_arr + (j_arr >= i_arr)
                u = rb.random(size=m)
                acc2 = val2 = 0
                for t in range(m):
                    v_, ac_ = local_step(e, int(i_arr[t]), int(j_arr[t]), float(u[t]), beta)
                    val2 += v_
                    acc2 += ac_
                assert (acc, val) == (acc2, val2)
                assert np.array_equal(a.c, e.c) and np.array_equal(a.b, e.b) and a.s == e.s


def test_long_chain_tv_n4_mixed() -> None:
    n, eps, beta = 4, 0.25, -3.0
    pi = _pi(n, eps, beta)
    assert 0.5 * np.abs(pi - 1 / pi.size).sum() > 0.15  # peso visiblemente no uniforme
    _, index = _orders(n)
    st = init_state(np.zeros((n, n), dtype=bool), eps)
    rng = np.random.default_rng(20261107)
    counts = np.zeros(pi.size)
    for _ in range(40_000):
        sweep_mixed(st, beta, rng, p_r=0.3)
        counts[index[st.c.tobytes()]] += 1
    tv = 0.5 * float(np.abs(counts / counts.sum() - pi).sum())
    assert tv < 0.02, tv
    _check_state(st.c, st.b, st.s, eps)
