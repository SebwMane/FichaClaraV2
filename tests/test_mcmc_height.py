"""Pruebas de mcmc_height (R18-1, OMEGA_R18_WL_ALTURA_PREREGISTRO.md §5).

- Nucleo sesgado exacto con N <= 5 y un W arbitrario: filas que suman 1, balance detallado y pi P = pi para
  pi prop. a exp(-W(h)) sobre los ordenes con etiquetado natural, irreducibilidad.
- Altura, depth, codepth, B, S4 y numero de relaciones incrementales == recalculo completo despues de cada propuesta
  (aceptada, rechazada o invalida) en cadenas con N = 8..40.
- W constante: cadena identica bit a bit a mcmc_general.sweep con beta = 0.
- Reponderacion: con N = 5 el promedio reponderado por e^{W(h)} de una cadena sesgada larga reproduce el promedio
  exacto de la medida de conteo.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality import mcmc_general
from omega.causality import mcmc_height as mh
from omega.causality.orders import is_partial_order, kr_order

Bool = npt.NDArray[np.bool_]


def natural_orders(n: int) -> list[Bool]:
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    out = []
    for bits in itertools.product((False, True), repeat=len(pairs)):
        c = np.zeros((n, n), dtype=bool)
        for (i, j), x in zip(pairs, bits, strict=True):
            c[i, j] = x
        if is_partial_order(c):
            out.append(c)
    return out


def copy_state(st: mh.HeightState) -> mh.HeightState:
    return mh.HeightState(st.c.copy(), st.b.copy(), st.depth.copy(), st.codepth.copy(), st.h, st.s4, st.nrel)


def natural(c: Bool) -> Bool:
    perm = np.argsort(c.sum(axis=0), kind="stable")
    out: Bool = c[np.ix_(perm, perm)]
    return out


def assert_exact(st: mh.HeightState) -> None:
    c = st.c
    assert st.h == mcmc_general.height(c)
    assert np.array_equal(st.depth, mh.depth_full(c))
    assert np.array_equal(st.codepth, mh.codepth_full(c))
    cf = c.astype(np.int64)
    assert np.array_equal(st.b, cf @ cf)
    assert st.nrel == int(c.sum())
    assert math.isclose(st.s4, mcmc_general.action4(c, st.b), abs_tol=1e-9)


# ---------------------------------------------------------------- nucleo exacto


@pytest.mark.parametrize(("n", "hmax", "w"), [(4, 3, (0.3, -1.1)), (5, 4, (0.7, -0.4, 1.3)), (5, 3, (-0.8, 0.9))])
def test_exact_kernel(n: int, hmax: int, w: tuple[float, ...]) -> None:
    bias = mh.Bias(np.asarray(w, dtype=np.float64), hmax)
    orders = natural_orders(n)
    index = {c.tobytes(): k for k, c in enumerate(orders)}
    m = len(orders)
    npairs = n * (n - 1) // 2
    p = np.zeros((m, m))
    for a, c in enumerate(orders):
        st = mh.init_state(c)
        for i, j in itertools.combinations(range(n), 2):
            prop = mh.propose(st, i, j)
            assert np.array_equal(st.c, c)  # propose no modifica el estado
            if prop is None:
                p[a, a] += 1.0 / npairs
                continue
            new = copy_state(st)
            mh.apply(new, prop)
            assert_exact(new)  # incluye h' predicho == altura del orden propuesto
            assert prop.h_new == mcmc_general.height(new.c)
            bidx = index[new.c.tobytes()]
            q = mh.accept_prob(bias, st.h, prop.h_new)
            p[a, bidx] += q / npairs
            p[a, a] += (1.0 - q) / npairs
    assert np.allclose(p.sum(axis=1), 1.0)
    hs = np.array([mcmc_general.height(c) for c in orders])
    pi = np.exp(-bias.full(n)[hs])
    pi /= pi.sum()
    flow = pi[:, None] * p
    assert np.allclose(flow, flow.T, atol=1e-14)
    assert np.allclose(pi @ p, pi, atol=1e-14)
    # irreducible: el soporte de P (sin la diagonal) es conexo
    reach = {0}
    frontier = [0]
    while frontier:
        a = frontier.pop()
        for bb in np.flatnonzero(p[a] > 0):
            if int(bb) not in reach:
                reach.add(int(bb))
                frontier.append(int(bb))
    assert len(reach) == m
    assert set(hs.tolist()) == set(range(1, n + 1))  # se cubren alturas fuera de la ventana por ambos lados


# ---------------------------------------------------------------- incremental en cadenas


def _starts(n: int, rng: np.random.Generator) -> list[Bool]:
    return [np.zeros((n, n), dtype=bool), np.triu(np.ones((n, n), dtype=bool), k=1), natural(kr_order(n, rng))]


@pytest.mark.parametrize(("n", "sweeps"), [(8, 30), (13, 12), (21, 4), (40, 1)])
def test_incremental_height_along_chains(n: int, sweeps: int) -> None:
    rng = np.random.default_rng(1000 + n)
    hmax = 6
    seen_h: set[int] = set()
    count = {"steps": 0}

    def check(st: mh.HeightState) -> None:
        assert_exact(st)
        seen_h.add(st.h)
        count["steps"] += 1

    for c0 in _starts(n, rng):
        bias = mh.Bias(rng.normal(0.0, 1.5, size=hmax - mh.HMIN + 1), hmax)
        st = mh.init_state(c0)
        for _ in range(sweeps):
            mh.sweep(st, bias, rng, check=check)
        # tambien con Wang-Landau activo (W cambia durante el barrido)
        wl = mh.WLAcc(0.5, np.zeros(bias.w.size, dtype=np.int64))
        mh.sweep(st, bias, rng, wl, check=check)
    assert count["steps"] > 0
    assert len(seen_h) >= 3


def test_constant_w_matches_mcmc_general() -> None:
    n = 24
    rng0 = np.random.default_rng(7)
    c0 = natural(kr_order(n, rng0))
    st = mh.init_state(c0)
    gs = mcmc_general.init_state(c0)
    bias = mh.Bias(np.full(5, 3.0), 6)
    r1, r2 = np.random.default_rng(11), np.random.default_rng(11)
    for _ in range(15):
        a1 = mh.sweep(st, bias, r1)
        a2 = mcmc_general.sweep(gs, 0.0, r2)
        assert a1 == a2
        assert np.array_equal(st.c, gs.c)
        assert math.isclose(st.s4, gs.s, abs_tol=1e-9)
    assert r1.random() == r2.random()


def _trips_offline(hs: list[int]) -> tuple[int, int]:
    ext = [mh.sector(h) for h in hs if mh.sector(h) >= 0]
    runs = [x for k, x in enumerate(ext) if k == 0 or x != ext[k - 1]]
    crossings = max(len(runs) - 1, 0)
    first_lh = next((k for k in range(1, len(runs)) if runs[k] == 1), None)
    trips = 0 if first_lh is None else sum(1 for k in range(first_lh + 1, len(runs)) if runs[k] == 0)
    return crossings, trips


def test_trip_counting() -> None:
    n = 9
    rng = np.random.default_rng(3)
    bias = mh.Bias(np.array([0.0, -0.5, -1.5, -3.0]), 5)  # favorece alturas altas
    for c0 in _starts(n, rng)[:2]:
        st = mh.init_state(c0)
        hs = [st.h]

        def rec(s: mh.HeightState, hs: list[int] = hs) -> None:
            hs.append(s.h)

        for _ in range(300):
            mh.sweep(st, bias, rng, check=rec)
        cr, tr = _trips_offline(hs)
        assert (st.crossings, st.trips) == (cr, tr)
        assert st.trips >= 3


# ---------------------------------------------------------------- reponderacion


def test_reweighting_reproduces_counting_measure() -> None:
    n, hmax = 5, 4
    orders = natural_orders(n)
    npairs = n * (n - 1) // 2
    rel_exact = float(np.mean([c.sum() / npairs for c in orders]))
    h_exact = float(np.mean([mcmc_general.height(c) for c in orders]))
    bias = mh.Bias(np.array([0.0, -1.2, -2.0]), hmax)  # sesgo fuerte hacia h alto (h = 5 recortado a W(4))
    rng = np.random.default_rng(20261019)
    st = mh.init_state(np.zeros((n, n), dtype=bool))
    rel, hh = [], []
    for _ in range(12000):
        mh.sweep(st, bias, rng)
        rel.append(st.nrel / npairs)
        hh.append(st.h)
    hv = np.asarray(hh)
    w = np.exp(bias.full(n)[hv])
    nb = 40
    for x, exact in ((np.asarray(rel), rel_exact), (hv.astype(np.float64), h_exact)):
        r = float(np.sum(w * x) / np.sum(w))
        z = (w * (x - r)).reshape(nb, -1).sum(axis=1)
        se = math.sqrt(nb * float(z.var(ddof=1))) / float(w.sum())
        assert abs(r - exact) < 4 * se + 1e-3, (r, exact, se)
    # el sesgo cambio de verdad la distribucion de h muestreada
    assert float(hv.mean()) > h_exact + 0.3


def test_wang_landau_learns_ln_g() -> None:
    n, hmax = 5, 4
    orders = natural_orders(n)
    hs = np.array([mcmc_general.height(c) for c in orders])
    lng = np.log([np.count_nonzero(hs == h) for h in range(mh.HMIN, hmax + 1)])
    bias = mh.Bias(np.zeros(hmax - mh.HMIN + 1), hmax)
    st = mh.init_state(np.zeros((n, n), dtype=bool))
    # con N = 5 un barrido son 10 propuestas: se mira la planitud cada 2000 barridos (R18-1 usa 1 con N >= 32)
    log = mh.wang_landau(st, bias, np.random.default_rng(5), max_sweeps=200000, lnf_min=1e-4, check_every=2000)
    assert log["reached_lnf_min"]
    assert len(log["halvings"]) == 14
    # contabilidad: cada visita dentro de la ventana sumo el ln f vigente
    total = sum(e["lnf_antes"] * sum(e["hist"]) for e in log["halvings"]) + log["lnf_final"] * sum(log["hist_final"])
    assert math.isclose(float(bias.w.sum()), total, rel_tol=1e-12)
    for e in log["halvings"]:
        hist = np.asarray(e["hist"])
        assert np.min(hist) >= mh.FLATNESS * hist.mean()
    # W aprende ln g(h) salvo una constante, con el error de saturacion propio de Wang-Landau
    learned = bias.w - bias.w[0]
    assert np.allclose(learned, lng - lng[0], atol=0.3), (learned, lng - lng[0])
