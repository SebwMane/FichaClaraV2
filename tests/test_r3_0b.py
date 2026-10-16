"""Pruebas rapidas de tools/r3_0b.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import r3_0b as R  # noqa: E402


def _bits(m):
    return {i for i in range(m.bit_length()) if m >> i & 1}


def _enum(preds, edges=True):
    if preds and isinstance(preds[0], int):
        preds = R.masks_to_tuples(preds)
    return R.enumerate_downsets(preds, 10**9, edges)


def test_chains_count():
    for w, n in ((1, 7), (2, 5), (3, 4)):
        assert _enum(R.poset_j1(w, n), False)["count"] == (n + 1) ** w


def test_antichain_hypercube():
    k = 6
    e = _enum([()] * k)
    assert e["count"] == 2**k
    a = R.cover_graph(e["count"], e["eu"], e["ev"])
    assert (np.asarray(a.sum(axis=1)).ravel() == k).all()
    assert e["eu"].size == k * 2 ** (k - 1)


def test_j1_2_grid():
    n = 6
    e = _enum(R.poset_j1(2, n))
    a = R.cover_graph(e["count"], e["eu"], e["ev"])
    assert e["count"] == (n + 1) ** 2
    assert e["eu"].size == 2 * n * (n + 1)
    deg = np.asarray(a.sum(axis=1)).ravel()
    assert sorted(set(deg)) == [2, 3, 4]
    assert (deg == 2).sum() == 4 and (deg == 4).sum() == (n - 1) ** 2 and (deg == 3).sum() == 4 * (n - 1)
    assert abs(a - a.T).nnz == 0 and set(a.data) == {1.0}


def test_j2_cuts():
    w, nc = 3, 6
    preds = R.poset_j2(w, nc)
    # reconstruir cortes (c_i) desde los down-sets
    assert _enum(preds, False)["count"] < (nc + 1) ** w
    # enumerar a mano los down-sets y comprobar |c_i - c_j| <= 2
    seen, stack = {0}, [0]
    while stack:
        d = stack.pop()
        for x in range(len(preds)):
            if not d >> x & 1 and preds[x] & d == preds[x]:
                d2 = d | 1 << x
                if d2 not in seen:
                    seen.add(d2)
                    stack.append(d2)
    assert len(seen) == _enum(preds, False)["count"]
    for d in seen:
        c = [sum(d >> (k * w + i) & 1 for k in range(nc)) for i in range(w)]
        assert max(c) - min(c) <= 2


def test_closure():
    direct = [0, 1 << 0, 1 << 1, 1 << 2, 0b1]  # 0<1<2<3 ; 0<4
    cl = R.transitive_closure(direct)
    assert _bits(cl[3]) == {0, 1, 2} and _bits(cl[2]) == {0, 1} and _bits(cl[4]) == {0}
    # orden de etiquetas no topologico
    cl2 = R.transitive_closure([1 << 1, 1 << 2, 0])  # 2<1<0
    assert _bits(cl2[0]) == {1, 2}
    # la enumeracion con directos o clausura coincide
    assert _enum(direct, False)["count"] == _enum(cl, False)["count"]


def test_determinism_and_prefix():
    key = (R.MASTER_R3, 8, 0, 0)
    assert R.make_poset("J3", 0, 40, key) == R.make_poset("J3", 0, 40, key)
    assert R.make_poset("J4", 0.1, 40, key) == R.make_poset("J4", 0.1, 40, key)
    big, small = (R.poset_j4(n, 0.1, R.rng_from_key(key)) for n in (40, 30))
    mask = (1 << 30) - 1
    assert [b & mask for b in big[:30]] == small
    assert R.make_poset("J2", 3, 5, key) == R.covers(R.poset_j2(3, 5))
    # J3 transitivo y estricto
    p = R.poset_j3(50, R.rng_from_key(key))
    assert R.transitive_closure(p) == p
    assert all(not (m >> i & 1) for i, m in enumerate(p))


def test_cap_and_choose():
    r = R.enumerate_downsets([()] * 12, 100)
    assert not r["complete"]
    ch = R.choose_size("J1", 2, (R.MASTER_R3, 2, 0, 0), 10_000)
    assert ch["found"] and ch["n_param"] == 99  # (n+1)^2 >= 1e4 minimo n=99


def test_covers_equivalent_to_closure():
    key = (R.MASTER_R3, 9, 0, 0)
    cl = R.poset_j4(14, 0.3, R.rng_from_key(key))
    cv = R.covers(cl)
    assert _enum(R.masks_to_tuples(cl), False)["count"] == R.enumerate_downsets(cv, 10**9, False)["count"]
    rng = R.rng_from_key(key)
    direct = [0] + [R._row_to_int(rng.random(j) < 0.3) for j in range(1, 14)]
    assert R.transitive_closure(direct) == cl
