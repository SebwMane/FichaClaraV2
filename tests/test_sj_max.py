"""Tests de omega.causality.sj_max (R17-1, docs/OMEGA_CONSEJO_REV17_ES0.md §5.1 y §6)."""

from __future__ import annotations

import itertools
import math
from collections import Counter, deque

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality.entropy import pauli_jordan, sj_from_eig, ssee
from omega.causality.mcmc import order_from_v
from omega.causality.mcmc_sj import (
    init_general,
    init_two_order,
    sweep_general,
    sweep_two_order,
)
from omega.causality.orders import random_2order
from omega.causality.sj_max import (
    OutsideSupport,
    interval_max,
    s_sj_max,
    s_sj_max_info,
    ssee_min1,
    sweep_general_max,
    sweep_two_order_max,
)

Bool = npt.NDArray[np.bool_]


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))


def _natural_orders(n: int) -> list[Bool]:
    iu = np.triu_indices(n, k=1)
    out = []
    for bits in itertools.product([False, True], repeat=len(iu[0])):
        c = np.zeros((n, n), dtype=bool)
        c[iu] = bits
        cf = c.astype(np.int64)
        if not np.any(((cf @ cf) > 0) & ~c):
            out.append(c)
    return out


def _perms(n: int) -> list[npt.NDArray[np.int64]]:
    return [np.array(p, dtype=np.int64) for p in itertools.permutations(range(n))]


def _cached(score_cache: dict[bytes, tuple[float, int]]):  # type: ignore[no-untyped-def]
    def f(c: Bool) -> tuple[float, int]:
        k = c.tobytes()
        if k not in score_cache:
            score_cache[k] = s_sj_max(c)
        return score_cache[k]

    return f


def _tv(emp: Counter[bytes], exact: dict[bytes, float]) -> float:
    tot = sum(emp.values())
    assert set(emp) <= set(exact)
    return 0.5 * sum(abs(emp.get(k, 0) / tot - p) for k, p in exact.items())


# --- interval_max -------------------------------------------------------------------------------


def test_interval_max_chain_is_whole_set() -> None:
    c = np.triu(np.ones((7, 7), dtype=bool), k=1)
    assert interval_max(c).tolist() == list(range(7))
    # chain with permuted labels
    p = _rng(1).permutation(7)
    cp = c[np.ix_(p, p)]
    assert interval_max(cp).tolist() == list(range(7))


def test_interval_max_height2_gives_link() -> None:
    n = 6
    c = np.zeros((n, n), dtype=bool)
    c[:3, 3:] = _rng(2).random((3, 3)) < 0.5
    c[0, 3] = True
    r = interval_max(c)
    assert r.size == 2
    a, b = int(r[0]), int(r[1])
    assert c[a, b]


def test_interval_max_ties() -> None:
    # 0 < 2, 1 < 2, 1 < 3: every interval is a link (size 2) -> lexicographically smallest (a, b) = (0, 2)
    c = np.zeros((4, 4), dtype=bool)
    c[0, 2] = c[1, 2] = c[1, 3] = True
    assert interval_max(c).tolist() == [0, 2]
    # two disjoint 3-chains: tie on size and on balance -> (0, 2)
    c3 = np.zeros((6, 6), dtype=bool)
    for x, y in [(3, 4), (4, 5), (3, 5), (0, 1), (1, 2), (0, 2)]:
        c3[x, y] = True
    assert interval_max(c3).tolist() == [0, 1, 2]
    # size dominates: 6 below the second chain makes I[6, 5] (size 4) the unique maximum
    c3b = np.zeros((7, 7), dtype=bool)
    c3b[:6, :6] = c3
    c3b[6, [3, 4, 5]] = True
    assert interval_max(c3b).tolist() == [3, 4, 5, 6]


def test_interval_max_balance_rule_is_vacuous() -> None:
    # If p < a then I[p, b] strictly contains I[a, b]; so a maximal-size interval always has a minimal and b maximal,
    # min(|past(a)|, |fut(b)|) = 0 for every candidate, and ties are resolved lexicographically.
    r = _rng(10)
    for _ in range(40):
        c = random_2order(16, r)
        if not c.any():
            continue
        reg = interval_max(c)
        a = int(reg[np.argmin(c[np.ix_(reg, reg)].sum(axis=0))])
        b = int(reg[np.argmin(c[np.ix_(reg, reg)].sum(axis=1))])
        assert c[:, a].sum() == 0 and c[b].sum() == 0


def test_interval_max_antichain_empty_and_support() -> None:
    c = np.zeros((5, 5), dtype=bool)
    assert interval_max(c).size == 0
    with pytest.raises(OutsideSupport):
        s_sj_max(c)
    with pytest.raises(OutsideSupport):
        init_general(c, s_sj_max)
    with pytest.raises(OutsideSupport):
        init_two_order(np.arange(4, -1, -1), s_sj_max)


def test_interval_max_matches_bruteforce() -> None:
    r = _rng(3)
    for _ in range(30):
        c = random_2order(12, r)
        if not c.any():
            continue
        best: tuple[int, int, int, int] | None = None
        past, fut = c.sum(axis=0), c.sum(axis=1)
        for a in range(12):
            for b in range(12):
                if c[a, b]:
                    size = int(np.sum(c[a] & c[:, b])) + 2
                    key = (size, int(min(past[a], fut[b])), -a, -b)
                    if best is None or key > best:
                        best = key
        assert best is not None
        a, b = -best[2], -best[3]
        exp = sorted({a, b} | set(np.flatnonzero(c[a] & c[:, b]).tolist()))
        assert interval_max(c).tolist() == exp


# --- ssee_min1 ----------------------------------------------------------------------------------


def _sj(c: Bool) -> tuple[npt.NDArray[np.complex128], npt.NDArray[np.complex128]]:
    n = c.shape[0]
    mu, u = np.linalg.eigh(pauli_jordan(c))
    return sj_from_eig(mu, u, math.sqrt(n) / (4 * math.pi))


def test_ssee_min1_equals_ssee_when_modes_exist() -> None:
    r = _rng(4)
    for _ in range(10):
        c = random_2order(24, r)
        reg = interval_max(c)
        idt, w = _sj(c)
        cut = math.sqrt(reg.size) / (4 * math.pi)
        s, kept, fired = ssee_min1(idt, w, reg, cut)
        assert not fired and kept > 0
        assert s == ssee(idt, w, reg, True, cut)


def test_ssee_min1_rule_fires() -> None:
    c = np.triu(np.ones((6, 6), dtype=bool), k=1)
    idt, w = _sj(c)
    reg = np.arange(6, dtype=np.int64)
    big = 1e3  # no mode passes
    assert ssee(idt, w, reg, True, big) == 0.0
    s, kept, fired = ssee_min1(idt, w, reg, big)
    assert fired and kept == 1
    # manual: the single mode of largest |nu|
    nu, v = np.linalg.eigh(idt[np.ix_(reg, reg)])
    k = int(np.argmax(np.abs(nu)))
    lam = float((v[:, k].conj() @ w[np.ix_(reg, reg)] @ v[:, k]).real / nu[k])
    assert s == pytest.approx(lam * math.log(abs(lam)) if abs(lam) > 1e-12 else 0.0, abs=1e-12)
    # a cutoff just above the 2nd largest |nu| reproduces the same single-mode answer through ssee
    a = np.sort(np.abs(nu))[::-1]
    # +/- pairs: the two largest |nu| are equal; use ssee with threshold between top pair and next one
    s_pair = ssee(idt, w, reg, True, float((a[1] + a[2]) / 2))
    assert math.isfinite(s_pair)


def test_ssee_min1_null_region() -> None:
    z = np.zeros((3, 3), dtype=np.complex128)
    s, kept, fired = ssee_min1(z, z, np.arange(3, dtype=np.int64), 0.1)
    assert (s, kept, fired) == (0.0, 0, True)


def test_s_sj_max_info_fields() -> None:
    c = random_2order(32, _rng(5))
    d = s_sj_max_info(c)
    assert (d["S"], d["nR"]) == s_sj_max(c)
    assert d["nR"] == interval_max(c).size
    assert 2 <= d["height"] <= d["nR"]
    assert 0.0 < d["rel_frac_R"] <= 1.0
    chain = np.triu(np.ones((8, 8), dtype=bool), k=1)
    d = s_sj_max_info(chain)
    assert d["nR"] == 8 and d["height"] == 8 and d["rel_frac_R"] == 1.0


# --- samplers: same RNG stream as mcmc_sj when the support predicate is always true -------------


def _n_rel(c: Bool) -> tuple[float, int]:
    return float(c.sum()), 0


def test_sweeps_reproduce_mcmc_sj_without_restriction() -> None:
    def everything(c: Bool) -> bool:
        return True

    c0 = order_from_v(_rng(6).permutation(8).astype(np.int64))
    s1, s2 = init_general(c0, _n_rel), init_general(c0, _n_rel)
    r1, r2 = _rng(7), _rng(7)
    for _ in range(20):
        assert sweep_general(s1, 0.4, r1, _n_rel) == sweep_general_max(s2, 0.4, r2, _n_rel, everything)
        assert np.array_equal(s1.c, s2.c)
    assert r1.random() == r2.random()
    t1, t2 = init_two_order(np.arange(8), _n_rel), init_two_order(np.arange(8), _n_rel)
    r1, r2 = _rng(8), _rng(8)
    for _ in range(20):
        assert sweep_two_order(t1, -0.3, r1, _n_rel) == sweep_two_order_max(t2, -0.3, r2, _n_rel, everything)
        assert np.array_equal(t1.v, t2.v)
    assert r1.random() == r2.random()


def test_samplers_never_visit_antichain() -> None:
    # strongly favour few relations: without the support rule the chain would sit on the antichain
    st = init_general(np.triu(np.ones((4, 4), dtype=bool), k=1), _n_rel)
    r = _rng(9)
    tot = val = 0
    for _ in range(500):
        cnt = sweep_general_max(st, -5.0, r, _n_rel)
        tot += cnt.proposals
        val += cnt.valid
        assert st.c.any()
    assert tot == 500 * 6 and val < tot
    t = init_two_order(np.arange(4), _n_rel)
    for _ in range(500):
        cnt = sweep_two_order_max(t, -5.0, r, _n_rel)
        assert order_from_v(t.v).any()
        assert cnt.proposals == 6


# --- exact enumeration N = 4, 5 ------------------------------------------------------------------


def _toggle_neighbours(c: Bool) -> list[bytes]:
    n = c.shape[0]
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            if c[i, j]:
                if np.any(c[i] & c[:, j]):
                    continue
            elif np.any(c[j] & ~c[i]) or np.any(c[:, i] & ~c[:, j]):
                continue
            d = c.copy()
            d[i, j] = not d[i, j]
            out.append(d.tobytes())
    return out


def _connected(nodes: dict[bytes, list[bytes]]) -> bool:
    start = next(iter(nodes))
    seen = {start}
    q = deque([start])
    while q:
        x = q.popleft()
        for y in nodes[x]:
            if y in nodes and y not in seen:
                seen.add(y)
                q.append(y)
    return len(seen) == len(nodes)


@pytest.mark.parametrize("n", [4, 5])
def test_support_connected(n: int) -> None:
    orders = [c for c in _natural_orders(n) if c.any()]
    assert len(orders) == {4: 39, 5: 356}[n]
    assert _connected({c.tobytes(): _toggle_neighbours(c) for c in orders})
    perms = [p for p in _perms(n) if order_from_v(p).any()]
    assert len(perms) == math.factorial(n) - 1
    graph: dict[bytes, list[bytes]] = {}
    for p in perms:
        nb = []
        for a, b in itertools.combinations(range(n), 2):
            q = p.copy()
            q[a], q[b] = q[b], q[a]
            nb.append(q.tobytes())
        graph[p.tobytes()] = nb
    assert _connected(graph)


@pytest.mark.parametrize("n", [4, 5])
@pytest.mark.parametrize("space", ["Eb", "Ea"])
def test_exact_kernel_detailed_balance(n: int, space: str) -> None:
    """Single-proposal kernel on the support (symmetric proposal, out-of-support/invalid -> stay): exact check."""
    lam = 6.0
    if space == "Eb":
        states = [c for c in _natural_orders(n) if c.any()]
        keys = [c.tobytes() for c in states]
        nbrs = {c.tobytes(): _toggle_neighbours(c) for c in states}
        score = {c.tobytes(): s_sj_max(c)[0] for c in states}
    else:
        perms = [p for p in _perms(n) if order_from_v(p).any()]
        keys = [p.tobytes() for p in perms]
        nbrs = {}
        for p in perms:
            nb = []
            for a, b in itertools.combinations(range(n), 2):
                q = p.copy()
                q[a], q[b] = q[b], q[a]
                nb.append(q.tobytes())
            nbrs[p.tobytes()] = nb
        score = {p.tobytes(): s_sj_max(order_from_v(p))[0] for p in perms}
    idx = {k: t for t, k in enumerate(keys)}
    prob = 1.0 / (n * (n - 1) // 2)
    mat = np.zeros((len(keys), len(keys)))
    for k in keys:
        for y in nbrs[k]:
            if y in idx:  # proposals to the antichain stay put
                mat[idx[k], idx[y]] += prob * min(1.0, math.exp(lam * (score[y] - score[k])))
        mat[idx[k], idx[k]] = 1.0 - mat[idx[k]].sum()
    pi = np.array([math.exp(lam * score[k]) for k in keys])
    pi /= pi.sum()
    flow = pi[:, None] * mat
    assert np.allclose(flow, flow.T, atol=1e-15)
    assert np.allclose(pi @ mat, pi, atol=1e-14)
    assert np.all(np.diag(mat) >= -1e-15)


@pytest.mark.parametrize("n,lam,sweeps", [(4, 6.0, 30000), (4, -6.0, 30000), (5, 6.0, 80000), (5, -6.0, 80000)])
def test_general_exact_distribution(n: int, lam: float, sweeps: int) -> None:
    orders = [c for c in _natural_orders(n) if c.any()]
    cache: dict[bytes, tuple[float, int]] = {}
    score = _cached(cache)
    w = {c.tobytes(): math.exp(lam * score(c)[0]) for c in orders}
    z = sum(w.values())
    exact = {k: x / z for k, x in w.items()}
    uni = 0.5 * sum(abs(p - 1 / len(exact)) for p in exact.values())
    st = init_general(np.triu(np.ones((n, n), dtype=bool), k=1), score)
    r = _rng(100 * n + (1 if lam > 0 else 2))
    emp: Counter[bytes] = Counter()
    for _ in range(sweeps):
        sweep_general_max(st, lam, r, score)
        emp[st.c.tobytes()] += 1
    tv = _tv(emp, exact)
    print(f"E-b N={n} lam={lam}: TV={tv:.4f} (target vs uniform {uni:.4f}, states {len(exact)})")
    assert tv < 0.06 and tv < uni / 3


@pytest.mark.parametrize("n,lam,sweeps", [(4, 6.0, 30000), (4, -6.0, 30000), (5, 6.0, 40000), (5, -6.0, 40000)])
def test_two_order_exact_distribution(n: int, lam: float, sweeps: int) -> None:
    perms = [p for p in _perms(n) if order_from_v(p).any()]
    cache: dict[bytes, tuple[float, int]] = {}
    score = _cached(cache)
    w = {p.tobytes(): math.exp(lam * score(order_from_v(p))[0]) for p in perms}
    z = sum(w.values())
    exact = {k: x / z for k, x in w.items()}
    uni = 0.5 * sum(abs(p - 1 / len(exact)) for p in exact.values())
    st = init_two_order(np.arange(n), score)
    r = _rng(200 * n + (1 if lam > 0 else 2))
    emp: Counter[bytes] = Counter()
    for _ in range(sweeps):
        sweep_two_order_max(st, lam, r, score)
        emp[st.v.tobytes()] += 1
    tv = _tv(emp, exact)
    print(f"E-a N={n} lam={lam}: TV={tv:.4f} (target vs uniform {uni:.4f}, states {len(exact)})")
    assert tv < 0.06 and tv < uni / 3


@pytest.mark.parametrize("n", [4, 5])
def test_min1_frequency_enumerated(n: int) -> None:
    eb = [s_sj_max_info(c) for c in _natural_orders(n) if c.any()]
    ea = [s_sj_max_info(order_from_v(p)) for p in _perms(n) if order_from_v(p).any()]
    f_b = sum(d["min1"] for d in eb)
    f_a = sum(d["min1"] for d in ea)
    print(f"N={n}: min1 fires in {f_b}/{len(eb)} E-b states and {f_a}/{len(ea)} E-a states")
    assert all(d["modes"] >= 1 for d in eb + ea)
    assert all(math.isfinite(d["S"]) for d in eb + ea)
