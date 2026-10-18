"""Pruebas de omega/dynamics/sq.py (R1-1)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from p1d3_panel import erdos_renyi, ring, torus_lattice  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.dynamics.sq import (  # noqa: E402
    SWEEPS,
    Graph,
    Rnd,
    choose_target,
    edge_flags,
    exact_distance_set,
    run_dynamics,
    s_edge,
    t_edge,
)


def _adj(edges: list[tuple[int, int]], n: int) -> list[set[int]]:
    adj: list[set[int]] = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def test_s_edge_c4_true() -> None:
    adj = _adj([(0, 1), (1, 2), (2, 3), (3, 0)], 4)
    assert all(s_edge(adj, u, v) for u, v in [(0, 1), (1, 2), (2, 3), (3, 0)])


def test_s_edge_tree_false() -> None:
    adj = _adj([(0, 1), (0, 2), (1, 3), (1, 4), (2, 5)], 6)
    assert not any(s_edge(adj, u, v) for u, v in [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5)])


def test_s_edge_triangle_false_and_pendant_on_c4() -> None:
    assert not s_edge(_adj([(0, 1), (1, 2), (0, 2)], 3), 0, 1)
    adj = _adj([(0, 1), (1, 2), (2, 3), (3, 0), (0, 4)], 5)
    assert not s_edge(adj, 0, 4) and s_edge(adj, 0, 1)


def test_s_edge_z2_torus_all_true() -> None:
    g = Graph.from_csr(torus_lattice(6, 2))
    assert all(s_edge(g.adj, *divmod(c, g.n)) for c in g.el)


def test_t_edge_triangle() -> None:
    adj = _adj([(0, 1), (1, 2), (0, 2), (2, 3)], 4)
    assert t_edge(adj, 0, 1) and t_edge(adj, 1, 2) and not t_edge(adj, 2, 3)


@pytest.mark.parametrize("o", ["SYNC", "ASYNC", "RAND"])
def test_bs_z3_torus_fixed_point(o: str) -> None:
    g = Graph.from_csr(torus_lattice(5, 3))
    before = g.edge_set()
    d, a = SWEEPS[o](g, "B-s", "R2", Rnd(rng_from_key((1, 2, 3))))
    assert (d, a) == (0, 0) and g.edge_set() == before
    res = run_dynamics(g, "B-s", "R2", o, rng_from_key((1, 2, 4)), T=5)
    assert res["absorbed"] and res["sweeps"] == 0


@pytest.mark.parametrize("o", ["SYNC", "ASYNC"])
def test_bt_z3_torus_deletes_all(o: str) -> None:
    g = Graph.from_csr(torus_lattice(5, 3))
    e0 = g.n_edges
    d, a = SWEEPS[o](g, "B-t", "R2", Rnd(rng_from_key((1, 2, 3))))
    assert d == e0 and a > 0


def test_bt_rand_deletes_some() -> None:
    g = Graph.from_csr(torus_lattice(5, 3))
    d, _ = SWEEPS["RAND"](g, "B-t", "R2", Rnd(rng_from_key((1, 2, 3))))
    assert d > 0


def test_r2_picks_distance_exactly_two() -> None:
    n = 40
    adj = _adj([(i, i + 1) for i in range(n - 1)], n)  # camino 0..39
    rnd = Rnd(rng_from_key((7,)))
    assert exact_distance_set(adj, 0, 2) == {2} and exact_distance_set(adj, 0, 3) == {3}
    for v in (0, n - 1):
        for _ in range(20):
            t = choose_target(adj, v, "R2", rnd)
            assert t == (2 if v == 0 else n - 3)
    adj2 = _adj([(0, 1), (1, 2), (1, 3), (3, 4), (2, 5), (5, 6)], 8)  # 7 aislado
    dist2 = {2, 3}
    assert exact_distance_set(adj2, 0, 2) == dist2
    assert {choose_target(adj2, 0, "R2", rnd) for _ in range(60)} == dist2
    assert {choose_target(adj2, 0, "R3", rnd) for _ in range(60)} == {4, 5}
    t7 = {choose_target(adj2, 7, "R2", rnd) for _ in range(200)}  # grado 0: uniforme no local
    assert 7 not in t7 and len(t7) > 3


def test_fallback_when_no_candidate() -> None:
    adj = _adj([(0, 1)], 6)  # componente K2: sin vertices a distancia 2
    rnd = Rnd(rng_from_key((3,)))
    ts = {choose_target(adj, 0, "R2", rnd) for _ in range(100)}
    assert ts <= {2, 3, 4, 5} and len(ts) > 1


def test_a3_formula_equals_set_based_s() -> None:
    a = erdos_renyi(300, 900, rng_from_key((5, 5)))
    g = Graph.from_csr(a)
    u, v, s, t = edge_flags(a)
    m = sparse.csr_array(a)
    m.data[:] = 1
    a3 = (m @ m @ m).tocsr()
    deg = np.asarray(m.sum(axis=1)).ravel()
    for x, y, sv, tv in zip(u.tolist(), v.tolist(), s.tolist(), t.tolist()):
        assert sv == (s_edge(g.adj, x, y)) == (a3[x, y] - deg[x] - deg[y] + 1 > 0)
        assert tv == t_edge(g.adj, x, y)
    assert s.any() and (~s).any() and t.any()


def test_dense_graph_flags_match() -> None:
    a = erdos_renyi(120, 1500, rng_from_key((6, 6)))
    g = Graph.from_csr(a)
    u, v, s, t = edge_flags(a)
    assert all(sv == s_edge(g.adj, x, y) for x, y, sv in zip(u.tolist(), v.tolist(), s.tolist()))


@pytest.mark.parametrize("o", ["SYNC", "ASYNC", "RAND"])
def test_dynamics_deterministic_and_simple(o: str) -> None:
    def go() -> tuple[set[int], dict]:
        g = Graph.from_csr(ring(200, 2))
        r = run_dynamics(g, "B-s", "R3", o, rng_from_key((9, 9)), T=6)
        assert all(u not in g.adj[u] for u in range(g.n))
        assert all(v in g.adj[u] for u in range(g.n) for v in g.adj[u] if u != v) and len(g.el) == len(set(g.el))
        return g.edge_set(), r

    (e1, r1), (e2, r2) = go(), go()
    assert e1 == e2 and r1["trajectory"] == r2["trajectory"]
