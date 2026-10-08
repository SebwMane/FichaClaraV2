"""Tests de tools/cf1_growth.py (docs/OMEGA_CF1_PRERREGISTRO.md)."""

from __future__ import annotations

import numpy as np
import pytest

from omega.c0.references import rng_from_key
from tools import cf1_growth as cf
from tools.cf1_growth import Grower, hasse_digraph, boolean_digraph, relabel, seed_R20, seed_S, variant_config


def rngs(tag: int = 0):
    return {k: rng_from_key((cf.MASTER_CF1, 99, tag, 0, k)) for k in (0, 1, 3)}


def make(variant="V1", w=2, tag=0, c=None):
    up, down = seed_S(w)
    return Grower(up, down, rngs(tag), **variant_config(variant, c))


def assert_hasse(g: Grower) -> None:
    """Cada arista de cobertura es una cobertura real (sin atajos transitivos) y el orden es aciclico."""
    for b in range(g.n):
        for a in g.down[b]:
            for c in g.down[b]:
                if c != a:
                    assert a not in g.below(c), f"{a}<{b} no es cobertura (pasa por {c})"
    assert len(g._ranks()) == g.n


# ---------------------------------------------------------------- regla ER
def test_er_candidate_rule_hand_built():
    # 0 < 1,2 < 3 < 4
    up = [[1, 2], [3], [3], [4], []]
    down = [[], [0], [0], [1, 2], [3]]
    g = Grower(up, down, rngs())
    # up<down: 3 (1<2) y 4 (0<1); la raiz (down=0) y los atomos (up=down=1) no
    assert set(g.cand) == {3, 4}
    assert g.pend == {}


def test_er_candidates_after_closure_s2_and_audit():
    g = make("V1", 2)
    assert g.closure() and g.n == 4
    assert g.cand == [3]
    g.run([], None, 150)
    g.audit()  # cand == {p : up<down} exacto


def test_er_sin_candidatos_stops():
    g = Grower([[], [], []], [[], [], []], rngs())
    assert g.step() == "ER_SIN_CANDIDATOS"
    st, _ = g.run([], None, 10)
    assert st == "ER_SIN_CANDIDATOS"


def test_el_extends_anywhere_and_asynchronous_draws():
    g = make("V3", 2)
    assert g.cand == []  # EL no usa la lista de candidatos
    st, _ = g.run([], None, 120)
    assert st == "COMPLETO" and g.n >= 120
    g.audit()
    ga = make("V1A", 2)
    ga.run([], None, 200)
    ga.audit()
    g1 = make("V1", 2)
    g1.run([], None, 200)
    assert ga.n_extend > g1.n_extend  # ASYNC extiende aun con pendientes


# ---------------------------------------------------------------- CAUSAL-LOCAL
def hand_local_example():
    # componente A: raiz 0, atomos 1,2,3; 1 tiene cubiertas 4,5; 2 tiene 6,7; 4 tiene 16,17
    # componente B: raiz 8, atomos 9,10; componente C: cadena 11<12<13, 13 con cubiertas 14,15
    n = 18
    up = [[] for _ in range(n)]
    down = [[] for _ in range(n)]

    def cov(a, b):
        up[a].append(b)
        down[b].append(a)

    for a, b in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 6), (2, 7), (4, 16), (4, 17), (8, 9), (8, 10),
                 (11, 12), (12, 13), (13, 14), (13, 15)]:
        cov(a, b)
    return up, down


def test_causal_local_eligibility_small_example():
    up, down = hand_local_example()
    g = Grower(up, down, rngs())
    assert set(g.pend) == {0, 1, 2, 4, 8, 13}
    assert g.blockers[0] == set() and g.blockers[8] == set() and g.blockers[13] == set()
    assert g.blockers[1] == {0} and g.blockers[2] == {0}
    assert g.blockers[4] == {0, 1}  # 0 y 1 son antecesores estrictos pendientes de 4
    z_elig = {z for z, _ in g.eligible_instances()}
    assert z_elig == {0, 8, 13}
    # 13 tiene rango 2 > rango de 0: CAUSAL-RANK no lo elegiria primero, CAUSAL-LOCAL si (no esta sobre 0)
    assert g.rank[13] == 2 and g.rank[0] == 0
    g.audit()
    gr = Grower(up, down, rngs(), sched="RANK")
    assert {z for z, _ in gr.eligible_instances()} == {0, 1, 2, 4, 8, 13}
    # la rank elige (rango minimo, -aridad, z, S): z=0
    gr.complete()
    assert gr.down[-1] == [1, 2]


def test_local_completes_blockers_first():
    up, down = hand_local_example()
    g = Grower(up, down, rngs())
    for _ in range(40):
        if not g.pend:
            break
        z_elig = {z for z, _ in g.eligible_instances()}
        assert all(not (g.below(z) & set(g.pend)) for z in z_elig)  # ningun pendiente en el cono pasado
        g.complete()
        g.audit()


# ---------------------------------------------------------------- crecimiento
def grid_lists(a, b):
    idx = {(i, j): i * b + j for i in range(a) for j in range(b)}
    up = [[] for _ in idx]
    down = [[] for _ in idx]
    for (i, j), k in idx.items():
        for d in ((1, 0), (0, 1)):
            q = (i + d[0], j + d[1])
            if q in idx:
                up[k].append(idx[q])
                down[idx[q]].append(k)
    return up, down


def test_product_grid_is_fixed_point_degree_le4():
    up, down = grid_lists(6, 5)
    g = Grower(up, down, rngs())
    assert g.pend == {}  # una rejilla producto no tiene instancias C_flag
    assert max(len(u) + len(d) for u, d in zip(g.up, g.down)) <= 4
    g.audit()


def test_growth_from_s2_is_a_valid_poset_and_matches_full_recompute():
    # NOTA (hallazgo): la regla literal ER + C_flag NO conserva una rejilla producto desde S2: la extension anade UNA cubierta y
    # la compleccion crea duplicados, de modo que el grado de subida supera 4. Aqui solo se exige validez y mantenimiento exacto.
    g = make("V1", 2)
    st, _ = g.run([], None, 300)
    assert st == "COMPLETO" and g.n == 300
    assert_hasse(g)
    g.audit()


@pytest.mark.parametrize("variant,w,c", [("V1", 2, None), ("V1", 3, None), ("V1R", 3, None), ("V1A", 2, None),
                                         ("V3", 2, None), ("V4", 6, 3), ("V5", 3, None)])
def test_incremental_matches_full_recompute(variant, w, c):
    g = make(variant, w, c=c)
    for i in range(0, 120, 8):
        g.run([], None, g.n + 8)
        g.audit()


# ---------------------------------------------------------------- CAP
def test_cap_blocks_completions_and_candidates():
    g = make("V4", 6, c=2)
    g.run([], None, 300)
    assert all(len(g.up[i]) <= 2 for i in range(1, g.n))  # la raiz (up = 6) no recibe nada
    assert len(g.up[0]) == 6
    assert len(g.cap_blocked) > 0
    assert all(len(S) >= 2 for ss in g.pend.values() for S in ss)
    # candidatos: |up| < |down| y |up| < c
    assert all(len(g.up[p]) < len(g.down[p]) and len(g.up[p]) < 2 for p in g.cand)
    g.audit()
    for c in (2, 3, 4):
        h = make("V4", 6, tag=c, c=c)
        h.run([], None, 300)
        assert max(len(h.up[i]) for i in range(1, h.n)) <= c
        h.audit()


def test_pair_rule_only_arity_two():
    g = make("V5", 3)
    g.run([], None, 150)
    assert all(len(S) == 2 for ss in g.pend.values() for S in ss)
    assert all(len(d) <= 2 for d in g.down[4:])


# ---------------------------------------------------------------- determinismo
def test_determinism_and_seed_sensitivity():
    a = make("V1", 3, tag=5)
    b = make("V1", 3, tag=5)
    c = make("V1", 3, tag=6)
    for g in (a, b, c):
        g.run([], None, 400)
    assert a.up == b.up and a.down == b.down
    assert a.up != c.up
    sp = [s for s in cf.make_specs((0,)) if s["id"] == "V1_w2_s0"][0]
    r1, _ = cf.run_spec(sp, (50, 100, 200), (0.7, 0.05), do_measure=False)
    r2, _ = cf.run_spec(sp, (50, 100, 200), (0.7, 0.05), do_measure=False)
    assert r1 == r2
    assert [s["n"] for s in r1["snapshots"]] == [50, 100, 200]


# ---------------------------------------------------------------- semillas
def test_r20_seed_and_v7():
    up, down = seed_R20(rng_from_key((cf.MASTER_CF1, 8, 20, 0, 2)))
    g = Grower(up, down, rngs())
    mins = [i for i in range(20) if not g.down[i]]
    assert len(mins) >= 2 and all(g.rank[i] == 0 for i in mins)
    assert_hasse(g)
    g.run([], None, 150)
    g.audit()


def test_v2_coalescence_covers_one_maximal_of_each_component():
    sp = [s for s in cf.make_specs((0,)) if s["id"] == "V2_w22_s0"][0]
    r, _ = cf.run_spec(sp, (40, 100, 200), (0.7, 0.05), do_measure=False)
    assert r["coalescence_element"] == 40 and r["n_final"] >= 200
    a, b = r["coalescence_down"]
    assert a < 20 <= b
    sp = [s for s in cf.make_specs((0,)) if s["id"] == "V2_w23_s0"][0]
    r, _ = cf.run_spec(sp, (40, 100, 200), (0.7, 0.05), do_measure=False)
    assert r["coalescence_element"] == 40 and r["coalescence_down"][0] < 20 <= r["coalescence_down"][1]


# ---------------------------------------------------------------- V0
@pytest.mark.parametrize("sched", ["LOCAL", "RANK"])
def test_v0_closure_boolean_and_relabel(sched):
    for w in (3, 4, 5, 7):
        g = make("V1R" if sched == "RANK" else "V1", w)
        assert g.closure(5000) and g.n == 2 ** w
        import networkx as nx
        assert nx.is_isomorphic(hasse_digraph(g.up), boolean_digraph(w))
    for w in (3, 4, 5):
        perm = [int(x) for x in np.random.default_rng(w).permutation(w + 1)]
        pu, pd = relabel(*seed_S(w), perm)
        h = Grower(pu, pd, rngs(), sched=sched)
        assert h.closure(5000) and h.n == 2 ** w
        import networkx as nx
        assert nx.is_isomorphic(hasse_digraph(h.up), boolean_digraph(w))
