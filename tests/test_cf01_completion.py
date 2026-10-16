"""Tests de tools/cf01_completion.py (docs/OMEGA_CF0.md §8)."""

from __future__ import annotations

import itertools

from tools.cf01_completion import FLAG, Poset, boolean_check, run_one


def _boolean(w):
    subs = [frozenset(c) for r in range(w + 1) for c in itertools.combinations(range(w), r)]
    idx = {s: i for i, s in enumerate(subs)}
    up = [set() for _ in subs]
    down = [set() for _ in subs]
    for s in subs:
        for e in range(w):
            if e not in s:
                up[idx[s]].add(idx[s | {e}])
                down[idx[s | {e}]].add(idx[s])
    return up, down


def test_c2_w2_gives_b2():
    r, _ = run_one(2, 2, "ALTA", 0)
    assert r["terminated"] and r["n"] == 4 and r["maximal"] == 1 and r["boolean"]


def test_boolean_checker_b3_true():
    up, down = _boolean(3)
    assert boolean_check(up, down, 3)["isomorphic"]


def test_boolean_checker_non_boolean():
    # 8 elementos con rangos 1,3,3,1 pero cubiertas distintas (un atomo sin cubrir a nadie del nivel 2 correcto)
    up, down = _boolean(3)
    # romper una arista de cobertura y crear otra equivocada: cambia la estructura
    a, b = 1, 4
    assert b in up[a]
    up[a].discard(b)
    down[b].discard(a)
    assert not boolean_check(up, down, 3)["isomorphic"]
    # tamano incorrecto
    up2, down2 = _boolean(2)
    assert not boolean_check(up2, down2, 3)["isomorphic"]


def test_arity3_single_element_for_3_corner():
    # w=3 con los tres pares ya cerrados: una esquina de 3 cuadrados. La regla de aridad 3 crea UN techo.
    P = Poset(3, 3, "ALTA", 0)
    for pair in [(1, 2), (1, 3), (2, 3)]:
        P.new(frozenset(pair))
    P.inst, P.by_arity = {}, {}
    for z in range(P.n):
        P.set_z(z)
    assert len(P.by_arity[3]) == 1  # (los pares en z=1,2,3 son aridad 2; ALTA da prioridad al cubo)
    P.cap = P.n + 1
    P.run()
    assert P.n == 8 and P.down[7] == {4, 5, 6}  # un unico elemento, que cubre los tres sub-supremos
    P.cap = 5000
    assert P.run()
    assert P.n == 8 and P.maximal() == 1
    assert boolean_check(P.up, P.down, 3)["isomorphic"]


def test_deterministic():
    a, _ = run_one(3, 3, "RONDA", 2, 400)
    b, _ = run_one(3, 3, "RONDA", 2, 400)
    assert a == b


def test_incremental_matches_full_recompute():
    P = Poset(4, FLAG, "BAJA", 1, 300)
    P.run()
    inc = {zs: kk for m in P.by_arity.values() for zs, kk in m.items()}
    assert inc == P.full_instances()


def test_causal_w3_flag_gives_b3():
    for o in ("CAUSAL", "RONDA-CAUSAL"):
        r, _ = run_one(3, FLAG, o, 0)
        assert r["terminated"] and r["n"] == 8 and r["boolean"], o


def test_causal_deterministic():
    for o in ("CAUSAL", "RONDA-CAUSAL", "CAUSAL-INV"):
        a, _ = run_one(3, 3, o, 1, 300)
        b, _ = run_one(3, 3, o, 1, 300)
        assert a == b


def test_old_schedulers_reproduce_stored_results():
    import json
    from pathlib import Path
    path = Path(__file__).resolve().parent.parent / "results" / "cf01" / "runs.jsonl"
    stored = {(r["w"], r["k_code"], r["order"], r["perm"]): r
              for r in map(json.loads, path.read_text().splitlines())}
    for key in [(2, 2, "ALTA", 0), (3, 3, "ALTA", 0), (3, 99, "ALTA", 2), (4, 4, "ALTA", 0), (2, 99, "BAJA", 1),
                (3, 3, "RONDA", 0), (3, 2, "RONDA", 3), (3, 99, "BAJA", 0)]:
        w, k, o, p = key
        assert run_one(w, k, o, p)[0] == stored[key], key
