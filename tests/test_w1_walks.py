"""Pruebas rapidas de tools/w1_walks.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import w1_walks as W  # noqa: E402
from a0_rc2 import cycle, path  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402


def _csr(a):
    return a.indptr, a.indices.astype(np.int32)


def test_lazy_walk_stays_and_moves_on_path():
    ip, ix = _csr(path(5))
    start = np.full(20000, 2, dtype=np.int32)
    tr = W.lazy_walks(ip, ix, start, 1, rng_from_key((1, 2, 3)))
    nxt = tr[1]
    assert abs((nxt == 2).mean() - 0.5) < 0.02
    assert abs((nxt == 1).mean() - 0.25) < 0.02 and abs((nxt == 3).mean() - 0.25) < 0.02
    # solo vecinos o el propio vertice; en el extremo (grado 1) va al unico vecino
    tr = W.lazy_walks(ip, ix, np.zeros(2000, dtype=np.int32), 1, rng_from_key((4,)))
    assert set(np.unique(tr[1])) == {0, 1}
    # pasos largos: siempre transiciones validas
    tr = W.lazy_walks(ip, ix, np.full(50, 2, dtype=np.int32), 300, rng_from_key((5,)), chunk=7)
    assert np.abs(np.diff(tr.astype(int), axis=0)).max() <= 1


def test_l1_on_two_cycle_and_determinism():
    g = cycle(2)  # un solo arco 0-1
    ip, ix = _csr(g)
    K, T = 20000, 10
    tr = W.lazy_walks(ip, ix, np.zeros(K, dtype=np.int32), T, rng_from_key((7,)))
    assert (tr[1:] >= 0).all()
    visits = (tr[1:] == 0).sum(axis=0)
    assert abs(visits.mean() - T / 2) < 0.1  # P(en 0 en t>=1) = 1/2
    tr2 = W.lazy_walks(ip, ix, np.zeros(K, dtype=np.int32), T, rng_from_key((7,)))
    assert (tr == tr2).all()


def test_laws_on_handmade_trajectories_path():
    # camino 0-1-2-3-4; T=4 ; M=1; W1 = 2,3,4,3,2 ; W2 = 2,1,0,1,2 ; W3 = 2,2,3,3,2
    w1 = [2, 3, 4, 3, 2]
    w2 = [2, 1, 0, 1, 2]
    w3 = [2, 2, 3, 3, 2]
    traj = np.array([w1, w2, w3], dtype=np.int32).T[:, :, None]  # (5,3,1)
    C = W.law_counts(traj, (2, 4))
    # L1: visitas de W1 a 2 en pasos 1..t -> t=2: 0 ; t=4: 1
    assert C["L1"][0].tolist() == [0.0, 1.0]
    # L2: W1[t]==W2[t], t>=1 -> t=2: 0 ; t=4: 1 (t=4)
    assert C["L2"][0].tolist() == [0.0, 1.0]
    # L3: rangos hasta t=2: {2,3,4} y {2,1,0} -> {2}: 1 ; hasta 4: {2,3,4} y {0,1,2} -> {2}: 1
    assert C["L3"][0].tolist() == [1.0, 1.0]
    # L4 t=2: {2,3,4}∩{0,1,2}∩{2,3} = {2}: 1; t=4 igual
    assert C["L4"][0].tolist() == [1.0, 1.0]
    # interseccion mayor: dos paseos sobre 0,1,2
    traj = np.array([[0, 1, 2, 1, 0], [0, 1, 0, 1, 2], [0, 0, 0, 1, 1]], dtype=np.int32).T[:, :, None]
    C = W.law_counts(traj, (2, 4))
    # t=2: R1={0,1,2}, R2={0,1}, R3={0}; t=4: R1=R2={0,1,2}, R3={0,1}
    assert C["L3"][0].tolist() == [2.0, 3.0]
    assert C["L4"][0].tolist() == [1.0, 2.0]


def test_comb_generators_sizes_and_degrees():
    L = 7
    a = W.comb1(L)
    assert a.shape[0] == L * (L + 1)
    deg = np.diff(a.indptr)
    assert (deg[:L] == 3).all()  # ciclo (2) + diente (1)
    tooth = deg[L:].reshape(L, L)
    assert (tooth[:, :-1] == 2).all() and (tooth[:, -1] == 1).all()
    assert a.nnz // 2 == L + L * L  # aristas
    s = 4
    b = W.comb2(s)
    assert b.shape[0] == s * s * (s + 1)
    degb = np.diff(b.indptr)
    assert (degb[: s * s] == 5).all()  # Z2 (4) + diente (1)
    assert b.nnz // 2 == 2 * s * s + s * s * s


def test_diluted_and_cylinder_and_specs():
    g = W.diluted_z3(10, 0.6, rng_from_key((3,)))
    assert g.shape[0] == 1000 and 0.5 * 3000 < g.nnz // 2 < 0.7 * 3000
    ids = [s["id"] for s in W.specs(False)]
    assert len(ids) == len(set(ids)) == 33
    assert W.specs(False)[0]["geom"] == "Z1"


def test_classification_rules():
    assert W.classify(0.2, 0.01, 1000) == "CRECE"
    assert W.classify(0.10, 0.01, 1000) == "CRECE"
    assert W.classify(0.05, 0.01, 1000) == "ACOTADA"
    assert W.classify(0.07, 0.01, 1000) == "AMBIGUA"
    assert W.classify(0.2, 0.031, 1000) == "INSUFICIENTE"
    assert W.classify(0.2, 0.01, 159) == "INSUFICIENTE"  # 159//16 = 9 < 10
    assert W.classify(0.2, 0.01, 160) == "CRECE"
    assert W.classify(None, None, 1000) == "INSUFICIENTE"


def test_exponent_recovers_slope():
    Ts = (100, 400, 1600)
    C = np.tile(np.array([10.0, 20.0, 40.0]), (50, 1)) * (1 + 0.01 * np.random.default_rng(0).random((50, 3)))
    idx = np.random.default_rng(1).integers(0, 50, size=(200, 50))
    a, se, _ = W.exponent(C, Ts, idx)
    assert abs(a - 0.5) < 0.05 and se < 0.01


def _items(cls_by_d):
    return [(d, c) for d, c in cls_by_d]


def test_resolve_clean_and_critical_ambiguous():
    items = [(1, "CRECE"), (2, "CRECE"), (3, "ACOTADA"), (4, "ACOTADA")]
    r = W.resolve(items)
    assert r["status"] == "RESOLUBLE" and r["dc"] == [2.0, 3.0]  # cualquier corte entre 2 y 3
    items = [(1, "CRECE"), (2, "AMBIGUA"), (2, "CRECE"), (3, "ACOTADA"), (4, "ACOTADA")]
    r = W.resolve(items)  # d=2 AMBIGUA: d*=2 valido; 2.5 no (AMBIGUA con d<d*)
    assert r["dc"] == [2.0, 2.0]
    items = [(1, "CRECE"), (2, "AMBIGUA"), (3, "ACOTADA")]
    assert W.resolve(items)["dc"] == [2.0, 2.0]


def test_resolve_insufficient_and_no_cut():
    r = W.resolve([(1, "CRECE"), (2, "ACOTADA"), (3, "INSUFICIENTE"), (4, "INSUFICIENTE")])
    assert r["status"] == "NO_RESUELTA"
    # INSUFICIENTE exactamente en d* se permite
    r = W.resolve([(1, "CRECE"), (2, "CRECE"), (3, "INSUFICIENTE"), (4, "ACOTADA")])
    assert r["status"] == "RESOLUBLE" and r["dc"] == [3.0, 3.0]
    r = W.resolve([(1, "ACOTADA"), (2, "CRECE"), (3, "ACOTADA")])
    assert r["status"] == "NO_RESOLUBLE"


def test_stability_and_universality():
    r2 = {"status": "RESOLUBLE", "dc": [2.0, 2.0]}
    assert W.stability(r2, {"status": "RESOLUBLE", "dc": [2.0, 3.0]}) == "ESTABLE"
    assert W.stability(r2, {"status": "RESOLUBLE", "dc": [4.0, 4.0]}) == "INESTABLE-N"
    assert W.stability(r2, {"status": "NO_RESUELTA", "dc": None}) == "INESTABLE-N"
    dc = [2.0, 2.0]
    u = W.universality(dc, [("a", 1, "CRECE"), ("b", 3, "ACOTADA"), ("c", 2, "ACOTADA"), ("p", None, "AMBIGUA")])
    assert u["status"] == "UNIVERSAL"  # d_g=2 critico acepta cualquier clase
    u = W.universality(dc, [("a", 1, "ACOTADA")])
    assert u["status"] == "NO_UNIVERSAL"
    u = W.universality(dc, [("a", 3, "AMBIGUA")])
    assert u["status"] == "INDETERMINADA"
    u = W.universality(dc, [("a", 3, "AMBIGUA"), ("b", 1, "ACOTADA")])
    assert u["status"] == "NO_UNIVERSAL"


def _law(status, dc, stab="ESTABLE", univ="UNIVERSAL"):
    return {"res": {"status": status, "dc": dc}, "stab": stab, "univ": univ}


def test_outcome_order():
    nores = _law("NO_RESUELTA", None, "NA", "NA")
    assert W.outcome({"L1": nores, "L2": nores})["outcome"] == "W1-D"
    # resoluble pero INESTABLE-N no cuenta
    assert W.outcome({"L1": _law("RESOLUBLE", [2, 2], "INESTABLE-N", "NA")})["outcome"] == "W1-D"
    # disjuntos -> W1-C aunque haya una universal
    o = W.outcome({"L1": _law("RESOLUBLE", [2, 2]), "L4": _law("RESOLUBLE", [3, 3], univ="NO_UNIVERSAL"), "L3": nores})
    assert o["outcome"] == "W1-C" and o["disjoint_pairs"] == [("L1", "L4")]
    # solapan + una universal -> W1-B
    o = W.outcome({"L1": _law("RESOLUBLE", [2, 2]), "L2": _law("RESOLUBLE", [2, 3], univ="NO_UNIVERSAL")})
    assert o["outcome"] == "W1-B"
    # solapan, ninguna universal -> W1-B-
    o = W.outcome({"L1": _law("RESOLUBLE", [2, 2], univ="NO_UNIVERSAL"), "L2": _law("RESOLUBLE", [2, 2], univ="INDETERMINADA")})
    assert o["outcome"] == "W1-B-"


def test_geometry_class_seed_rule():
    assert W.geometry_class(["CRECE", "CRECE"]) == "CRECE"
    assert W.geometry_class(["CRECE", "ACOTADA"]) == "AMBIGUA"
    assert W.geometry_class(["CRECE", "INSUFICIENTE"]) == "INSUFICIENTE"


def test_run_cell_smoke_small():
    spec = {"id": "x", "geom": "Z2", "layer": 1, "d": 2, "d_g": None, "size_idx": 0, "seed": 0, "gen": "torus2", "param": 16}
    r = W.run_cell(spec)
    assert r["N_giant"] == 256 and r["T_max"] >= 1
    assert set(r["laws"]) == set(W.LAWS)
