"""Pruebas rapidas de tools/w6_calib.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import w6_calib as W  # noqa: E402


def test_grid_box_degrees():
    for w, side in ((2, 6), (3, 5), (4, 4)):
        a = W.lattice_open(side, w)
        deg = np.asarray(a.sum(axis=1)).ravel()
        assert a.shape[0] == side**w
        assert deg.min() == w and deg.max() == 2 * w
        assert int((deg == 2 * w).sum()) == (side - 2) ** w  # interior
        assert int(a.nnz // 2) == w * side ** (w - 1) * (side - 1)


def test_size_choice_closest():
    assert W.size_param("torus", 4, 640_000) == 28
    assert W.size_param("torus", 2, 10_000) == 100
    assert W.size_param("fcc", 3, 640_000) == 54


def test_rgg_nested_subset():
    d = 2
    x = np.random.default_rng(0).random((4000, d))
    r = W.rgg_radius(d, n3=4000)
    small = W.rgg_subbox(x, 1 / 8, d)
    mid = W.rgg_subbox(x, 1 / 1, d)
    s = (1 / 8) ** (1 / d)
    assert (small < s).all() and small.shape[0] < mid.shape[0]
    # el conjunto de puntos pequeno es subconjunto del grande y el RGG pequeno es el inducido
    idx = np.flatnonzero(np.all(x < s, axis=1))
    assert np.array_equal(x[idx], small)
    big = W.rgg_open(x, r).toarray()
    sub = W.rgg_open(small, r).toarray()
    assert np.array_equal(big[np.ix_(idx, idx)], sub)
    # grado medio en el bulk ~ 8
    assert 6.5 < big.sum() / x.shape[0] < 8.5


def test_rgg_radius_formula():
    assert abs(np.pi * W.rgg_radius(2) ** 2 * W.N3_RGG - W.RGG_K) < 1e-9


def _m(fam, dd, rk, w5=True):
    return {"family": fam, "abs_ddelta": dd, "rel_dk": rk, "w5_pair1": w5, "w5_pair2": w5}


def test_threshold_formula():
    th = W.thresholds([_m("a", 0.001, 0.01), _m("b", 0.01, 0.02)])
    assert th["tau_d"] == 0.02 and th["tau_k"] == 0.05 and th["setter_d"] is None and th["setter_k"] is None
    th = W.thresholds([_m("a", 0.04, 0.03), _m("b", 0.02, 0.08)])
    assert abs(th["tau_d"] - 0.06) < 1e-12 and abs(th["tau_k"] - 0.12) < 1e-12
    assert th["setter_d"] == "a" and th["setter_k"] == "b"


def test_w6_rule_cases():
    assert W.w6_rule(_m("x", 0.01, 0.02), 0.03, 0.07)["valid"]
    assert not W.w6_rule(_m("x", 0.05, 0.02), 0.03, 0.07)["W6.2"]
    assert not W.w6_rule(_m("x", 0.01, 0.20), 0.03, 0.07)["W6.3"]
    assert not W.w6_rule(_m("x", 0.01, 0.02, w5=False), 0.03, 0.07)["W6.4"]
    assert not W.w6_rule(_m("x", None, 0.02), 0.03, 0.07)["W6.2"]
    assert W.w6_rule(_m("x", 0.03, 0.07), 0.03, 0.07)["valid"]  # frontera inclusiva


def test_calibrate_verdict_and_fragility():
    ctrl = [_m("c1", 0.01, 0.03), _m("c2", 0.01, 0.20)]  # c2 fija tau_k = 0.30
    j3 = [_m(f"j{i}", 0.01, 0.17) for i in range(3)]
    cal = W.calibrate(ctrl, j3)
    assert cal["thresholds"]["setter_k"] == "c2"
    assert not cal["P1"] and cal["verdict"] == "W6 INVALIDA"  # J3 pasa W6.3 con tau_k=0.30
    # sin c2, tau_k = 0.05 y J3 pasa a invalido: el veredicto cambia -> potencia fragil
    assert not cal["P2"] and cal["P2_label"] == "potencia fragil"
    j3b = [_m(f"j{i}", 0.01, 0.40) for i in range(3)]
    cal = W.calibrate(ctrl, j3b)
    assert cal["P1"] and cal["P2"] and cal["verdict"] == "W6 VALIDA" and cal["P3"]
    cal = W.calibrate([_m("c1", 0.1, 0.01)], j3b)
    assert not cal["P3"] and cal["P3_label"] == "W6.2 no informativa"
