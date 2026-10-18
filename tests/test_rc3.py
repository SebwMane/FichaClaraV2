"""Pruebas rapidas de tools/rc3.py."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

import rc3 as R  # noqa: E402

from omega.c0.references import rng_from_key  # noqa: E402
from omega.diagnostics.sampled_growth import sampled_ball_profile  # noqa: E402


def _deg(a):
    return np.asarray(a.sum(axis=1)).ravel()


def _R(adj, k=(1, 2, 3)):
    return R.half_mass_radius(sampled_ball_profile(adj, rng_from_key(k)))


def test_R_path_and_cycle():
    from a0_rc2 import cycle, path
    n = 201
    # ciclo: |B_r| = 2r+1 -> cruza n/2=100.5 en r=50 (m=101, m(49)=99): R = 49 + 1.5/2 = 49.75
    assert abs(_R(cycle(n)) - 49.75) < 1e-9
    # camino: R de media masa de orden n/4 (m(r) medio ~ 1 + 2r - boundary): solo cota
    assert 30 < _R(path(n)) < 100


def test_R_clique_half():
    assert abs(_R(R.clique_graph(100)) - 0.5) < 0.02


def test_R_not_crossed_returns_rcap():
    prof = {"n_giant": 1000, "mean": [3.0, 5.0, 6.0], "r_cap": 3}
    assert R.half_mass_radius(prof) == 3.0
    prof2 = {"n_giant": 10, "mean": [3.0, 8.0], "r_cap": 2}
    assert abs(R.half_mass_radius(prof2) - (1 + (5 - 3) / 5)) < 1e-12


def test_delta_cycles_about_one():
    from a0_rc2 import cycle
    # R(C_n) ~ n/4 (R_MAX=200 del observador congelado satura R: se usan n pequenos con R < 200)
    r1, r8 = _R(cycle(96)), _R(cycle(768))
    d = R.scaling_delta(r1, r8, 96, 768)
    assert abs(d - 1.0) < 0.02
    # saturacion documentada: C_8000 y C_1000 ambos con R = r_cap = 200 -> delta = 0
    assert R.scaling_delta(_R(cycle(1000)), _R(cycle(8000)), 1000, 8000) == 0.0


def test_fcc_bcc_degrees():
    f = R.fcc_lattice(4)
    assert f.shape[0] == 4 * 4**3 and (_deg(f) == 12).all()
    b = R.bcc_lattice(4)
    assert b.shape[0] == 2 * 4**3 and (_deg(b) == 8).all()
    assert R.fcc_lattice(14).shape[0] == 10976 and R.bcc_lattice(17).shape[0] == 9826


def test_other_generators():
    a = R.triangular_open(10)
    assert a.shape[0] == 100 and _deg(a).max() == 6 and _deg(a).min() == 2
    kt = R.k_tree(50, 4, rng_from_key((1, 1)))
    assert kt.nnz // 2 == 10 + 4 * (50 - 5)
    k3 = R.k_tree(40, 3, rng_from_key((1, 2)))
    assert k3.nnz // 2 == 3 * 40 - 6
    t = R.tree_of_cubes(7, 5)
    assert t.shape[0] == 875 and t.nnz // 2 == 7 * 300 + 6 * 25
    d4 = R.torus_with_diagonals(6, 4, 0.2, rng_from_key((1, 3)))
    assert d4.nnz // 2 == 6**4 * 4 + int(round(0.2 * 6**4 * 4))
    tc = R.cartesian(R.binary_tree(15), R.cycle(10))
    assert tc.shape[0] == 150


def test_exclusion_logic():
    f = R.pair_exclusions
    assert f(0.05, ("CONEXO",) * 2, ("COHERENTE",) * 2) == ["X1"]
    assert f(0.125, ("CONEXO",) * 2, ("COHERENTE",) * 2) == []  # estricto
    assert f(0.3, ("RAMIFICADO",) * 2, ("COHERENTE",) * 2) == ["X2"]
    assert f(0.3, ("RAMIFICADO", "CONEXO"), ("COHERENTE",) * 2) == []
    assert f(0.3, ("DOS_EXTREMOS",) * 2, ("COHERENTE",) * 2) == ["X3"]
    assert f(0.9, ("CONEXO",) * 2, ("COHERENTE",) * 2) == ["X3"]
    assert f(0.75, ("CONEXO",) * 2, ("COHERENTE",) * 2) == []
    assert f(0.3, ("CONEXO",) * 2, ("INCOHERENTE",) * 2) == ["X4"]
    assert f(0.3, ("CONEXO",) * 2, ("INCOHERENTE", "SIN_VENTANA")) == []
    assert f(None, ("NO_EVALUABLE",) * 2, ("NO_EVALUABLE",) * 2) == []
    assert f(0.01, ("RAMIFICADO",) * 2, ("INCOHERENTE",) * 2) == ["X1", "X2", "X4"]


def _row(fam, grp, d, new, exc):
    return {"id": fam, "family": fam, "group": grp, "d": d, "new": new, "excluded": exc}


def test_metrics_and_verdict():
    rows = [_row(f"v{i}", "V", 2 + i % 3, True, False) for i in range(20)]
    rows += [_row(f"x{i}", "X", None, True, True) for i in range(10)]
    rows += [_row("h", "E", None, False, True)]
    m = R.metrics([r for r in rows if r["new"]])
    assert m["FE"] == 0 and m["EP"] == 1 and m["verdict"] == "RC3-VALIDA"
    rows[0]["excluded"] = True  # 1/20 = 0.05 y 1 por d: sigue VALIDA
    assert R.metrics(rows)["verdict"] == "RC3-VALIDA"
    rows[3]["excluded"] = True  # d=2 dos exclusiones: FE=0.10 -> PARCIAL
    assert R.metrics(rows)["verdict"] == "RC3-PARCIAL"
    for i in range(20, 24):
        rows[i]["excluded"] = False  # EP = 0.6
    assert R.metrics(rows)["verdict"] == "RC3-INVALIDA"
    # familia con 2 instancias, 2 sin excluir -> no VALIDA aunque FE y EP globales pasen
    rows2 = [_row(f"v{i}", "V", 2, True, False) for i in range(20)]
    rows2 += [{**_row("f", "X", None, True, i < 0), "id": f"f{i}"} for i in range(2)]
    rows2 += [_row(f"y{i}", "X", None, True, True) for i in range(40)]
    m2 = R.metrics(rows2)
    assert m2["EP"] > 0.9 and not m2["valid_conditions"]["families_all_but_one"] and m2["verdict"] == "RC3-PARCIAL"
    assert R.metrics([_row("h", "E", None, True, True)])["verdict"] == "NO_EVALUABLE"


def test_specs_pairs():
    sp = R.specs(smoke=True)
    assert len(sp) == sum(2 if f[4] else 1 for f in R.FAMS)
    assert all(math.isclose(1, 1) and len(s["args"]) == 2 for s in sp)
