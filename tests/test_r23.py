"""Pruebas rápidas de tools/r23.py (R23, piloto exploratorio). Semillas propias (880023), nunca las entropías del
experimento. No dependen de las referencias de R20 con N0 = 8000: la prueba de punta a punta escribe un conjunto de
referencias planas pequeño en tmp_path y lo inyecta con --refs."""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SEED = 880023


def _m() -> ModuleType:
    tools = str(ROOT / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("r23")


def _rng(*key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(SEED, spawn_key=key)))


# --- 1. generadores: exactamente n elementos y orden parcial ------------------------------------------------


@pytest.mark.parametrize("n", [7, 200, 333])
def test_generadores_n_y_orden_parcial(n: int) -> None:
    m = _m()
    from omega.causality.orders import is_partial_order

    assert len(m.FAMILIAS) == 8 and len(set(m.NOMBRES)) == 8
    for i, fam in enumerate(m.NOMBRES):
        c = m.generar(fam, n, _rng(1, i, n), _rng(2, i, n))
        assert c.shape == (n, n) and c.dtype == np.bool_
        assert not np.any(np.diag(c))  # irreflexivo
        assert not np.any(c & c.T)  # antisimétrico
        assert is_partial_order(c)  # y transitivo
        if n >= 200:
            assert np.count_nonzero(c) > 0


def test_generadores_deterministas() -> None:
    m = _m()
    for i, fam in enumerate(m.NOMBRES):
        a = m.generar(fam, 120, _rng(3, i), _rng(4, i))
        b = m.generar(fam, 120, _rng(3, i), _rng(4, i))
        assert np.array_equal(a, b)


# --- 2. el control nuevo tiene N = |Y_q| ----------------------------------------------------------------------


def test_control_mismo_n(monkeypatch: pytest.MonkeyPatch) -> None:
    m = _m()
    monkeypatch.setattr(m, "ENT_PAD", SEED)
    for fam in ("plana", "D5_gemelos", "KR"):
        rows, _ = m.parent_rows(fam, 0, 240)
        assert [r["qi"] for r in rows] == [0, 1, 2, 3]
        assert rows[0]["N"] == 240  # q = 1: sin dilución
        for r in rows:
            assert r["obs"]["N"] == r["ctrl"]["N"] == r["N"]
        assert len({r["N"] for r in rows}) > 1  # la dilución cambia el tamaño


# --- 3. u y AUC valen 0.5 en promedio con datos nulos ---------------------------------------------------------


def test_u_rank_empates_y_nan() -> None:
    m = _m()
    nul = np.array([1.0, 2.0, 2.0, 3.0])
    assert m.u_rank(2.0, nul) == pytest.approx((1 + 0.5 * 2) / 4)
    assert m.u_rank(0.5, nul) == 0.0 and m.u_rank(9.0, nul) == 1.0
    assert math.isnan(m.u_rank(math.nan, nul))
    assert m.u_rank(2.0, np.array([1.0, np.nan, 3.0])) == 0.5  # nulos sin D quedan fuera


def test_u_auc_nulo_sintetico() -> None:
    m = _m()
    r = _rng(5)
    nul = r.exponential(size=100)
    # exacto: cada nulo contra el conjunto nulo (empate consigo mismo por la mitad) promedia 0.5
    assert np.mean([m.u_rank(float(d), nul) for d in nul]) == pytest.approx(0.5, abs=1e-12)
    # brazos nulos: AUC (media de u sobre 30 padres frente a 100 nulos) promedia 0.5 sobre repeticiones con conjuntos
    # nulos nuevos (condicional a un conjunto nulo fijo, el AUC esperado se aparta de 0.5 en ~0.03: §4)
    aucs = []
    for _ in range(600):
        nul_k = r.exponential(size=100)
        aucs.append(np.mean([m.u_rank(float(d), nul_k) for d in r.exponential(size=30)]))
    assert abs(float(np.mean(aucs)) - 0.5) < 0.01
    # bootstrap: la estimación es la media de los evaluables; los nan no cuentan
    x = np.array([0.2, 0.8, np.nan, 0.5])
    idx = m.r20.boot_indices(4, _rng(6), 200)
    est, lo, hi = m.boot_ci(x, idx, "media")
    assert est == pytest.approx(0.5)
    assert math.isnan(lo) or lo <= est <= hi


# --- 4. la regla de componentes por par reproduce la máscara completa ----------------------------------------


def _synth(r: np.random.Generator, n: int, p: int = 14) -> np.ndarray:
    z = r.normal(size=(n, p)) @ (np.eye(p) + 0.3 * r.normal(size=(p, p)))
    z[:, 3] = np.nan  # r_40 fuera del vector común
    return np.asarray(z + 5.0, dtype=np.float64)


def test_componentes_por_par_mascara_completa() -> None:
    m = _m()
    r = _rng(7)
    fit, cal, nul = _synth(r, 100), _synth(r, 100), _synth(r, 100)
    pb = m.PairBand(fit, cal, nul)
    assert pb.base == m.COMUNES and len(pb.base) == 13
    vy, vz = _synth(r, 1)[0], _synth(r, 1)[0]
    pr = pb.pair(vy, vz)
    # máscara completa, calculada aparte con r14b.Band y r20.conformal_radius
    comps = m.COMUNES
    b = m.r14b.Band(fit, fit[:0], comps, m.r14b.comp_floors())
    rad = m.r20.conformal_radius(b.dist(cal[:, list(comps)], sub=True))[0]
    dn = b.dist(nul[:, list(comps)], sub=True) / rad
    assert pr["evaluable"] and pr["dim"] == 13 and pr["S"] == comps
    assert pr["D_Y"] == pytest.approx(float(b.dist(vy)[0]) / rad, rel=1e-12)
    assert pr["D_Z"] == pytest.approx(float(b.dist(vz)[0]) / rad, rel=1e-12)
    assert pr["u_Y"] == pytest.approx(m.u_rank(float(b.dist(vy)[0]) / rad, dn))
    assert pr["n_nulo"] == 100
    # una componente que falta en un solo brazo sale para los dos (y para el radio y la nula)
    vz2 = vz.copy()
    vz2[1] = np.nan  # r_10
    pr2 = pb.pair(vy, vz2)
    assert pr2["dim"] == 12 and 1 not in pr2["S"]
    b2 = m.r14b.Band(fit, fit[:0], pr2["S"], m.r14b.comp_floors())
    rad2 = m.r20.conformal_radius(b2.dist(cal[:, list(pr2["S"])], sub=True))[0]
    assert pr2["D_Y"] == pytest.approx(float(b2.dist(vy)[0]) / rad2, rel=1e-12)
    # menos de 10 componentes: no evaluable
    vy3 = vy.copy()
    vy3[[0, 1, 2, 4]] = np.nan
    pr3 = pb.pair(vy3, vz)
    assert pr3["dim"] == 9 and not pr3["evaluable"] and math.isnan(pr3["D_Y"]) and math.isnan(pr3["u_Z"])


# --- clasificación (§6) -----------------------------------------------------------------------------------


def test_etiquetas_y_lectura() -> None:
    m = _m()
    mg = m.MARGEN
    ok = {"q=1": 0.0, "q=1/8": 0.0}
    et = m.etiquetas((0.96, 1.0), (0.5, 0.7), (-1.0, -mg - 0.01))
    assert et == {"visible_al_inicio": True, "cerca": True, "lejos": False, "mejora": True, "igual": False}
    assert m.lectura(et, ok)["lectura"] == m.COMPATIBLE
    et = m.etiquetas((0.96, 1.0), (0.96, 1.0), (-0.1, 0.1))
    lec = m.lectura(et, ok)
    assert lec["lectura"] == m.TAMANO and "lejos" in lec["motivo"] and lec["tambien_cumple"] == [m.SIN]
    et = m.etiquetas((0.96, 1.0), (0.96, 1.0), (-0.5, 0.0))
    assert m.lectura(et, ok)["lectura"] == m.SIN
    et = m.etiquetas((0.90, 1.0), (0.5, 0.7), (-1.0, -0.5))
    assert m.lectura(et, ok)["lectura"] == m.NO_APLICA
    et = m.etiquetas((0.96, 1.0), (0.6, 0.9), (-1.0, -0.5))
    assert m.lectura(et, ok)["lectura"] == m.INDET
    et = m.etiquetas((0.96, 1.0), (0.5, 0.7), (-1.0, -0.5))
    assert m.lectura(et, {"q=1": 0.0, "q=1/8": 0.21})["lectura"] == m.INDET
    assert m.lectura(et, {"q=1": 0.0, "q=1/8": 0.20})["lectura"] == m.COMPATIBLE
    et = m.etiquetas((math.nan, math.nan), (0.5, 0.7), (-1.0, -0.5))
    assert not et["visible_al_inicio"]
    # igual es inclusivo en ±ln 1.3; mejora es estricto
    assert m.etiquetas((1, 1), (1, 1), (-mg, mg))["igual"]
    assert not m.etiquetas((1, 1), (1, 1), (-1.0, -mg))["mejora"]


def test_semillas_distintas() -> None:
    m = _m()
    keys = [(m.K_PTS, 0, 8000, 0), (m.K_AUX, 0, 8000, 0), (m.K_DIL, 0, 8000, 0, 1), (m.K_SEL_Y, 0, 8000, 0, 1),
            (m.K_CTRL_PTS, 0, 8000, 0, 1), (m.K_CTRL_AUX, 0, 8000, 0, 1), (m.K_SEL_Z, 0, 8000, 0, 1), (m.K_BOOT, 0)]
    draws = {float(m.rng_for(m.ENT_PAD, *k).random()) for k in keys}
    assert len(draws) == len(keys) == 8
    assert m.EXPERIMENT_ID == 334 and m.ENT_NULO == 20261118 and m.ENT_PAD == 20261119


# --- punta a punta, pequeño ---------------------------------------------------------------------------------


def _fake_refs(m: ModuleType, refs: Path, n0: int, nref: int) -> None:
    """Referencias planas pequeñas con semillas de prueba, con la forma de results/omega_r20/refs."""
    for c, conj in enumerate(("ajuste", "calibracion")):
        for qi, q in enumerate(m.QS):
            for s in range(nref):
                r = _rng(9, c, qi, s)
                n = n0 if q >= 1.0 else int(r.binomial(n0, q))
                cm, _ = m.r14b.GEN["plana"](n, r, r)
                row = {"conjunto": conj, "qi": qi, "q": q, "s": s, **m.r14b.measure(cm, _rng(10, c, qi, s))}
                m.r19._write(refs / conj / f"q{qi}_s{s:03d}.json", {"row": row})


def test_end_to_end_small(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    m = _m()
    monkeypatch.chdir(ROOT)
    out = tmp_path / "out"
    refs = tmp_path / "refs"
    n0, npad, nnull, nref = 320, 3, 12, 19
    _fake_refs(m, refs, n0, nref)

    def run_args(**kw: object) -> argparse.Namespace:
        base = dict(part="all", family=None, parents=None, workers=1, budget_min=None, n0=n0, nparents=npad,
                    nnull=nnull, out=str(out), allow_dirty=True)
        base.update(kw)
        return argparse.Namespace(**base)

    m.cmd_run(run_args(part="parents", family="plana,KR", parents="0-1", workers=2))  # pool con fork
    assert m.parent_path(out, "KR", 1).exists() and not m.parent_path(out, "KR", 2).exists()
    m.cmd_run(run_args())
    m.cmd_run(run_args())  # reanuda: nada por correr
    with pytest.raises(SystemExit):
        m.cmd_run(run_args(nnull=nnull + 1))  # otra configuración en el mismo out
    man = json.loads((out / "manifest.json").read_text())
    assert man["official"] is False and man["n0"] == n0 and man["experiment_id"] == 334
    for fam in m.NOMBRES:
        for p in range(npad):
            d = json.loads(m.parent_path(out, fam, p).read_text())
            assert d["spec"]["official"] is False
            assert d["rows"][0]["N"] == n0
            for row in d["rows"]:
                assert row["obs"]["N"] == row["ctrl"]["N"] == row["N"]
    assert len(list((out / "nulos").glob("*.json"))) == 4 * nnull
    m.cmd_summarize(argparse.Namespace(out=str(out), refs=str(refs), nref=nref))
    s = json.loads((out / "summary.json").read_text())
    assert set(s["familias"]) == set(m.NOMBRES)
    validas = {m.COMPATIBLE, m.TAMANO, m.SIN, m.NO_APLICA, m.INDET}
    for fam, r in s["familias"].items():
        assert r["lectura"] in validas
        assert len(r["por_q"]) == 4 and len(r["D_Y"]) == npad
        for c in r["por_q"]:
            for k in ("AUC_Y", "AUC_Z"):
                if c[k] is not None:
                    assert 0.0 <= c[k] <= 1.0
    assert s["blanco"]["0"]["dim_base"] >= 10 and s["nref"] == nref and len(s["refs_sha256"]) == 64
    assert not s["procedencia"]["todas_oficiales"]
    assert (out / "tabla.csv").read_text().count("\n") == 1 + 8 * 4
    assert (out / "entorno.json").exists() and any((out / "logs").glob("run_*.log"))
    m.cmd_figure(argparse.Namespace(out=str(out), fig=str(tmp_path / "r23.png")))
    assert (tmp_path / "r23.png").stat().st_size > 10000
