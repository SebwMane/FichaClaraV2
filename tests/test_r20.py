"""Pruebas rápidas de tools/r20.py (R20). Semillas propias (880002), nunca las entropías del experimento. Usan un
módulo sellado de prueba escrito en tmp_path (con IDENTIFICABLES), nunca el sellado real."""

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

DUMMY = '''
import numpy as np


def _plana(n, rng):
    """Diamante 4D plano por rechazo y relación de Minkowski (solo numpy)."""
    pts = np.empty((0, 4))
    while pts.shape[0] < n:
        c = rng.uniform(-1.0, 1.0, size=(max(64, 4 * n), 4))
        ok = c[np.sqrt((c[:, 1:] ** 2).sum(axis=1)) <= 1.0 - np.abs(c[:, 0])]
        pts = np.vstack([pts, ok])
    x = pts[:n]
    dt = x[None, :, 0] - x[:, None, 0]
    dr = np.sqrt(((x[None, :, 1:] - x[:, None, 1:]) ** 2).sum(axis=2))
    return (dt > 0) & (dt >= dr)


FAMILIAS = [
    {"nombre": "F_gemelos", "tipo": "se_desvanece", "estructuras": ["G1", "G2"]},
    {"nombre": "F_region", "tipo": "autosemejante", "estructuras": ["R1"]},
]
ESTRUCTURAS = {
    "G1": {"tipo": "copias", "k": 1, "frac": 0.2},
    "G2": {"tipo": "copias", "k": 2, "frac": 0.2},
    "R1": {"tipo": "region", "generar": _plana},
}
IDENTIFICABLES = {
    "G1": {"identificable": True, "medidas": {}, "semillas": None, "N": 8000},
    "G2": {"identificable": False, "medidas": {}, "semillas": None, "N": 8000},
    "R1": {"identificable": True, "medidas": {}, "semillas": None, "N": 8000},
}
'''


def _m() -> ModuleType:
    tools = str(ROOT / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("r20")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(880002, spawn_key=(seed,)))


def _dummy(tmp_path: Path, text: str = DUMMY) -> Path:
    p = tmp_path / "sellado_dummy.py"
    p.write_text(text)
    return p


# --- Radio conforme y regla de cobertura -------------------------------------------------------------------


def test_conformal_rank_and_ties() -> None:
    m = _m()
    assert m.conformal_k(100) == 96 and m.conformal_k(19) == 19 and m.conformal_k(20) == 20
    d = np.arange(100, dtype=float)[::-1].copy()  # el 96.º menor es 95
    r, k, j = m.conformal_radius(d)
    assert (r, k) == (95.0, 96) and d[j] == 95.0
    # empates en el corte: se ordenan por índice; el valor no cambia y el índice elegido es determinista
    d = np.concatenate([np.full(50, 1.0), np.full(50, 2.0)])
    d[[3, 70, 90]] = 2.0  # 47 unos y 53 doses
    r, k, j = m.conformal_radius(d)
    order = [i for i in range(100) if d[i] == 1.0] + [i for i in range(100) if d[i] == 2.0]
    assert r == 2.0 and j == order[95]
    with pytest.raises(ValueError):
        m.conformal_radius(np.arange(18, dtype=float))  # k = 19 > 18
    assert m.min_validas(100) == 87 and m.min_validas(20) == 18


def test_conformal_coverage_exchangeable() -> None:
    """Cobertura marginal de un punto nuevo intercambiable = 96/101 y riesgo del gate 87/100 ≈ 0.0119 (Beta-Binomial)."""
    m = _m()
    r = _rng(1)
    reps = 20000
    cal = r.random((reps, 100))
    radius = np.sort(cal, axis=1)[:, m.conformal_k(100) - 1]
    new = r.random(reps)
    cov = float(np.mean(new <= radius))
    assert abs(cov - 96 / 101) < 0.006
    val = r.random((reps, 100))
    cnt = np.sum(val <= radius[:, None], axis=1)
    assert abs(float(np.mean(cnt < 87)) - 0.0119) < 0.004
    # la función del ejecutor da el mismo radio que el orden estadístico
    for i in range(5):
        assert m.conformal_radius(cal[i])[0] == radius[i]


def test_conformal_band_on_vectors() -> None:
    """Banda completa (Ledoit-Wolf con las de ajuste, radio con las de calibración) sobre vectores gaussianos de 14
    componentes: la cobertura media de validación ≈ 96/101; una componente con nan se excluye de todo."""
    m = _m()
    r = _rng(2)
    p = len(m.r14b.COMPONENTS)
    a = r.normal(size=(p, p)) / math.sqrt(p)
    covs = []
    for _ in range(150):
        fit, cal, val = (r.normal(size=(100, p)) @ a + 3.0 for _ in range(3))
        b = m.ConformalBand(fit, cal, val)
        covs.append(b.cubiertas / 100)
        assert b.ok == (b.cubiertas >= 87) and b.k == 96
    assert abs(float(np.mean(covs)) - 96 / 101) < 0.012
    fit, cal, val = (r.normal(size=(30, p)) + 3.0 for _ in range(3))
    val[4, 2] = np.nan  # r_10 no medible en una validación: la componente sale para todos
    b = m.ConformalBand(fit, cal, val)
    assert 2 not in b.comps and len(b.comps) == p - 1 and b.n_val == 30 and b.minimo == math.ceil(0.87 * 30)
    v = fit[0].copy()
    assert math.isfinite(b.D(v))
    v[2] = np.nan
    assert math.isfinite(b.D(v))  # componente fuera de la banda: no importa
    v[5] = np.nan
    assert math.isnan(b.D(v))  # componente de la banda no medible: D no definido


# --- Control emparejado por N -----------------------------------------------------------------------------


@pytest.mark.parametrize("k", [1, 2])
def test_paired_copias_sizes(k: int) -> None:
    m = _m()
    n0, frac = 400, 0.2
    _, base, _ = m.r19.copias_generar(n0, k, frac, _rng(10), _rng(11))
    assert base.shape[0] == n0 - k * round(frac * n0 / k)
    for n in (n0, 233, 201, 101, 50, 49):
        cc, info = m.copias_control_n(base, n, k, frac, _rng(100 + n), _rng(200 + n))
        nc = k * round(frac * n / k)
        assert cc.shape == (n, n) and info["n_copias"] == nc and info["n_base"] == n - nc
        assert int(info["es_copia"].sum()) == nc and info["copias_verificadas"] == nc
        if n == n0:
            assert info["n_base"] == base.shape[0]  # q = 1: base del mismo tamaño, copias nuevas
    assert m.r19.n_copias_de(1001, 0.10, 2) == (100, 901)
    with pytest.raises(RuntimeError):
        m.copias_control_n(base[:10, :10], 400, k, frac, _rng(1), _rng(2))


def test_paired_region_same_n(tmp_path: Path) -> None:
    m = _m()
    mod = m.load_sellado(_dummy(tmp_path))
    gen = mod.ESTRUCTURAS["R1"]["generar"]
    for n in (1, 37, 120):
        assert m.region_control_n(gen, n, _rng(n)).shape == (n, n)
    with pytest.raises(TypeError):
        m.region_control_n(lambda n, r: np.zeros((n, n)), 10, _rng(1))


# --- Valores no finitos ----------------------------------------------------------------------------------------


def test_nonfinite_L_median_bootstrap() -> None:
    m = _m()
    lr = m.log_ratio(np.array([0.0, 1.0, 0.0, np.nan, 2.0]), np.array([1.0, 0.0, 0.0, 1.0, 1.0]))
    assert lr[0] == -np.inf and lr[1] == np.inf and np.isnan(lr[2]) and np.isnan(lr[3]) and lr[4] == pytest.approx(math.log(2))
    # cuantiles extendidos: iguales a numpy con datos finitos
    x = _rng(3).normal(size=31)
    for p in (0.025, 0.5, 0.975):
        assert m.ext_quantile_sorted(np.sort(x), p)[0] == pytest.approx(np.percentile(x, p * 100))
    assert m.ext_median(np.array([-np.inf, -np.inf, 1.0, 2.0])) == -np.inf
    assert m.ext_median(np.array([-np.inf, 1.0, 2.0, np.inf])) == 1.5
    assert m.ext_median(np.array([1.0, np.inf, np.inf])) == np.inf
    assert math.isnan(m.ext_median(np.array([-np.inf, np.inf])))
    assert math.isnan(m.ext_median(np.array([1.0, np.nan, 2.0])))
    # 30 realizaciones con 12 en −∞: mediana finita, IC con −∞ abajo, orden respetado
    a = np.concatenate([np.full(12, -np.inf), -2.0 + 0.1 * _rng(4).normal(size=18)])
    idx = m.boot_indices(30, _rng(5))
    med, lo, hi = m.median_ci_ext(a, idx)
    assert math.isfinite(med) and lo == -np.inf and lo <= med <= hi < -m.MARGEN
    assert m.evaluable(med, lo, hi) and m.clase("se_desvanece", lo, hi) == "separa"
    assert m.clase("autosemejante", lo, hi) == "fuera" and m.clase("autosemejante", -np.inf, 0.0) == "ambiguo"
    # +∞ en la mayoría: mediana +∞
    b = np.concatenate([np.full(20, np.inf), np.zeros(10)])
    med, lo, hi = m.median_ci_ext(b, idx)
    assert med == np.inf and hi == np.inf and m.clase("se_desvanece", lo, hi) == "contrario"
    # nan (0/0 o no medible) en una sola realización: no evaluable, sin excluirla
    c = a.copy()
    c[7] = np.nan
    med, lo, hi = m.median_ci_ext(c, idx)
    assert not m.evaluable(med, lo, hi) and m.clase("se_desvanece", lo, hi) == "no evaluable"
    # los mismos índices para todos los q (realizaciones completas)
    assert idx.shape == (m.N_BOOT, 30)


# --- Veredicto y población -------------------------------------------------------------------------------------


def test_verdict_and_population() -> None:
    m = _m()
    fams = [{"nombre": "A", "tipo": "se_desvanece", "estructuras": ["a1", "a2", "a3"]},
            {"nombre": "B", "tipo": "se_desvanece", "estructuras": ["b1", "b2", "b3"]},
            {"nombre": "C", "tipo": "autosemejante", "estructuras": ["c1", "c2", "c3"]},
            {"nombre": "D", "tipo": "autosemejante", "estructuras": ["d1", "d2", "d3"]}]
    ident = {s: True for f in fams for s in f["estructuras"]}
    assert m.poblacion(fams, ident)["ok"]
    ident2 = {**ident, "d1": False, "d2": False}  # D queda con 1 identificable
    pob = m.poblacion(fams, ident2)
    assert not pob["ok"] and pob["familias_suficientes"]["autosemejante"] == ["C"]
    assert m.poblacion(fams, {**ident, "a3": False, "b1": False})["ok"]  # 2 por familia bastan
    pasa = {"a1": ("se_desvanece", -1.8, -0.9, True), "b1": ("se_desvanece", -np.inf, -0.4, True),
            "c1": ("autosemejante", -0.2, 0.1, True), "d1": ("autosemejante", -0.1, 0.2, True)}
    assert m.verdict(True, True, pasa)[0] == "PASA"
    assert m.verdict(False, True, pasa) == ("INDETERMINADO", "calibración")
    assert m.verdict(True, False, pasa) == ("INDETERMINADO", "población insuficiente")
    assert m.verdict(True, True, {})[0] == "INDETERMINADO"
    assert m.verdict(True, True, {**pasa, "b1": ("se_desvanece", -0.2, 0.3, True)})[0] == "FALLA"
    assert m.verdict(True, True, {**pasa, "c1": ("autosemejante", 0.4, np.inf, True)})[0] == "FALLA"
    assert m.verdict(True, True, {**pasa, "c1": ("autosemejante", -0.9, -0.4, True)})[0] == "FALLA"
    assert m.verdict(True, True, {**pasa, "b1": ("se_desvanece", -0.6, 0.0, True)})[0] == "INDETERMINADO"
    # votante no evaluable: nunca PASA; FALLA si otra contradice
    ne = {**pasa, "a1": ("se_desvanece", math.nan, math.nan, False)}
    assert m.verdict(True, True, ne) == ("INDETERMINADO", "alguna votante no evaluable")
    assert m.verdict(True, True, {**ne, "c1": ("autosemejante", 0.5, 0.9, True)})[0] == "FALLA"
    # una no evaluable con IC finito guardado tampoco cuenta
    assert m.verdict(True, True, {**pasa, "a1": ("se_desvanece", -2.0, -1.0, False)})[0] == "INDETERMINADO"


# --- Sellado e interfaz ----------------------------------------------------------------------------------------


def test_validate_and_structures(tmp_path: Path) -> None:
    m = _m()
    mod = m.load_sellado(_dummy(tmp_path))
    warn = m.validate_sellado(mod)
    assert any("familias (§3" in w for w in warn) and any("candidatas" in w for w in warn)
    sts = {s["nombre"]: s for s in m.structures(mod)}
    assert sts["G1"]["vota"] and sts["R1"]["vota"] and not sts["G2"]["vota"] and sts["G2"]["estado"] == m.NO_IDENT
    assert sts["G1"]["idx"] == 100 and sts["G2"]["idx"] == 101 and sts["G2"]["k"] == 2  # índice de semilla aparte de k
    assert sts["G1"]["ent_puntos"] == m.ENT_SELLADO == 20261115 and sts["plana"]["ent_puntos"] == m.ENT_EST == 20261113
    assert not sts["D5_gemelos"]["vota"] and sts["plana"]["estado"] == "calibracion"
    assert m.ENT_REF == 20261116 and m.EXPERIMENT_ID == 333
    bad = DUMMY.replace('    "R1": {"identificable": True, "medidas": {}, "semillas": None, "N": 8000},\n', "")
    with pytest.raises(SystemExit):
        m.validate_sellado(m.load_sellado(_dummy(tmp_path, bad)))
    bad2 = DUMMY.split("IDENTIFICABLES")[0]
    with pytest.raises(SystemExit):
        m.validate_sellado(m.load_sellado(_dummy(tmp_path, bad2)))


def test_end_to_end_small(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    m = _m()
    monkeypatch.chdir(ROOT)
    out = tmp_path / "out"
    sel = _dummy(tmp_path)
    args = argparse.Namespace(part="all", struct=None, real=None, s=None, sellado=str(sel), workers=1, budget_min=None,
                              n0=200, reals=3, nref=19, out=str(out), allow_dirty=True)
    m.cmd_run(args)
    m.cmd_run(args)  # reanuda: nada por correr
    m.cmd_summarize(argparse.Namespace(out=str(out), sellado=None))
    s = json.loads((out / "summary.json").read_text())
    assert s["veredicto"] == "INDETERMINADO"  # una familia por tipo: población insuficiente o calibración
    assert s["votan"] == ["G1", "R1"] and s["no_identificables"] == ["G2"] and "D5_gemelos" in s["no_votan"]
    assert s["familias"]["F_gemelos"]["excluidas_no_identificables"] == ["G2"]
    assert not s["poblacion"]["ok"]
    for qi in range(4):
        b = s["bandas"][str(qi)]
        assert b["k"] == 19 and b["validacion"]["n"] == 19 and b["validacion"]["minimo"] == 17
    t = s["estructuras"]["G2"]
    assert len(t["L_q"]) == 3 and len(t["por_q"]) == 4 and t["estado"] == m.NO_IDENT
    # control emparejado por N en cada (realización, q)
    for name in ("G1", "G2", "R1", "plana"):
        for real in range(3):
            d = json.loads(m.struct_path(out, name, real).read_text())
            for row in d["rows"]:
                assert row["obs"]["N"] == row["ctrl"]["N"] == row["N"]
                if row["ctrl_defecto"]:
                    k = 1 if name == "G1" else 2
                    assert row["ctrl_defecto"]["n_copias"] == k * round(0.2 * row["N"] / k)
    assert s["sellado_sha256"] == m.r19.sha256(sel) and s["manifest"]["n0"] == 200
    assert (out / "entorno.json").exists() and any((out / "logs").glob("run_*.log"))
    assert (out / "tabla_familias.csv").exists()
    m.cmd_figure(argparse.Namespace(out=str(out), fig=str(tmp_path / "r20.png")))
    assert (tmp_path / "r20.png").stat().st_size > 10000
