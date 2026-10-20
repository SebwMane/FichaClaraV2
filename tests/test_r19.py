"""Pruebas rápidas de tools/r19.py (R19-(b)). Semillas propias (880001), nunca la entropía del experimento. Usan un
módulo sellado de prueba escrito en tmp_path (una familia de copias y una de región), nunca el sellado real."""

from __future__ import annotations

import argparse
import importlib
import json
import math
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

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
    {"nombre": "F_tripletes", "tipo": "se_desvanece", "estructuras": ["T_tripletes"]},
    {"nombre": "F_region", "tipo": "autosemejante", "estructuras": ["R_plana"]},
]
ESTRUCTURAS = {
    "T_tripletes": {"tipo": "copias", "k": 2, "frac": 0.2},
    "R_plana": {"tipo": "region", "generar": _plana},
}
'''


def _m() -> ModuleType:
    tools = str(ROOT / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    return importlib.import_module("r19")


def _rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(880001, spawn_key=(seed,)))


def _dummy(tmp_path: Path) -> Path:
    p = tmp_path / "sellado_dummy.py"
    p.write_text(DUMMY)
    return p


def _check_copies(m: ModuleType, cp: np.ndarray, info: dict[str, Any], k: int) -> None:
    """Cada copia tiene el pasado y el futuro de su original (el elemento no copia con el mismo origen), no se
    relaciona con él, y cada original tiene exactamente k copias."""
    origen, es_copia = info["origen"], info["es_copia"]
    assert int(es_copia.sum()) == info["n_copias"]
    orig_label = {int(o): e for e, o in enumerate(origen) if not es_copia[e]}
    per: dict[int, int] = {}
    for e in np.flatnonzero(es_copia):
        o = orig_label[int(origen[e])]
        assert np.array_equal(cp[e], cp[o]) and np.array_equal(cp[:, e], cp[:, o])
        assert not cp[e, o] and not cp[o, e]
        per[int(origen[e])] = per.get(int(origen[e]), 0) + 1
    assert set(per.values()) <= {k} and len(per) == info["n_originales"]
    assert not np.any(cp & cp.T)


def test_copies_generator(tmp_path: Path) -> None:
    m = _m()
    for k, frac, n in ((1, 0.10, 400), (2, 0.20, 401), (3, 0.15, 300)):
        cp, base, info = m.copias_generar(n, k, frac, _rng(1), _rng(2))
        nc = k * round(frac * n / k)
        assert cp.shape == (n, n) and info["n_copias"] == nc and base.shape[0] == n - nc
        _check_copies(m, cp, info, k)
        # pares con igual pasado y futuro: al menos los de los grupos (original + k copias)
        assert m.r14b.twin_pairs(cp) >= info["n_originales"] * (k + 1) * k // 2
    # D5 como copias con k = 1, frac = 0.10: 800 gemelos sobre 7200
    assert m.n_copias_de(8000, 0.10, 1) == (800, 7200)


def test_k2_not_divisible() -> None:
    m = _m()
    # regla fijada: n_originales = round(frac·N/k), n_copias = k·n_originales
    assert m.n_copias_de(1001, 0.10, 2) == (100, 901)  # 100.1/2 = 50.05 -> 50
    assert m.n_copias_de(1015, 0.10, 2) == (102, 913)  # 101.5/2 = 50.75 -> 51
    assert m.n_copias_de(1000, 0.10, 3) == (99, 901)  # 33.3 -> 33
    assert m.n_copias_de(999, 0.25, 1) == (250, 749)
    cp, _, info = m.copias_generar(1015 // 3, 2, 0.10, _rng(3), _rng(4))  # N = 338: 33.8/2 = 16.9 -> 17 originales
    assert info["n_originales"] == 17 and info["n_copias"] == 34 and cp.shape == (338, 338)
    _check_copies(m, cp, info, 2)
    with pytest.raises(ValueError):
        m.add_copies(np.zeros((10, 10), dtype=bool), 3, 2, _rng(5))


def test_control_keeps_fraction_each_q() -> None:
    m = _m()
    n0, k, frac = 800, 2, 0.2
    _, base, _ = m.copias_generar(n0, k, frac, _rng(10), _rng(11))
    for qi, q in enumerate(m.QS):
        sizes = []
        for s in range(40 if q < 1 else 3):
            cc, info = m.copias_control(base, n0, q, k, frac, _rng(100 + 10 * qi + s), _rng(500 + 10 * qi + s))
            n = cc.shape[0]
            sizes.append(n)
            # se cuentan las copias después de construir: el número es k·round(frac·N/k) y la fracción ≈ frac
            _check_copies(m, cc, info, k)
            assert int(info["es_copia"].sum()) == k * round(frac * n / k)
            assert abs(int(info["es_copia"].sum()) / n - frac) <= k / n
            # la base del control es un subconjunto de la base de X con el orden inducido
            nb = info["n_base"]
            assert nb == n - info["n_copias"] and nb <= base.shape[0]
        if q >= 1.0:
            assert sizes == [n0] * 3
        else:
            assert abs(np.mean(sizes) - n0 * q) < 4 * math.sqrt(n0 * q * (1 - q) / len(sizes)) + 1


def test_region_control_is_binomial_and_matches_dilution(tmp_path: Path) -> None:
    m = _m()
    mod = m.load_sellado(_dummy(tmp_path))
    gen = mod.ESTRUCTURAS["R_plana"]["generar"]

    def empty(n: int, rng: np.random.Generator) -> np.ndarray:
        return np.zeros((n, n), dtype=bool)

    n0, q = 400, 0.25
    sizes = np.array([m.region_control(empty, n0, q, _rng(1000 + s), _rng(2000 + s)).shape[0] for s in range(400)])
    assert abs(sizes.mean() - 100) < 2.0
    assert 50 < sizes.var(ddof=1) < 105  # Binomial: 75
    assert m.region_control(empty, n0, 1.0, _rng(1), _rng(2)).shape == (n0, n0)
    # control de un sprinkling ~ Φ_q del sprinkling: fracción de relaciones y tamaño parecidos en promedio
    fo: list[float] = []
    fc: list[float] = []
    no: list[int] = []
    nc: list[int] = []
    for s in range(30):
        _, cd = m.r14b.dilute(gen(n0, _rng(3000 + s)), 0.5, _rng(4000 + s))
        cc = m.region_control(gen, n0, 0.5, _rng(5000 + s), _rng(6000 + s))
        for arr, nn, c in ((fo, no, cd), (fc, nc, cc)):
            n = c.shape[0]
            arr.append(np.count_nonzero(c) / (n * (n - 1) / 2))
            nn.append(n)
    se = math.sqrt(np.var(fo, ddof=1) / 30 + np.var(fc, ddof=1) / 30)
    assert abs(np.mean(fo) - np.mean(fc)) < 4 * se
    assert abs(np.mean(no) - np.mean(nc)) < 15
    with pytest.raises(TypeError):
        m.region_control(lambda n, r: np.zeros((n, n)), 10, 1.0, _rng(1), _rng(2))


def test_L_and_verdict_logic() -> None:
    m = _m()
    g = m.MARGEN
    assert m.clase("se_desvanece", -2.0, -0.5) == "separa"
    assert m.clase("se_desvanece", -0.1, 0.4) == "contrario"
    assert m.clase("se_desvanece", -0.5, 0.0) == "ambiguo"
    assert m.clase("autosemejante", -0.1, 0.2) == "equivalente"
    assert m.clase("autosemejante", 0.5, 0.9) == "fuera"
    assert m.clase("autosemejante", 0.1, 0.5) == "ambiguo"
    assert m.clase("nula", -g, g) == "equivalente"
    pasa = {"a": ("se_desvanece", -1.8, -0.9), "b": ("se_desvanece", -1.2, -0.4), "c": ("autosemejante", -0.2, 0.1)}
    assert m.verdict(True, pasa) == "PASA"
    assert m.verdict(False, pasa) == "INDETERMINADO"  # calibración
    assert m.verdict(True, {}) == "INDETERMINADO"
    assert m.verdict(True, {**pasa, "b": ("se_desvanece", -0.2, 0.3)}) == "FALLA"
    assert m.verdict(True, {**pasa, "c": ("autosemejante", 0.4, 0.8)}) == "FALLA"
    assert m.verdict(True, {**pasa, "c": ("autosemejante", -0.9, -0.4)}) == "FALLA"
    assert m.verdict(True, {**pasa, "b": ("se_desvanece", -0.6, 0.0)}) == "INDETERMINADO"
    assert m.verdict(True, {**pasa, "c": ("autosemejante", 0.1, 0.5)}) == "INDETERMINADO"
    # L sintético: D_obs = D_ctrl·e^{−1.5} con ruido -> IC por debajo de −ln 1.3
    r = _rng(20)
    dctl = np.exp(r.normal(1.0, 0.2, size=20))
    lq = np.log(dctl * np.exp(-1.5 + 0.1 * r.normal(size=20)) / dctl)
    med, lo, hi = m.median_ci(lq, _rng(21), 500)
    assert lo <= med <= hi < -g and abs(med + 1.5) < 0.1
    assert m.min_cobertura(20) == 16 and m.min_cobertura(4) == 4


def test_validate_and_structures(tmp_path: Path) -> None:
    m = _m()
    mod = m.load_sellado(_dummy(tmp_path))
    warn = m.validate_sellado(mod)
    assert any("familias de tipo se_desvanece" in w for w in warn)  # una sola familia de cada tipo
    sts = {s["nombre"]: s for s in m.structures(mod)}
    assert sts["D5_gemelos"]["tipo"] == "copias" and sts["D5_gemelos"]["k"] == 1 and sts["D5_gemelos"]["frac"] == 0.10
    assert sts["T_tripletes"]["grupo"] == "sellada" and sts["T_tripletes"]["ent_puntos"] == m.ENT_SELLADO
    assert sts["plana"]["ent_puntos"] == m.ENT_EST and sts["R_plana"]["k"] == 101
    bad = ModuleType("bad")
    bad.FAMILIAS = [{"nombre": "x", "tipo": "se_desvanece", "estructuras": ["nope"]}]  # type: ignore[attr-defined]
    bad.ESTRUCTURAS = {}  # type: ignore[attr-defined]
    with pytest.raises(SystemExit):
        m.validate_sellado(bad)


def test_end_to_end_small(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    m = _m()
    monkeypatch.chdir(ROOT)
    out = tmp_path / "out"
    sel = _dummy(tmp_path)
    args = argparse.Namespace(part="all", struct=None, real=None, s=None, sellado=str(sel), workers=1, budget_min=None,
                              n0=240, reals=3, nref=6, nnuevas=4, out=str(out), allow_dirty=True)
    m.cmd_run(args)
    m.cmd_run(args)  # reanuda: nada por correr
    m.cmd_summarize(argparse.Namespace(out=str(out), sellado=None))
    s = json.loads((out / "summary.json").read_text())
    assert s["veredicto"] in ("PASA", "FALLA", "INDETERMINADO")
    assert s["votan"] == ["R_plana", "T_tripletes"] and "D5_gemelos" in s["no_votan"]
    t = s["estructuras"]["T_tripletes"]
    assert len(t["L_q"]) == 3 and len(t["por_q"]) == 4
    assert abs(t["por_q"][3]["frac_copias_ctrl_media"] - 0.2) < 0.1
    assert s["sellado_sha256"] == m.sha256(sel) and s["manifest"]["n0"] == 240
