"""Tests de A-E1' (docs/OMEGA_CONSEJO_REV17_ES0.md §7.2): omega.causality.sj_multicut y los ajustes de
tools/a_e1p_auditoria.py."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest

from omega.causality.entropy import (
    causal_complement,
    interval_balanced,
    pauli_jordan,
    sj_from_eig,
    ssee,
)
from omega.causality.orders import kr_order
from omega.causality.sj_fast import evaluate
from omega.causality.sj_multicut import multicut
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond

TOL = 1e-9
CS = (0.5, 1.0, 2.0)
Bool = npt.NDArray[np.bool_]


def _tool() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "tools" / "a_e1p_auditoria.py"
    spec = importlib.util.spec_from_file_location("a_e1p_auditoria", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["a_e1p_auditoria"] = mod
    spec.loader.exec_module(mod)
    return mod


def _rng(seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(20261103, spawn_key=(328, 99, seed))))


def _orders() -> list[tuple[str, Bool]]:
    return [
        ("geo2d_96", causal_matrix(sprinkle_diamond(96, 2, _rng(0)))),
        ("geo2d_200", causal_matrix(sprinkle_diamond(200, 2, _rng(1)))),
        ("geo3d_128", causal_matrix(sprinkle_diamond(128, 3, _rng(2)))),
        ("geo4d_128", causal_matrix(sprinkle_diamond(128, 4, _rng(3)))),
        ("kr_96", kr_order(96, _rng(4))),
        ("cadena_64", np.triu(np.ones((64, 64), dtype=bool), k=1)),
    ]


def _reference(c: Bool, region: npt.NDArray[np.int64], cm: float) -> float:
    n = c.shape[0]
    mu, u = np.linalg.eigh(pauli_jordan(c))
    idt, w = sj_from_eig(mu, u, cm * math.sqrt(n) / (4 * math.pi))
    return ssee(idt, w, region, True, cm * math.sqrt(region.size) / (4 * math.pi))


@pytest.mark.parametrize("name,c", _orders())
def test_multicut_equals_reference(name: str, c: Bool) -> None:
    res = multicut(c, CS, c_perp=1.0)
    region = interval_balanced(c, c.shape[0] // 8)
    assert np.array_equal(res.region, region)
    assert region.size > 0, name
    n = c.shape[0]
    mu_ref = np.abs(np.linalg.eigvalsh(pauli_jordan(c)))
    for cr in res.cuts:
        assert abs(cr.s - _reference(c, region, cr.c_mult)) <= TOL, (name, cr.c_mult)
        assert cr.n_kept_global == int(np.sum(mu_ref >= cr.c_mult * math.sqrt(n) / (4 * math.pi)))
        assert 0 <= cr.n_kept_region <= region.size
    assert abs(res.cuts[1].s - evaluate(c).s) <= TOL  # mismo valor que el evaluador R17-3 en c = 1
    rp = causal_complement(c, region)
    assert res.n_perp == rp.size
    if rp.size:
        mu, u = np.linalg.eigh(pauli_jordan(c))
        idt, w = sj_from_eig(mu, u, math.sqrt(n) / (4 * math.pi))
        ref = ssee(idt, w, rp, True, math.sqrt(rp.size) / (4 * math.pi))
        assert res.s_perp is not None and abs(res.s_perp - ref) <= TOL


def test_multicut_subset_path_large_n() -> None:
    """N >= SUBSET_MIN_N usa el eigh por subconjunto: mismo valor que la referencia en los tres cortes."""
    c = causal_matrix(sprinkle_diamond(440, 2, _rng(5)))
    res = multicut(c, CS)
    for cr in res.cuts:
        assert abs(cr.s - _reference(c, res.region, cr.c_mult)) <= TOL


def test_no_region() -> None:
    c = np.zeros((64, 64), dtype=bool)
    res = multicut(c, CS)
    assert res.region.size == 0 and res.cuts == []


def test_fit_recovers_known_slope() -> None:
    tool = _tool()
    rng = np.random.default_rng(7)
    rows = []
    for n in tool.NS:
        for _ in range(30):
            nr = int(rng.integers(int(0.95 * n / 8), int(1.05 * n / 8) + 1))
            s = 0.4 * float(tool.xlog(np.array([nr], dtype=float))[0]) + 1.5 + rng.normal(0, 0.01)
            rows.append({"N": n, "has_region": True, "n_R": nr, "S_c1.0": s})
    f = tool.area_fit(rows, tool.NS, 1.0)
    assert abs(f["a"] - 0.4) < 3 * f["se"] + 1e-3
    assert f["se"] < 0.01 and f["r2"] > 0.99
    assert abs(f["da"] - abs(f["a"] - 1 / 3)) < 1e-15
    exact = tool.ols(np.array([1.0, 2.0, 3.0, 4.0]), np.array([2.0, 2.5, 3.0, 3.5]))
    assert abs(exact["a"] - 0.5) < 1e-12 and exact["se"] < 1e-12
    assert tool.ols(np.array([1.0, 1.0, 1.0]), np.array([1.0, 2.0, 3.0]))["a"] is None


def test_power_fit_exponent() -> None:
    tool = _tool()
    rows = [{"N": n, "t_fast_s": 1e-6 * n**2.5} for n in tool.NS]
    pf = tool.power_fit(rows, "t_fast_s")
    assert abs(pf["exponente"] - 2.5) < 1e-9
    assert abs(tool.ms_at(pf, 1000) - 1e-3 * 1000**2.5) < 1e-6 * 1000**2.5
