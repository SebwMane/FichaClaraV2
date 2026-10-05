"""Pruebas ligeras de las funciones puras de tools/o04b_lite.py y tools/l3b_stability.py (bloque L, Consejo Rev. 2)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from omega.config.seeds import SeedKey, make_rng
from omega.io.passport import array_digest
from omega.types import FloatArray
from tools import l3b_stability as L
from tools import o04b_lite as O


def _path_adj(n: int) -> FloatArray:
    a = np.zeros((n, n), dtype=np.float64)
    for i in range(n - 1):
        a[i, i + 1] = a[i + 1, i] = 1.0
    return a


def _sym(n: int, key: int) -> FloatArray:
    rng = make_rng(SeedKey(20261005, (9, key)))
    v = np.triu(rng.random((n, n)), k=1)
    return np.asarray(v + v.T, dtype=np.float64)


# ------------------------------------------------------------------ L-3b: observables


def test_r_orig_known_cases() -> None:
    geom = L.make_geometry(_path_adj(4))
    w = np.zeros((4, 4))
    w[0, 1] = w[1, 0] = 0.5   # arista del input
    w[0, 3] = w[3, 0] = 0.5   # no arista
    assert L.r_orig(w, geom.edge_upper) == pytest.approx(0.5)
    assert L.r_orig(_path_adj(4), geom.edge_upper) == pytest.approx(1.0)
    assert L.r_orig(np.zeros((4, 4)), geom.edge_upper) == 0.0  # suma 0 -> 0


def test_rho_spearman_known_cases() -> None:
    n = 6
    a = _path_adj(n)
    geom = L.make_geometry(a)
    d = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :]).astype(np.float64)
    w_dec = np.where(d > 0, 0.9 ** d, 0.0)
    w_inc = np.where(d > 0, 0.1 * d, 0.0)
    assert L.rho_spearman(w_dec, geom) == pytest.approx(1.0)
    assert L.rho_spearman(w_inc, geom) == pytest.approx(-1.0)
    const = np.ones((n, n)) - np.eye(n)
    assert math.isnan(L.rho_spearman(const, geom))  # W constante -> NaN (cuenta como fallo)
    # la distancia 5 (extremos) queda fuera: solo cuentan los pares con d0 <= 4
    assert int(geom.pair_mask.sum()) == 15 - 1


def test_rho_spearman_binary_input_is_below_threshold() -> None:
    """Hecho medido (informe de L-3b): rho_S(A, -d0) de un input binario sobre todos los pares d0<=4 queda < 0.5."""
    a, geom = L.cached_input("T3", 0)
    assert 0.3 < L.rho_spearman(a, geom) < 0.5


def test_giant_fraction() -> None:
    w = np.zeros((6, 6))
    for i, j in ((0, 1), (1, 2), (3, 4)):
        w[i, j] = w[j, i] = 1.0
    w[2, 3] = w[3, 2] = 0.05  # < 0.1 max -> no cuenta
    assert L.giant_fraction(w) == pytest.approx(3 / 6)
    assert L.giant_fraction(np.zeros((6, 6))) == 0.0


@pytest.mark.parametrize(
    ("r", "rs", "giant", "cls", "expected"),
    [
        (0.6, 0.6, 0.6, "otro", True),
        (0.5, 0.5, 0.5, "otro", True),            # umbrales inclusivos
        (0.49, 0.9, 0.9, "otro", False),          # R_orig
        (0.9, 0.49, 0.9, "otro", False),          # rho_S
        (0.9, 0.9, 0.49, "otro", False),          # gigante
        (0.9, float("nan"), 0.9, "otro", False),  # NaN falla
        (0.9, 0.9, 0.9, "cliques_solapadas", False),
        (0.9, 0.9, 0.9, "clique_única", False),
        (0.9, 0.9, 0.9, "multi_clique", False),
        (0.9, 0.9, 0.9, "vacío", False),
        (0.9, 0.9, 0.9, "uniforme", False),
    ],
)
def test_geometric_persistent_rule(r: float, rs: float, giant: float, cls: str, expected: bool) -> None:
    assert L.is_geometric_persistent(r, rs, giant, cls) is expected


def test_noisy_start_symmetric_clipped_and_deterministic() -> None:
    a = L.build_input("T3", 0)
    w1 = L.noisy_start(a, 1e-2, make_rng(SeedKey(1, (2,))))
    w2 = L.noisy_start(a, 1e-2, make_rng(SeedKey(1, (2,))))
    assert np.array_equal(w1, w2) and np.array_equal(w1, w1.T)
    assert float(w1.min()) >= 0.0 and float(w1.max()) <= 1.0 and float(np.abs(np.diag(w1)).max()) == 0.0
    assert np.array_equal(L.noisy_start(a, 0.0, make_rng(SeedKey(1, (2,)))), a)


def test_grid_and_task_counts() -> None:
    cfg = L.build_config(False)
    cells = L.cell_grid(cfg)
    assert sum(1 for c in cells if c["block"] == "s0") == 28 and sum(1 for c in cells if c["block"] == "omega_b") == 20
    assert [c["cell_idx"] for c in cells] == list(range(48))
    tasks = L.build_tasks(cfg)
    assert len(tasks) == 5 * 48 * 7 and len({t["id"] for t in tasks}) == len(tasks)
    smoke = L.build_tasks(L.build_config(True))
    assert len(smoke) == 2 * 4 * 1


def test_inputs_have_expected_structure() -> None:
    assert L.build_input("T3", 0).sum() / 2 == 648
    assert L.build_input("T3_decorated", 0).sum() / 2 == 1404
    assert L.build_input("K8_union", 0).sum() / 2 == 27 * 28
    assert L.build_input("ER", 1).sum() == L.build_input("RGG3", 1).sum()  # misma m


# ------------------------------------------------------------------ L-3b: decision global


def _row(inp: str, eps_idx: int, s: int, geo: bool, stable: bool | None = None, *, cell: int = 0, status: str = "CONVERGED",
         drift: float | None = None) -> dict[str, Any]:
    return {"id": f"{inp}|s0|c{cell}|e{eps_idx}|s{s}", "input": inp, "block": "s0", "cell_idx": cell, "eps_idx": eps_idx, "s": s,
            "geo": geo, "stable": (geo if stable is None else stable) if geo else None, "status": status,
            "max_steps_final": status == "MAX_STEPS", "drift_R_orig": drift, "kkt_ok": True}


def _cell(inp: str, eps_idx: int, flags: list[bool], **kw: Any) -> list[dict[str, Any]]:
    return [_row(inp, eps_idx, s, f, **kw) for s, f in enumerate(flags)]


def test_global_decision_yes_when_nothing_persists() -> None:
    rows = _cell("T3", 0, [False] * 3) + _cell("T3", 1, [False] * 3) + _cell("RGG3", 0, [False] * 3)
    assert L.global_decision(rows)["decision"] == "SÍ (incompatible)"


@pytest.mark.parametrize("inp", ["T3", "RGG3"])
def test_global_decision_no(inp: str) -> None:
    rows = _cell(inp, 0, [True, True, False]) + _cell("T3_decorated", 0, [True] * 3)
    assert L.global_decision(rows)["decision"] == "NO"


def test_global_decision_no_decorated_only() -> None:
    rows = _cell("T3_decorated", 0, [True, True, False]) + _cell("T3", 0, [True, False, False])
    assert L.global_decision(rows)["decision"] == "NO-decorado"


def test_vote_needs_two_thirds_and_stability() -> None:
    assert L.global_decision(_cell("T3", 0, [True, False, False]))["decision"] == "SÍ (incompatible)"
    rows = [_row("T3", 0, 0, True, stable=False), _row("T3", 0, 1, True, stable=False), _row("T3", 0, 2, True, stable=True)]
    assert L.global_decision(rows)["decision"] == "SÍ (incompatible)"  # geometricos pero inestables


def test_global_decision_metastable_by_eps_1e3() -> None:
    rows = _cell("RGG3", 0, [False] * 3) + _cell("RGG3", 1, [True, True, False])
    g = L.global_decision(rows)
    assert g["decision"] == "METAESTABLE" and len(g["metastable_cells_eps_1e-3_only"]) == 1


def test_global_decision_metastable_by_drift() -> None:
    rows = _cell("T3", 0, [False] * 3)
    rows.append(_row("T3", 0, 3, False, status="MAX_STEPS", drift=-0.02))
    assert L.global_decision(rows)["decision"] == "METAESTABLE"
    rows[-1]["drift_R_orig"] = 0.005
    assert L.global_decision(rows)["decision"] == "SÍ (incompatible)"
    ctrl = [_row("ER", 0, 0, False, status="MAX_STEPS", drift=0.5)]  # los controles no votan la deriva
    assert L.global_decision(_cell("T3", 0, [False] * 3) + ctrl)["decision"] == "SÍ (incompatible)"


def test_global_decision_control_alarm_does_not_change_decision() -> None:
    rows = _cell("T3", 0, [False] * 3) + _cell("ER", 0, [True, False, False])
    g = L.global_decision(rows)
    assert g["decision"] == "SÍ (incompatible)" and g["control_geometric_persistent_alarm"]["alarm"] is True


def test_global_decision_counts_non_kkt_without_vote() -> None:
    rows = _cell("T3", 0, [False] * 3)
    rows[0]["kkt_ok"] = False
    g = L.global_decision(rows)
    assert g["n_converged_not_kkt_no_vote"] == 1 and g["decision"] == "SÍ (incompatible)"


# ------------------------------------------------------------------ O4B-1


def test_verify_npz_digests(tmp_path: Path) -> None:
    arrays = {"w0": _sym(5, 1)[np.triu_indices(5, 1)], "w_final": _sym(5, 2)[np.triu_indices(5, 1)]}
    path = tmp_path / "x.npz"
    np.savez_compressed(path, **arrays)  # type: ignore[arg-type]
    good = {k: array_digest(v) for k, v in arrays.items()}
    ok, why, loaded = O.verify_npz(path, good, O.file_sha256(path))
    assert ok and why == "" and set(loaded) == {"w0", "w_final"}
    bad = dict(good, w_final="0" * 64)
    ok, why, _ = O.verify_npz(path, bad)
    assert not ok and "w_final" in why
    ok, why, _ = O.verify_npz(path, {"w0": good["w0"]})
    assert not ok
    ok, _, _ = O.verify_npz(path, good, "f" * 64)  # sha256 del archivo
    assert not ok
    ok, why, _ = O.verify_npz(tmp_path / "no.npz", good)
    assert not ok and "ausente" in why


def _orow(i: int, cls: str, status: str = "CONVERGED", kkt_ok: bool | None = True, block: str = "omega_b") -> dict[str, Any]:
    return {"run_id": f"r{i}", "block": block, "state_class": cls, "status": status, "kkt_ok": kkt_ok}


def test_o04b_decision() -> None:
    ok_rows = [_orow(0, "uniforme"), _orow(1, "clique_única"), _orow(2, "cliques_solapadas"), _orow(3, "vacío")]
    assert O.decide(ok_rows)["decision"] == "EXITO"
    assert O.decide([*ok_rows, _orow(4, "otro")])["decision"] == "FRACASO"
    assert O.decide([*ok_rows, _orow(5, "multi_clique", kkt_ok=False)])["decision"] == "FRACASO"
    # un no convergido sin KKT no cuenta; un S0 "otro" tampoco vota
    assert O.decide([*ok_rows, _orow(6, "uniforme", status="MAX_STEPS", kkt_ok=None), _orow(7, "otro", block="s0")])["decision"] == "EXITO"
    assert O.decide([_orow(8, "otro", block="s0")])["decision"] == "SIN DATOS Omega-B"
