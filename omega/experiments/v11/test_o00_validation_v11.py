"""Omega-0 (Ω-1.1): validacion de las medidas nuevas sobre grafos de referencia (DESIGN §4, tablas de §1.1-§1.6).

Las expectativas se escriben en `summary.json` (stage "preregistered") ANTES de medir y no se tocan despues
(M§44); los resultados se anaden en el stage "final". Un valor que no se reproduce se informa, no se ajusta.
Smoke: anillo 100, T^2 10x10, T^3 6^3, K_10, octaedro, arbol pequeno, vacio. Full: tablas completas.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.certificate.taxonomy import assess_run
from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings11 import DistanceMode, Omega11Config
from omega.controls.erdos_renyi import erdos_renyi_gnp
from omega.controls.random_geometric import balanced_tree, rgg_torus
from omega.controls.small_world import watts_strogatz
from omega.experiments.reference_graphs import complete_graph, periodic_lattice
from omega.experiments.v11.gate import PassportWriter, check_expectation, derive_cfg, runs_root, write_summary
from omega.geometry.weyl import weyl_dimension
from omega.phases.scan import default_config
from omega.types import FloatArray, RunStatus

NAME = "o00_validation"
EXPERIMENT_ID = 1100
ENTRYPOINT = "omega.experiments.v11.test_o00_validation_v11:run"

Builder = Callable[[np.random.Generator], FloatArray]


def _octahedron() -> FloatArray:
    w = np.ones((6, 6)) - np.eye(6)
    for i, j in ((0, 1), (2, 3), (4, 5)):
        w[i, j] = w[j, i] = 0.0
    return np.asarray(w, dtype=np.float64)


Expect = tuple[str, str, Any, float | None]


def _ring(n: int) -> Builder:
    return lambda r: periodic_lattice((n,))


def _torus(shape: tuple[int, ...]) -> Builder:
    return lambda r: periodic_lattice(shape)


def _case(name: str, build: Builder, expect: list[Expect], *, weyl_only: bool = False,
          w_min: float | None = None, reps: int = 1) -> dict[str, Any]:
    return {"name": name, "build": build, "expect": expect, "weyl_only": weyl_only, "w_min": w_min, "reps": reps}


def _kappa_zero() -> Expect:
    return ("kappa_max_abs", "lt", 1e-9, None)


def _cases(mode: Literal["smoke", "full"]) -> list[dict[str, Any]]:
    """Casos y expectativas preregistradas: (metrica, op, objetivo, tolerancia)."""
    ring, torus = _ring, _torus
    tol_torus_local: list[Expect] = [("detour_fraction", "abs_within", 1.0, 1e-9)]
    if mode == "smoke":
        return [
            _case("ring_100", ring(100), [("d_weyl", "abs_within", 1.0, 0.25), _kappa_zero()]),
            _case("torus2_10x10", torus((10, 10)), [("d_weyl", "abs_within", 2.0, 0.6), _kappa_zero(), *tol_torus_local]),
            _case("torus3_6x6x6", torus((6, 6, 6)), [("d_weyl", "abs_within", 3.0, 0.8), _kappa_zero(), *tol_torus_local]),
            _case("complete_10", lambda r: complete_graph(10), [
                ("kappa_kn_error", "lt", 1e-9, None), ("d_weyl_status", "eq", "no_window", None),
                ("primary", "eq", "Ω-F1", None)]),
            _case("octahedron", lambda r: _octahedron(), [("betti", "eq", [1, 0, 1], None)]),
            _case("tree3_h3", lambda r: balanced_tree(3, 3), [("passes", "eq", False, None)]),
            _case("empty_30", lambda r: np.zeros((30, 30)), [("primary", "eq", "Ω-F0", None)]),
        ]
    return [
        _case("ring_300", ring(300), [("d_weyl", "abs_within", 1.01, 0.05), _kappa_zero()]),
        _case("ring_1000", ring(1000), [("d_weyl", "abs_within", 1.01, 0.05)], weyl_only=True),
        _case("torus2_17", torus((17, 17)), [("d_weyl", "abs_within", 2.13, 0.2), _kappa_zero(), *tol_torus_local]),
        _case("torus2_28", torus((28, 28)), [("d_weyl", "abs_within", 2.08, 0.2)], weyl_only=True),
        _case("torus2_32", torus((32, 32)), [("d_weyl", "abs_within", 2.07, 0.2)], weyl_only=True),
        _case("torus3_7", torus((7, 7, 7)), [("d_weyl", "abs_within", 2.88, 0.3), _kappa_zero(), *tol_torus_local]),
        _case("torus3_9", torus((9, 9, 9)), [("d_weyl", "abs_within", 2.98, 0.25), _kappa_zero(),
                                              ("homogeneity_cv", "abs_within", 0.0, 0.02), *tol_torus_local]),
        _case("torus3_10", torus((10, 10, 10)), [("d_weyl", "abs_within", 3.08, 0.25)], weyl_only=True),
        _case("torus3_15", torus((15, 15, 15)), [("d_weyl", "abs_within", 3.15, 0.25)], weyl_only=True),
        _case("rgg3_800_k12", lambda r: rgg_torus(800, 3, 12.0, r), [
            ("d_weyl", "abs_within", 2.96, 0.25), ("passes", "eq", True, None)], reps=3),
        _case("rgg2_800_k10", lambda r: rgg_torus(800, 2, 10.0, r), [("d_weyl", "abs_within", 2.0, 0.25)], weyl_only=True),
        _case("ws_800_k12_p0.05", lambda r: watts_strogatz(800, 12, 0.05, r), [("primary", "eq", "Ω-F3", None)]),
        _case("er_800_k12", lambda r: erdos_renyi_gnp(800, 12.0 / 799.0, r), [("codes", "contains", "Ω-F4", None)]),
        _case("tree3_h5", lambda r: balanced_tree(3, 5), [("passes", "eq", False, None)]),
        _case("slab_32x5x5", lambda r: periodic_lattice((32, 5, 5)), [("passes", "eq", False, None)]),
        _case("complete_100", lambda r: complete_graph(100), [
            ("kappa_kn_error", "lt", 1e-9, None), ("d_weyl_status", "eq", "no_window", None), ("primary", "eq", "Ω-F1", None)]),
        _case("empty_100", lambda r: np.zeros((100, 100)), [("primary", "eq", "Ω-F0", None)]),
        _case("uniform_300_wmin0.1", lambda r: r.random((300, 300)), [("d_weyl", "gt", 5.0, None)], weyl_only=True, w_min=0.1),
        _case("uniform_300_wmin0.5", lambda r: r.random((300, 300)), [("d_weyl", "gt", 5.0, None)], weyl_only=True, w_min=0.5),
        _case("uniform_300_wmin0.9", lambda r: r.random((300, 300)), [("d_weyl", "gt", 5.0, None)], weyl_only=True, w_min=0.9),
    ]


def _sym(w: FloatArray) -> FloatArray:
    u = np.triu(w, 1)
    return np.asarray(u + u.T, dtype=np.float64)


def _num(x: float) -> float | None:
    return float(x) if np.isfinite(x) else None


def _measure(w: FloatArray, cfg: Omega11Config, key: SeedKey, weyl_only: bool) -> dict[str, Any]:
    if weyl_only:
        d = weyl_dimension(w, cfg.base.graph.w_min, cfg.weyl)
        return {"n": int(w.shape[0]), "d_weyl": _num(d.value), "d_weyl_status": d.status, "d_weyl_plateau": bool(d.plateau)}
    cfg = evidence_cfg(cfg, w)  # presupuesto de coste solo en estados F1 garantizados (auditoria B1)
    ev = collect_run_evidence(w, RunStatus.CONVERGED, cfg, evidence_rng(key))
    a = assess_run(ev, cfg)
    kv = ev.curvature.edge_values
    alpha = cfg.curvature.idleness
    n = ev.n
    kn = 1.0 - abs(alpha - (1.0 - alpha) / (n - 1))
    hop = ev.suite.estimates.get(DistanceMode.HOP)
    return {
        "n": n,
        "d_weyl": _num(ev.d_weyl.value),
        "d_weyl_status": ev.d_weyl.status,
        "d_weyl_plateau": bool(ev.d_weyl.plateau),
        "d_vol": _num(ev.geometry.d_eff.value),
        "d_s": _num(ev.geometry.d_s.value),
        "d_hop": None if hop is None else _num(hop.value),
        "zeta": _num(ev.suite.resistance_exponent),
        "distance_sensitive": bool(ev.suite.distance_sensitive),
        "kappa_mean": _num(ev.curvature.mean),
        "kappa_max_abs": float(np.max(np.abs(kv))) if kv.size else None,
        "kappa_kn_error": float(np.max(np.abs(kv - kn))) if kv.size else None,
        "betti": [int(b) for b in ev.topo.clique_at_w_min.betti],
        "homogeneity_cv": _num(ev.homogeneity.cv),
        "isotropy_median": _num(ev.isotropy.median_ratio),
        "detour_fraction": _num(ev.locality.detour_fraction),
        "annulus_ok": bool(ev.annulus.ok),
        "b1_density_w_min": _num(ev.topo.at_w_min.b1_density),
        "codes": [c.value for c in a.codes],
        "primary": None if a.primary is None else a.primary.value,
        "passes": bool(a.passes),
    }


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-00 y devuelve la ruta de `summary.json` (sin prerrequisitos: es el primer paso)."""
    cases = _cases(mode)
    cfg = derive_cfg(cfg, experiment_id=EXPERIMENT_ID, replicates=max(3, cfg.base.seeds.replicates))
    prereg = [
        {"id": f"{c['name']}.{m}", "case": c["name"], "metric": m, "op": op, "target": t, "tol": tol}
        for c in cases for (m, op, t, tol) in c["expect"]
    ]
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID,
        "description": "validacion Omega-0 de las medidas nuevas (DESIGN §4: tablas §1.1-§1.6)",
        "note": (
            "ER k12: el diseno dice 'da F4'; con la precedencia F3 > F4 el codigo primario es F3 (la ER es no local), "
            "asi que se preregistra F4 en el conjunto de codigos."
        ),
        "expected": prereg,
        "expectations_source": "OMEGA_1_1_DESIGN §1.3 (Weyl), §1.6 (kappa), §1.2 (beta), §4 (O-00)",
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    table: list[dict[str, Any]] = []
    outcomes: dict[str, dict[str, Any]] = {}
    for idx, c in enumerate(cases):
        rows: list[dict[str, Any]] = []
        for rep in range(int(c["reps"])):
            key = seed_key(cfg.base.seeds, idx, rep)
            w = c["build"](make_rng(key))
            if c["name"].startswith("uniform"):
                w = _sym(w)
            c_cfg = derive_cfg(cfg, n=int(w.shape[0]), w_min=c["w_min"])
            m = _measure(w, c_cfg, key, bool(c["weyl_only"]))
            m["case"], m["rep"] = c["name"], rep
            if rep == 0 and not c["weyl_only"]:
                m["passport"] = writer.save(c_cfg, seed=key, w=w, termination="STATIC",
                                            init_distribution=f"reference/{c['name']}", results=m)
            rows.append(m)
            table.append(m)
        for metric, op, target, tol in c["expect"]:
            met = [check_expectation(op, r.get(metric), target, tol) for r in rows]
            outcomes[f"{c['name']}.{metric}"] = {
                "measured": [r.get(metric) for r in rows], "met": bool(all(met)), "n_reps": len(rows)}

    summary = {
        **header,
        "stage": "final",
        "table": table,
        "expectation_results": outcomes,
        "expectations_met": all(v["met"] for v in outcomes.values()),
        "n_expectations": len(outcomes),
        "n_expectations_met": sum(1 for v in outcomes.values() if v["met"]),
        "passports": writer.labels,
        "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode)


def _small_cfg(n: int = 24) -> Omega11Config:
    return Omega11Config(base=default_config(n, master_entropy=20240901, replicates=3))


def test_smoke(tmp_path: Path) -> None:
    path = run(_small_cfg(), tmp_path, mode="smoke")
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["stage"] == "final" and data["complete"] is True
    assert len(data["expected"]) >= 10
    res = data["expectation_results"]
    # Expectativas analiticas (exactas): kappa=0 en retículas, K_n, beta del octaedro, K_n -> F1, vacio -> F0.
    for key in ("ring_100.kappa_max_abs", "torus2_10x10.kappa_max_abs", "torus3_6x6x6.kappa_max_abs",
                "complete_10.kappa_kn_error", "complete_10.primary", "octahedron.betti", "empty_30.primary",
                "tree3_h3.passes"):
        assert res[key]["met"] is True, (key, res[key])
    assert len(data["passports"]) == 7 and all(p.startswith("Ω-EXP-") for p in data["passports"])
    assert len(list((tmp_path / NAME / "passports").glob("OMEGA-EXP-*.json"))) == 7


def test_smoke_preregisters_before_results(tmp_path: Path) -> None:
    """El stage preregistrado se escribe antes que los resultados (M§44): la expectativa no depende de la medida."""
    path = run(_small_cfg(), tmp_path, mode="smoke")
    final = json.loads(path.read_text(encoding="utf-8"))
    ids = [e["id"] for e in final["expected"]]
    assert len(ids) == len(set(ids))
    assert set(final["expectation_results"]) == set(ids)


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(200, master_entropy=20240901, replicates=3))
    data = json.loads(run(cfg, runs_root(), mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
