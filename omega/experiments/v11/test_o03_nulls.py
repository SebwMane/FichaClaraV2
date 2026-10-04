"""Omega-3 (Ω-1.1): bateria de nulos y tasa de falsos positivos por campo (DESIGN §1.10, §4).

Familias: ER (k 6/12/30), nulos emparejados del RGG3 (CM, recableado, ER, pesos barajados), Watts-Strogatz
(p en 0.01, 0.05, 0.1, 0.2, 1), RGG2, RGG3 (control positivo) y arbol balanceado. Se registran codigos y la tasa con que
cada bandera de corrida se cumple (falso positivo en los nulos). Preregistrado:
  * full (N=800, 5 semillas): RGG3 k12 pasa todos los campos de corrida en el 100%; ninguna otra familia (salvo RGG2,
    solo informe) pasa en ninguna corrida; WS p=0.05 tiene primario F3 en el 100%;
  * smoke (N=100, 1 semilla, ER 6/12, WS 0.05/1): ninguna familia nula (ER, nulos del RGG3, WS, arbol) pasa; RGG3/RGG2 solo informe (N=100 no tiene
    ventana 3D, DESIGN §0).
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.certificate.taxonomy import RUN_FLAG_KEYS, assess_run
from omega.config.seeds import SeedKey, make_rng, seed_key
from omega.config.settings11 import NullModel, Omega11Config
from omega.controls.battery import null_battery, null_seed_key
from omega.controls.erdos_renyi import erdos_renyi_gnp
from omega.controls.random_geometric import balanced_tree, rgg_torus
from omega.controls.small_world import watts_strogatz
from omega.experiments.v11.gate import PassportWriter, derive_cfg, require_prerequisites, runs_root, write_summary
from omega.phases.scan import default_config
from omega.types import FloatArray, RunStatus

NAME = "o03_nulls"
EXPERIMENT_ID = 1103
ENTRYPOINT = "omega.experiments.v11.test_o03_nulls:run"
FIELD_FLAGS = tuple(k for k in RUN_FLAG_KEYS if k not in ("converged", "empty", "dense_or_uniform", "small_world"))
WS_P = (0.01, 0.05, 0.1, 0.2, 1.0)

# Cada generador devuelve [(subfamilia, estado, SeedKey de evidencia, null_model | None)].
Gen = Callable[[SeedKey, int], list[tuple[str, FloatArray, SeedKey, NullModel | None]]]


def _generators(cfg: Omega11Config, mode: Literal["smoke", "full"]) -> tuple[int, int, list[tuple[str, Gen]]]:
    n = 100 if mode == "smoke" else 800
    reps = 1 if mode == "smoke" else 5
    h = 3 if mode == "smoke" else 5
    w_min = cfg.base.graph.w_min
    gens: list[tuple[str, Gen]] = []

    def er(kbar: float) -> Gen:
        def g(key: SeedKey, nn: int) -> list[tuple[str, FloatArray, SeedKey, NullModel | None]]:
            return [(f"er_k{kbar:g}", erdos_renyi_gnp(nn, kbar / (nn - 1), make_rng(key)), key, NullModel.ERDOS_RENYI)]
        return g

    def ws(p: float) -> Gen:
        def g(key: SeedKey, nn: int) -> list[tuple[str, FloatArray, SeedKey, NullModel | None]]:
            return [(f"ws_p{p:g}", watts_strogatz(nn, 12, p, make_rng(key)), key, NullModel.SMALL_WORLD)]
        return g

    def rgg(dim: int, kbar: float) -> Gen:
        def g(key: SeedKey, nn: int) -> list[tuple[str, FloatArray, SeedKey, NullModel | None]]:
            return [(f"rgg{dim}_k{kbar:g}", rgg_torus(nn, dim, kbar, make_rng(key)), key, NullModel.RANDOM_GEOMETRIC)]
        return g

    def derived(key: SeedKey, nn: int) -> list[tuple[str, FloatArray, SeedKey, NullModel | None]]:
        base = rgg_torus(nn, 3, 12.0, make_rng(key))
        nulls = null_battery(base, w_min, cfg.null_models, key)
        return [(f"rgg3_{k.value}", w, null_seed_key(key, k), k) for k, w in nulls.items()]

    def tree(key: SeedKey, nn: int) -> list[tuple[str, FloatArray, SeedKey, NullModel | None]]:
        return [(f"tree3_h{h}", balanced_tree(3, h), key, None)]

    for kb in (6.0, 12.0) if mode == "smoke" else (6.0, 12.0, 30.0):
        gens.append((f"er_k{kb:g}", er(kb)))
    gens.append(("rgg3_k12", rgg(3, 12.0)))
    gens.append(("rgg3_derived_nulls", derived))
    for p in (0.05, 1.0) if mode == "smoke" else WS_P:
        gens.append((f"ws_p{p:g}", ws(p)))
    gens.append(("rgg2_k10", rgg(2, 10.0)))
    gens.append((f"tree3_h{h}", tree))
    return n, reps, gens


def _expected(mode: Literal["smoke", "full"], names: list[str]) -> list[dict[str, Any]]:
    informative = {"rgg3_k12", "rgg2_k10"}
    out: list[dict[str, Any]] = []
    for nm in names:
        if mode == "full" and nm == "rgg3_k12":
            out.append({"family": "rgg3_k12", "metric": "pass_fraction", "target": 1.0})
        elif nm not in informative:
            out.append({"family": nm, "metric": "pass_fraction", "target": 0.0})
    if mode == "full":
        out.append({"family": "ws_p0.05", "metric": "primary_F3_fraction", "target": 1.0})
    return out


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-03 y devuelve `summary.json`. En full exige O-00..O-02 full completos."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    n, reps, gens = _generators(cfg, mode)
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=reps)
    expected = _expected(mode, [g[0] for g in gens])
    header: dict[str, Any] = {
        "experiment": EXPERIMENT_ID, "description": "bateria de nulos y falsos positivos por campo (DESIGN §1.10)",
        "n": n, "replicates": reps, "expected": expected,
        "note": "familias 'rgg3_k12' y 'rgg2_k10' son controles positivos; en smoke solo informe",
    }
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    collected: dict[str, list[dict[str, Any]]] = {}
    for gi, (_, gen) in enumerate(gens):
        for rep in range(reps):
            key = seed_key(cfg.base.seeds, gi, rep)
            for sub, w, ekey, null in gen(key, n):
                c_case = derive_cfg(cfg, n=int(w.shape[0]))
                ev = collect_run_evidence(w, RunStatus.CONVERGED, c_case, evidence_rng(ekey))
                a = assess_run(ev, c_case)
                row: dict[str, Any] = {
                    "codes": [c.value for c in a.codes], "primary": None if a.primary is None else a.primary.value,
                    "passes": a.passes, "flags": dict(a.flags),
                    "d_vol": ev.geometry.d_eff.value, "d_s": ev.geometry.d_s.value, "d_weyl": ev.d_weyl.value,
                    "consensus_dimension": ev.consensus_dimension,
                }
                if rep == 0:
                    row["passport"] = writer.save(c_case, seed=ekey, w=w, termination="STATIC", init_distribution=f"control/{sub}",
                                                  results={k: v for k, v in row.items() if k != "flags"}, null_model=null)
                collected.setdefault(sub, []).append(row)

    results: dict[str, Any] = {}
    for sub, rows in collected.items():
        m = len(rows)
        results[sub] = {
            "n_runs": m,
            "pass_fraction": sum(1 for r in rows if r["passes"]) / m,
            "primary_F3_fraction": sum(1 for r in rows if r["primary"] == "Ω-F3") / m,
            "code_fractions": {k: v / m for k, v in sorted(Counter(c for r in rows for c in r["codes"]).items())},
            "primary_fractions": {k: v / m for k, v in sorted(Counter(str(r["primary"]) for r in rows).items())},
            "field_true_rate": {f: sum(1 for r in rows if r["flags"][f]) / m for f in FIELD_FLAGS},
            "rows": [{k: v for k, v in r.items() if k != "flags"} for r in rows],
        }
    met: dict[str, bool] = {}
    for e in expected:
        fam = e["family"]
        # La familia preregistrada puede ser un grupo (rgg3_derived_nulls): se evalua cada subfamilia.
        subs = [s for s in results if s == fam or s.startswith("rgg3_") and fam == "rgg3_derived_nulls"]
        if fam == "rgg3_derived_nulls":
            subs = [s for s in results if s.startswith("rgg3_") and s != "rgg3_k12"]
        for s in subs:
            met[f"{s}.{e['metric']}"] = bool(results[s][e["metric"]] == e["target"])
    summary = {
        **header, "stage": "final", "results": results, "expectation_results": met,
        "expectations_met": all(met.values()), "passports": writer.labels, "complete": True,
    }
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(100, master_entropy=20240901, replicates=1))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    res = data["results"]
    assert {"er_k6", "er_k12", "rgg3_k12", "rgg2_k10", "ws_p0.05", "rgg3_degree_preserving_rewire",
            "rgg3_configuration_model", "rgg3_shuffled_weights", "tree3_h3"} <= set(res)
    assert all(set(r["field_true_rate"]) == set(FIELD_FLAGS) for r in res.values())
    assert data["expectations_met"] is True  # ningun nulo pasa a N=100
    assert len(data["passports"]) == len(res)


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(800, master_entropy=20240901, replicates=5))
    data = json.loads(run(cfg, runs_root(), mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
