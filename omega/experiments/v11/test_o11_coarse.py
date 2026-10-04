"""Omega-11 (Ω-1.1, paso 10): coarse-graining preregistrado y consistencia de clase (DESIGN §1.7, §4, riesgos R1 y R8).

Regla `heavy_edge_matching/max/v1` (congelada antes de este experimento), 3 replicas del rng, niveles hasta n' < 50 o 4 niveles.
Por familia se informan las clases por nivel, los codigos nuevos y `CoarseConsistency` (`class_consistency` de
`omega.coarse_graining.analysis`): consistente := >= 2 niveles medibles, la misma `dimension_class` en todos y ningun codigo de
fallo nuevo respecto al nivel 0; los niveles sin clase (`no_window`) cuentan como indeterminados, no como fallo.

Preregistrado (sin ajustar umbrales, M§44):
  * anillo (smoke N=200; full N=800 y XL): consistente, con clase 1 en todos los niveles medibles;
  * toro 2D (full N=784 y XL): consistente, con clase 2; en smoke (16², N=256) se informa pero no se afirma: sus niveles gruesos
    (n' ~ 140, 76) no tienen ventana de escala (R1: la ventana exige N grande), de modo que solo el nivel 0 es medible;
  * T³, RGG2, RGG3 y estados Ω: solo informe (R1: la coherencia en 3D no es certificable con N <= 800).
XL (full): anillo y T² en `finite_size.xl_sizes` y RGG3 en el primer tamano XL. O-11 es el paso 11 del STEP_ORDER (el ultimo).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Literal

import pytest

from omega.certificate.evidence import evidence_rng
from omega.certificate.taxonomy import assess_run
from omega.coarse_graining.analysis import class_consistency, hierarchy_evidence
from omega.config.seeds import make_rng, seed_key
from omega.config.settings11 import NullModel, Omega11Config
from omega.controls.random_geometric import rgg_torus
from omega.experiments.reference_graphs import periodic_lattice
from omega.experiments.v11.gate import (
    PassportWriter,
    assessment_row,
    check_expectation,
    derive_cfg,
    require_prerequisites,
    runs_root,
    write_summary,
)
from omega.experiments.v11.test_o07_dimensions import state_families
from omega.phases.finite_size import evidence_cfg, size_config
from omega.phases.scan import default_config
from omega.types import FloatArray, RunStatus

NAME = "o11_coarse"
EXPERIMENT_ID = 1111
ENTRYPOINT = "omega.experiments.v11.test_o11_coarse:run"


def _cases(cfg: Omega11Config, mode: Literal["smoke", "full"]) -> list[dict[str, Any]]:
    """[{name, w, claim (clase esperada | None), null_model, init}]."""
    cases: list[dict[str, Any]] = []

    def add(name: str, w: FloatArray, claim: int | None, init: str, null: NullModel | None = None) -> None:
        cases.append({"name": name, "w": w, "claim": claim, "init": init, "null_model": null})

    xl = tuple(cfg.finite_size.xl_sizes) if mode == "full" else ()
    if mode == "smoke":
        add("ring_200", periodic_lattice((200,)), 1, "control/ring")
        add("torus2d_16x16", periodic_lattice((16, 16)), None, "control/torus2d")
        return cases
    add("ring_800", periodic_lattice((800,)), 1, "control/ring")
    add("torus2d_28x28", periodic_lattice((28, 28)), 2, "control/torus2d")
    add("torus3d_9x9x9", periodic_lattice((9, 9, 9)), None, "control/torus3d")
    key = seed_key(cfg.base.seeds, 90, 0)
    add("rgg2_k10_800", rgg_torus(800, 2, 10.0, make_rng(key)), None, "control/rgg2_k10", NullModel.RANDOM_GEOMETRIC)
    key = seed_key(cfg.base.seeds, 91, 0)
    add("rgg3_k12_800", rgg_torus(800, 3, 12.0, make_rng(key)), None, "control/rgg3_k12", NullModel.RANDOM_GEOMETRIC)
    for n in xl:
        add(f"ring_{n}", periodic_lattice((n,)), 1, "control/ring")
        s = round(math.sqrt(n))
        add(f"torus2d_{s}x{s}", periodic_lattice((s, s)), 2, "control/torus2d")
    if xl:
        key = seed_key(cfg.base.seeds, 92, 0)
        add(f"rgg3_k12_{xl[0]}", rgg_torus(xl[0], 3, 12.0, make_rng(key)), None, "control/rgg3_k12", NullModel.RANDOM_GEOMETRIC)
    return cases


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-11 y devuelve `summary.json`. En full exige los pasos 1..10 del STEP_ORDER (el ultimo)."""
    if mode == "full":
        require_prerequisites(out_root, NAME)
    cfg = derive_cfg(cfg, experiment_id=EXPERIMENT_ID, replicates=1)
    cases = _cases(cfg, mode)
    expected: dict[str, Any] = {
        "consistent_with_class": {c["name"]: c["claim"] for c in cases if c["claim"] is not None},
        "informe_only": [c["name"] for c in cases if c["claim"] is None], "rule": "heavy_edge_matching/max/v1",
        "replicates": cfg.coarse.replicates,
    }
    header: dict[str, Any] = {"experiment": EXPERIMENT_ID, "description": "coarse-graining y consistencia de clase", "expected": expected}
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: dict[str, Any] = {}
    all_cases = [(c["name"], c["w"], c["claim"], c["init"], c["null_model"]) for c in cases]
    if mode == "full":  # estados Ω (N_Ω <= 200) tras los controles; solo informe
        for fam in state_families(cfg, 800, "full"):
            if fam["kind"] == "omega":
                all_cases.append((fam["name"], fam["w"], None, fam["init"], None))
    for ci, (name, w, claim, init, null) in enumerate(all_cases):
        key = seed_key(cfg.base.seeds, 300 + ci, 0)
        c_case = size_config(cfg, int(w.shape[0]))
        c_ev = evidence_cfg(c_case, w)
        levels = hierarchy_evidence(w, c_ev, evidence_rng(key))
        cc = class_consistency(levels, c_ev)
        row = assessment_row(levels[0], c_ev)
        results[name] = {
            "n": int(w.shape[0]), "claim_class": claim, "level_n": [ev.n for ev in levels],
            "classes": list(cc.classes), "codes": [None if c is None else c.value for c in cc.codes],
            "level_codes": [[c.value for c in assess_run(ev, c_ev).codes] for ev in levels],
            "consensus_dimension": [ev.consensus_dimension for ev in levels],
            "measured_levels": cc.measured_levels, "undetermined_levels": cc.undetermined_levels, "consistent": cc.consistent,
            "claim_met": None if claim is None else bool(
                cc.consistent and all(check_expectation("eq", c, claim) for c in cc.classes if c is not None)),
            "assessment_level0": row,
            "passport": writer.save(c_case, seed=key, w=w, termination="STATIC", init_distribution=init, results=row, null_model=null),
        }
    judged = {k: v["claim_met"] for k, v in results.items() if v["claim_met"] is not None}
    asserted = {k: v for k, v in judged.items() if mode == "full" or k.startswith("ring")}
    summary = {**header, "stage": "final", "results": results, "expectation_results": judged,
               "expectations_met": all(asserted.values()), "passports": writer.labels, "complete": True}
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(200, master_entropy=20240901, replicates=1))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    res = data["results"]
    assert set(res) == {"ring_200", "torus2d_16x16"}
    ring = res["ring_200"]
    assert ring["consistent"] is True and ring["claim_met"] is True and ring["level_n"][0] == 200
    assert {c for c in ring["classes"] if c is not None} == {1} and ring["measured_levels"] >= 3
    assert res["torus2d_16x16"]["claim_met"] is None and res["torus2d_16x16"]["classes"][0] == 2
    assert data["expectations_met"] is True and len(data["passports"]) == 2


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(800, master_entropy=20240901, replicates=1))
    out = runs_root()
    require_prerequisites(out, NAME)  # compuerta explicita: paso 11 del STEP_ORDER
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
