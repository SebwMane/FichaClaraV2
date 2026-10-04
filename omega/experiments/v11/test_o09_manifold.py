"""Omega-9 (Ω-1.1): homogeneidad, isotropia, localidad, anillos y MANIFOLD_PROXY en estados Ω y controles (DESIGN §1.4, §1.5, §4).

Catalogo de estados: `state_families` de O-07. Por familia se informan homogeneidad (cv de D_i), isotropia (MDS local, solo
diagnostico de una metrica ya medida), localidad (fraccion de aristas con desvio <= 3 saltos), conectividad de anillos y la
bandera `manifold_proxy_ok` (anillos ∧ curvatura ∧ b1_density <= 0.03; en el codigo es MANIFOLD_PROXY, nunca "manifold").
O-09 no tiene orden propio en el STEP_ORDER (§4): su modo full exige lo mismo que O-08.

Preregistrado (full, N~800; valores medidos en DESIGN §1.4-§1.5, sin ajustar umbrales, M§44):
  * T³ y RGG3 k12: homogeneity_ok, isotropy_ok, locality, annulus_ok y manifold_proxy_ok True;
  * Watts-Strogatz p=0.05: locality False (fraccion medida 0.949 < 0.97);
  * anillo y arbol 3-ario: manifold_proxy_ok False (conectividad de anillos 0.00);
  * estados Ω: solo informe.
En smoke (N=60) las expectativas se informan pero no se afirman (la ventana de escala exige N grande, R1).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

import pytest

from omega.certificate.evidence import collect_run_evidence, evidence_rng
from omega.certificate.taxonomy import assess_run
from omega.config.settings11 import Omega11Config
from omega.contracts import RunEvidence
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
from omega.phases.finite_size import evidence_cfg
from omega.phases.scan import default_config
from omega.types import RunStatus

NAME = "o09_manifold"
EXPERIMENT_ID = 1109
ENTRYPOINT = "omega.experiments.v11.test_o09_manifold:run"
PREREQ_STEP = "o08_topology"  # sin orden propio: lo mismo que O-08
MANIFOLD_FLAGS = ("homogeneity_ok", "isotropy_ok", "locality", "manifold_proxy_ok")


def manifold_row(ev: RunEvidence, cfg: Omega11Config) -> dict[str, Any]:
    flags = assess_run(ev, cfg).flags
    return {
        "homogeneity": {"mu": ev.homogeneity.mu, "sigma": ev.homogeneity.sigma, "cv": ev.homogeneity.cv,
                        "n_valid": ev.homogeneity.n_valid, "ok": ev.homogeneity.ok},
        "isotropy": {"median_ratio": ev.isotropy.median_ratio, "p10_ratio": ev.isotropy.p10_ratio, "radius": ev.isotropy.radius,
                     "k": ev.isotropy.k, "n_sources": ev.isotropy.n_sources, "ok": ev.isotropy.ok},
        "locality": {"detour_fraction": ev.locality.detour_fraction, "ok": ev.locality.ok},
        "annulus": {"fractions": {str(k): v for k, v in ev.annulus.fractions.items()}, "ok": ev.annulus.ok},
        "curvature_ok": ev.curvature.ok, "b1_density_w_min": ev.topo.at_w_min.b1_density, "clustering_ratio": ev.clustering_ratio,
        "flags": {k: bool(flags[k]) for k in (*MANIFOLD_FLAGS, "small_world")},
        "MANIFOLD_PROXY": bool(flags["manifold_proxy_ok"]),
    }


def run(cfg: Omega11Config, out_root: Path, *, mode: Literal["smoke", "full"]) -> Path:
    """Ejecuta O-09 y devuelve `summary.json`. En full exige lo mismo que O-08 (sin orden propio)."""
    if mode == "full":
        require_prerequisites(out_root, PREREQ_STEP)
    n = 60 if mode == "smoke" else 800
    cfg = derive_cfg(cfg, n=n, experiment_id=EXPERIMENT_ID, replicates=1)
    expected: dict[str, Any] = {
        "all_true": {"torus3d": list(MANIFOLD_FLAGS), "rgg3_k12": list(MANIFOLD_FLAGS)},
        "ws_p0.05": {"locality": False},
        "manifold_proxy_false": ["ring", "tree3"],
        "omega_states": "solo informe", "asserted": mode == "full",
    }
    header: dict[str, Any] = {"experiment": EXPERIMENT_ID, "description": "homogeneidad, isotropia, localidad, anillos y MANIFOLD_PROXY",
                              "n": n, "expected": expected, "prerequisites": PREREQ_STEP if mode == "full" else "sin compuerta (smoke)"}
    write_summary(out_root, NAME, {**header, "stage": "preregistered", "results": None, "complete": False}, mode=mode)

    writer = PassportWriter(out_root, NAME, ENTRYPOINT)
    results: dict[str, Any] = {}
    for fam in state_families(cfg, n, mode):
        w, c_case, key = fam["w"], fam["cfg"], fam["key"]
        c_ev = evidence_cfg(c_case, w)
        st = RunStatus(fam["termination"]) if fam["termination"] in {s.value for s in RunStatus} else RunStatus.CONVERGED
        ev = collect_run_evidence(w, st, c_ev, evidence_rng(key))
        row = assessment_row(ev, c_ev)
        results[fam["name"]] = {"kind": fam["kind"], "n": int(w.shape[0]), **manifold_row(ev, c_ev), "assessment": row,
                                "passport": writer.save(c_case, seed=key, w=w, termination=fam["termination"],
                                                        init_distribution=fam["init"], results=row, null_model=fam["null_model"])}
    met: dict[str, bool] = {}
    for name, flags in expected["all_true"].items():
        met[f"{name}.all_true"] = all(check_expectation("eq", results[name]["flags"][f], True) for f in flags)
    met["ws_p0.05.locality_false"] = check_expectation("eq", results["ws_p0.05"]["flags"]["locality"], False)
    for prefix in expected["manifold_proxy_false"]:
        name = next(k for k in results if k.startswith(prefix))
        met[f"{name}.manifold_proxy_false"] = check_expectation("eq", results[name]["MANIFOLD_PROXY"], False)
    summary = {**header, "stage": "final", "results": results, "expectation_results": met,
               "expectations_met": all(met.values()), "passports": writer.labels, "complete": True}
    return write_summary(out_root, NAME, summary, mode=mode)


def test_smoke(tmp_path: Path) -> None:
    cfg = Omega11Config(base=default_config(60, master_entropy=20240901, replicates=1))
    data = json.loads(run(cfg, tmp_path, mode="smoke").read_text(encoding="utf-8"))
    assert data["step"] == NAME and data["mode"] == "smoke" and data["complete"] is True
    res = data["results"]
    assert {"ring", "torus3d", "rgg3_k12", "ws_p0.05", "omega_u01_static"} <= set(res)
    for r in res.values():
        assert set(r["flags"]) >= set(MANIFOLD_FLAGS) and isinstance(r["MANIFOLD_PROXY"], bool)
    cv = res["torus3d"]["homogeneity"]["cv"]
    assert cv is None or cv == pytest.approx(0.0, abs=1e-9)  # toro: nodos equivalentes (None: sin ventana con N=64)
    assert len(data["passports"]) == len(res) and isinstance(data["expectations_met"], bool)


@pytest.mark.slow
def test_full_run() -> None:
    cfg = Omega11Config(base=default_config(800, master_entropy=20240901, replicates=1))
    out = runs_root()
    require_prerequisites(out, PREREQ_STEP)  # sin orden propio: lo mismo que O-08
    data = json.loads(run(cfg, out, mode="full").read_text(encoding="utf-8"))
    assert data["complete"] is True
