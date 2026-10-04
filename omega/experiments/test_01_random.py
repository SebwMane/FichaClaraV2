"""Experimento 1 (M§23): linea base aleatoria U(0,1) SIN dinamica.

Mide conectividad, clustering, componente gigante, distancias, D_eff y D_s sobre W0 aleatoria
y barre w_min (R2). Resultado esperado (ANALYSIS §2, R6): D_eff = no_window (no D~3 ingenuo);
la red es hiperdensa (E) porque rho(W>w_min) ~ 1-w_min. No es un criterio de exito: se registra.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import pytest

from omega.config.seeds import make_rng, seed_key
from omega.config.settings import FunctionalParams, OmegaConfig
from omega.io.passport import save_run
from omega.network.initialization import random_uniform_weights
from omega.phases.classification import classify_run
from omega.phases.scan import default_config, observe, run_summary, static_run, to_jsonable, with_experiment
from omega.phases.stability import seed_summary
from omega.types import FloatArray, RunResult, RunStatus

EXPERIMENT_ID = 1


def run(cfg: OmegaConfig, out_dir: Path, git_commit: str = "unknown") -> Path:
    """Ejecuta el Experimento 1; guarda pasaportes en out_dir/passports y devuelve summary.json."""
    cfg = with_experiment(cfg, EXPERIMENT_ID)
    out_dir = Path(out_dir)
    marker = FunctionalParams(alpha=0.0)  # sin dinamica: solo marcador de parametros
    results: list[RunResult] = []
    states: list[FloatArray] = []
    for rep in range(cfg.seeds.replicates):
        key = seed_key(cfg.seeds, 0, rep)
        w0 = random_uniform_weights(cfg.init.n, make_rng(key))
        res = static_run(cfg, marker, key, w0)
        save_run(res, cfg, out_dir / "passports", f"exp01_rep{rep:03d}", git_commit)
        results.append(res)
        states.append(w0)

    n = len(results)
    no_window = sum(1 for r in results if r.observables.geometry.d_eff.status == "no_window")
    sweep: dict[str, Any] = {}
    for wm in cfg.graph.w_min_sensitivity:
        c = dataclasses.replace(cfg, graph=dataclasses.replace(cfg.graph, w_min=wm))
        rows = []
        for w in states:
            o = observe(w, c)
            rows.append(
                {
                    "label": classify_run(o, RunStatus.CONVERGED, c.phases).label,
                    "binary_density": o.topology.binary_density,
                    "giant_fraction": o.topology.giant_fraction,
                    "clustering_binary": o.topology.clustering_binary,
                    "d_eff_status": o.geometry.d_eff.status,
                    "d_eff": o.geometry.d_eff.value,
                    "d_s_status": o.geometry.d_s.status,
                    "d_s": o.geometry.d_s.value,
                }
            )
        sweep[f"{wm:g}"] = {
            "no_window_fraction": sum(1 for r in rows if r["d_eff_status"] == "no_window") / n,
            "rows": rows,
        }

    def stat(fn: Any) -> dict[str, float]:
        m, s, k = seed_summary(results, fn)
        return {"mean": m, "std": s, "n": float(k)}

    summary = {
        "experiment": EXPERIMENT_ID,
        "description": "linea base aleatoria U(0,1) sin dinamica (M§23)",
        "n_nodes": cfg.init.n,
        "replicates": n,
        "w_min": cfg.graph.w_min,
        "expected": {"d_eff_status": "no_window", "phase": "E"},
        "d_eff_no_window_fraction": no_window / n,
        "expectation_no_window_met": no_window == n,
        "observables": {
            "mean_strength": stat(lambda r: r.observables.topology.mean_strength),
            "std_strength": stat(lambda r: r.observables.topology.std_strength),
            "binary_density": stat(lambda r: r.observables.topology.binary_density),
            "giant_fraction": stat(lambda r: r.observables.topology.giant_fraction),
            "clustering_weighted": stat(lambda r: r.observables.topology.clustering_weighted),
            "clustering_binary": stat(lambda r: r.observables.topology.clustering_binary),
            "path_length": stat(lambda r: r.observables.geometry.path_length),
            "d_s": stat(lambda r: r.observables.geometry.d_s.value),
        },
        "label_counts": {
            lab: sum(1 for r in results if r.assessment.label.value == lab)
            for lab in sorted({r.assessment.label.value for r in results})
        },
        "w_min_sweep": sweep,
        "runs": [run_summary(r) for r in results],
    }
    path = out_dir / "summary.json"
    path.write_text(json.dumps(to_jsonable(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path


def test_smoke(tmp_path: Path) -> None:
    cfg = default_config(24, master_entropy=101, replicates=2)
    path = run(cfg, tmp_path, git_commit="test")
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["replicates"] == 2 and len(data["runs"]) == 2
    assert data["d_eff_no_window_fraction"] == 1.0 and data["expectation_no_window_met"] is True
    assert set(data["w_min_sweep"]) == {"0.01", "0.05", "0.1", "0.2", "0.5"}
    assert data["label_counts"] == {"E": 2}
    assert len(list((tmp_path / "passports").glob("*.json"))) == 2


@pytest.mark.slow
def test_full_run(tmp_path: Path) -> None:
    cfg = default_config(200, master_entropy=20240901, replicates=10)
    data = json.loads(run(cfg, tmp_path).read_text(encoding="utf-8"))
    assert data["replicates"] == 10
