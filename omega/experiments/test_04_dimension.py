"""Experimento 4 (M§9-M§11, R2, R3, R6): dimension emergente con barrido de w_min.

Para cada estado (linea base aleatoria, estado final de la dinamica de cfg.functional y un toro
3D periodico de control con W=1) y cada w_min de cfg.graph.w_min_sensitivity se informan D_eff
(cascara, bola y saltos), D_s, sus estados/ventanas y la curva D_eff(r) (escalas, perfil y
pendientes locales). Una etiqueta de fase solo se acepta si es estable en >= 4/5 valores de w_min.
"""

from __future__ import annotations

import dataclasses
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from omega.config.seeds import make_rng, seed_key
from omega.config.settings import DynamicsConfig, FunctionalParams, OmegaConfig
from omega.experiments.reference_graphs import periodic_lattice
from omega.io.passport import save_run
from omega.network.initialization import random_uniform_weights
from omega.phases.classification import classify_run
from omega.phases.scan import (
    default_config,
    observe,
    point_params,
    simulate,
    static_run,
    to_jsonable,
    with_experiment,
    with_params,
)
from omega.types import DimensionEstimate, FloatArray, RunStatus

EXPERIMENT_ID = 4
MAX_STATES = 3
CONTROL_SHAPE = (7, 7, 7)


def _dim(e: DimensionEstimate, curve: bool = False) -> dict[str, Any]:
    out: dict[str, Any] = {
        "value": e.value, "stderr": e.stderr, "status": e.status,
        "window": None if e.window is None else list(e.window), "plateau": e.plateau, "method": e.method,
    }
    if curve:
        out.update({"scales": e.scales, "profile": e.profile, "local_slopes": e.local_slopes})
    return out


def _sweep(w: FloatArray, cfg: OmegaConfig, status: RunStatus = RunStatus.CONVERGED) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    labels: list[str] = []
    for wm in cfg.graph.w_min_sensitivity:
        c = dataclasses.replace(cfg, graph=dataclasses.replace(cfg.graph, w_min=wm))
        o = observe(w, c)
        label = classify_run(o, status, c.phases).label.value
        labels.append(label)
        rows[f"{wm:g}"] = {
            "label": label,
            "giant_fraction": o.topology.giant_fraction,
            "binary_density": o.topology.binary_density,
            "path_length": o.geometry.path_length,
            "d_eff_shell": _dim(o.geometry.d_eff, curve=True),  # curva D_eff(r)
            "d_eff_ball": _dim(o.geometry.d_eff_ball),
            "d_eff_hops": _dim(o.geometry.d_eff_hops),
            "d_s": _dim(o.geometry.d_s),
        }
    modal, count = Counter(labels).most_common(1)[0]
    return {
        "by_w_min": rows,
        "modal_label": modal,
        "label_agreement": count / len(labels),
        "label_stable_4_of_5": count / len(labels) >= 0.8,
    }


def run(cfg: OmegaConfig, out_dir: Path, git_commit: str = "unknown") -> Path:
    """Ejecuta el Experimento 4; devuelve out_dir/summary.json."""
    cfg = with_experiment(cfg, EXPERIMENT_ID)
    out_dir = Path(out_dir)
    n = cfg.init.n
    reps = min(cfg.seeds.replicates, MAX_STATES)
    states: dict[str, Any] = {}
    marker = FunctionalParams(alpha=0.0)
    for rep in range(reps):
        key = seed_key(cfg.seeds, 0, rep)
        w0 = random_uniform_weights(n, make_rng(key))
        base = static_run(cfg, marker, key, w0)
        save_run(base, cfg, out_dir / "passports", f"exp04_baseline_rep{rep:03d}", git_commit)
        states[f"baseline_rep{rep}"] = _sweep(w0, cfg)

        ekey = seed_key(cfg.seeds, 1, rep)
        ecfg = with_params(cfg, cfg.functional)
        evolved = simulate(ecfg, cfg.functional, ekey)
        save_run(evolved, ecfg, out_dir / "passports", f"exp04_evolved_rep{rep:03d}", git_commit)
        states[f"evolved_rep{rep}"] = _sweep(evolved.trajectory.w_final, cfg, evolved.trajectory.status)
    states[f"control_torus_{'x'.join(map(str, CONTROL_SHAPE))}"] = _sweep(periodic_lattice(CONTROL_SHAPE), cfg)

    summary = {
        "experiment": EXPERIMENT_ID,
        "description": "dimension emergente con barrido de w_min e informe D_eff(r) (M§9-M§11)",
        "n_nodes": n,
        "w_min_values": list(cfg.graph.w_min_sensitivity),
        "functional": dataclasses.asdict(cfg.functional),
        "states": states,
    }
    path = out_dir / "summary.json"
    path.write_text(json.dumps(to_jsonable(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path


def test_smoke(tmp_path: Path) -> None:
    cfg = default_config(20, master_entropy=104, replicates=1)
    cfg = dataclasses.replace(
        cfg, dynamics=DynamicsConfig(max_steps=20_000), functional=point_params(3.0, 0.0, 20)
    )
    data = json.loads(run(cfg, tmp_path, "test").read_text(encoding="utf-8"))
    assert set(data["states"]) == {"baseline_rep0", "evolved_rep0", "control_torus_7x7x7"}
    base = data["states"]["baseline_rep0"]
    assert set(base["by_w_min"]) == {"0.01", "0.05", "0.1", "0.2", "0.5"}
    assert all(v["d_eff_shell"]["status"] == "no_window" for v in base["by_w_min"].values())
    ctrl = data["states"]["control_torus_7x7x7"]["by_w_min"]["0.1"]
    assert ctrl["d_eff_shell"]["status"] == "ok" and 2.4 < ctrl["d_eff_shell"]["value"] < 3.3
    assert len(ctrl["d_eff_shell"]["scales"]) == len(ctrl["d_eff_shell"]["local_slopes"]) > 0
    assert data["states"]["evolved_rep0"]["modal_label"] == "E"


@pytest.mark.slow
def test_full_run(tmp_path: Path) -> None:
    cfg = default_config(200, master_entropy=20240901, replicates=10)
    cfg = dataclasses.replace(cfg, functional=point_params(1.0, 1.0, 200))
    data = json.loads(run(cfg, tmp_path).read_text(encoding="utf-8"))
    assert len(data["states"]) == 7
