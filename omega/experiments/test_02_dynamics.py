"""Experimento 2 (M§24): dinamica de S0 desde pesos aleatorios, registrando toda la evolucion.

Por cada punto y replica se guarda el pasaporte (con escalares por paso e instantaneas log2) y,
ademas, los observables geometricos en cada instantanea (D-4: la geometria solo se calcula ahi).
Puntos: el funcional de cfg y dos puntos reducidos de referencia (alpha_hat = 1.5 y 2.5, mismo
gamma_hat) para confirmar la biestabilidad preliminar de R1 (A por debajo, E por encima de ~2).
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from omega.config.convert import raw_to_reduced
from omega.config.seeds import seed_key
from omega.config.settings import DynamicsConfig, FunctionalParams, OmegaConfig
from omega.io.passport import save_run
from omega.network.weights import from_upper_triangle
from omega.phases.scan import (
    default_config,
    observe,
    point_params,
    run_summary,
    simulate,
    to_jsonable,
    with_experiment,
    with_params,
)
from omega.types import RunResult

EXPERIMENT_ID = 2
REFERENCE_ALPHA_HAT = (1.5, 2.5)
MAX_REPLICATES_PER_POINT = 3


def _snapshot_table(r: RunResult, cfg: OmegaConfig) -> list[dict[str, Any]]:
    """Observables en cada instantanea de la trayectoria (pasos {0,1,2,4,...} y final)."""
    rows: list[dict[str, Any]] = []
    n = cfg.init.n
    for step, vec in zip(r.trajectory.snapshot_steps, r.trajectory.snapshots, strict=True):
        o = observe(from_upper_triangle(vec, n), cfg)
        rows.append(
            {
                "step": int(step),
                "tau": float(step) * r.trajectory.dt,
                "mean_weight": o.topology.mean_weight,
                "mean_strength": o.topology.mean_strength,
                "cv_strength": o.topology.cv_strength,
                "binary_density": o.topology.binary_density,
                "giant_fraction": o.topology.giant_fraction,
                "clustering_binary": o.topology.clustering_binary,
                "frac_at_zero": o.topology.frac_at_zero,
                "frac_at_one": o.topology.frac_at_one,
                "d_eff": o.geometry.d_eff.value,
                "d_eff_status": o.geometry.d_eff.status,
                "d_s": o.geometry.d_s.value,
                "d_s_status": o.geometry.d_s.status,
            }
        )
    return rows


def run(cfg: OmegaConfig, out_dir: Path, git_commit: str = "unknown") -> Path:
    """Ejecuta el Experimento 2; devuelve out_dir/summary.json."""
    cfg = with_experiment(cfg, EXPERIMENT_ID)
    out_dir = Path(out_dir)
    n = cfg.init.n
    _, g_hat = raw_to_reduced(cfg.functional, n)
    points: list[tuple[str, FunctionalParams]] = [("config", cfg.functional)]
    points += [(f"ref_a{a:g}", point_params(a, g_hat, n)) for a in REFERENCE_ALPHA_HAT]

    reps = min(cfg.seeds.replicates, MAX_REPLICATES_PER_POINT)
    records: list[dict[str, Any]] = []
    for pi, (name, p) in enumerate(points):
        pcfg = with_params(cfg, p)
        for rep in range(reps):
            key = seed_key(cfg.seeds, pi, rep)
            res = simulate(pcfg, p, key)
            run_id = f"exp02_{name}_rep{rep:03d}"
            save_run(res, pcfg, out_dir / "passports", run_id, git_commit)
            action = np.asarray(res.trajectory.scalars["action"])
            records.append(
                {
                    "point": name,
                    "run_id": run_id,
                    "summary": run_summary(res),
                    "tau": res.trajectory.tau,
                    "dt": res.trajectory.dt,
                    "action_initial": float(action[0]),
                    "action_final": float(action[-1]),
                    "action_nonincreasing": bool(np.all(np.diff(action) <= 1e-9 * max(1.0, abs(action[0])))),
                    "snapshots": _snapshot_table(res, pcfg),
                }
            )
    summary = {
        "experiment": EXPERIMENT_ID,
        "description": "dinamica de S0 con trayectoria completa (M§24)",
        "n_nodes": n,
        "points": [{"name": nm, "params": dataclasses.asdict(p)} for nm, p in points],
        "replicates_per_point": reps,
        "records": records,
        "label_by_point": {
            nm: [r["summary"]["label"] for r in records if r["point"] == nm] for nm, _ in points
        },
    }
    path = out_dir / "summary.json"
    path.write_text(json.dumps(to_jsonable(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path


def test_smoke(tmp_path: Path) -> None:
    cfg = default_config(16, master_entropy=102, replicates=2)
    cfg = dataclasses.replace(
        cfg, dynamics=DynamicsConfig(max_steps=20_000), functional=point_params(0.5, 1.0, 16)
    )
    data = json.loads(run(cfg, tmp_path, git_commit="test").read_text(encoding="utf-8"))
    assert data["replicates_per_point"] == 2
    labels = data["label_by_point"]
    assert labels["config"] == ["A", "A"]  # alpha_hat=0.5 < 2
    assert labels["ref_a1.5"] == ["A", "A"] and labels["ref_a2.5"] == ["E", "E"]
    rec = data["records"][0]
    assert rec["snapshots"][0]["step"] == 0 and rec["snapshots"][-1]["step"] == rec["summary"]["steps"]
    assert rec["action_nonincreasing"] is True
    assert len(list((tmp_path / "passports").glob("*.npz"))) == 6


@pytest.mark.slow
def test_full_run(tmp_path: Path) -> None:
    cfg = default_config(200, master_entropy=20240901, replicates=10)
    cfg = dataclasses.replace(cfg, functional=point_params(1.0, 1.0, 200))
    data = json.loads(run(cfg, tmp_path).read_text(encoding="utf-8"))
    assert data["replicates_per_point"] == 3
