"""Experimento 3 (M§25-M§26, D-8, D-15): mapa de fases en (alpha_hat, gamma_hat).

Se barren DOS mallas con varias semillas por punto: la reducida de D-8 (alpha_hat in {0,0.25,..,4}
mas [1.5,2.5] a paso 0.02; gamma_hat in {0,0.1,0.3,1,3,10,30,100}) y la cruda de M§25
(alpha/beta, gamma/beta = 0, 0.1, ..., 1.0; se corta en 1.0 por coste, decision documentada).
Se registra EXPLICITAMENTE el resultado negativo si solo aparecen A/E (R1, M§20, M§44): no se
reajustan umbrales ni parametros. Se guarda el pasaporte de la replica 0 de cada punto.
"""

from __future__ import annotations

import dataclasses
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from omega.config.convert import raw_to_reduced
from omega.config.seeds import seed_key
from omega.config.settings import DynamicsConfig, OmegaConfig, ScanConfig
from omega.io.passport import save_run
from omega.phases.scan import (
    default_config,
    parameter_grid,
    phase_probabilities,
    point_params,
    run_summary,
    simulate,
    to_jsonable,
    with_experiment,
    with_params,
)
from omega.types import RunResult

EXPERIMENT_ID = 3
TRIVIAL_PHASES = frozenset({"A", "E"})


def default_reduced_scan() -> ScanConfig:
    """Malla reducida preregistrada (D-8)."""
    coarse = np.arange(0.0, 4.0 + 1e-9, 0.25)
    fine = np.arange(1.5, 2.5 + 1e-9, 0.02)
    alphas = sorted({round(float(a), 10) for a in np.concatenate([coarse, fine])})
    return ScanConfig(
        alpha_hat=tuple(alphas),
        gamma_hat=(0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0),
        scaling="reduced",
    )


def default_raw_scan() -> ScanConfig:
    """Malla cruda de M§25: alpha/beta y gamma/beta = 0, 0.1, ..., 1.0 (beta=1)."""
    vals = tuple(round(0.1 * i, 10) for i in range(11))
    return ScanConfig(alpha_hat=vals, gamma_hat=vals, scaling="raw")


def _scan_grid(
    cfg: OmegaConfig, scan: ScanConfig, out_dir: Path, tag: str, git_commit: str
) -> dict[str, Any]:
    gcfg = dataclasses.replace(cfg, scan=scan)
    results: list[RunResult] = []
    for idx, a, g in parameter_grid(scan):
        p = point_params(a, g, cfg.init.n, scan)
        pcfg = with_params(gcfg, p)
        for rep in range(cfg.seeds.replicates):
            res = simulate(pcfg, p, seed_key(cfg.seeds, idx, rep))
            results.append(res)
            if rep == 0:
                save_run(res, pcfg, out_dir / f"passports_{tag}", f"exp03_{tag}_pt{idx:04d}", git_commit)
    labels = Counter(r.assessment.label.value for r in results)
    probs = phase_probabilities(results)
    nontrivial = sorted(set(labels) - TRIVIAL_PHASES)
    only_ae = set(labels) <= TRIVIAL_PHASES
    return {
        "scaling": scan.scaling,
        "alpha_grid": list(scan.alpha_hat),
        "gamma_grid": list(scan.gamma_hat),
        "n_points": len(parameter_grid(scan)),
        "n_runs": len(results),
        "label_counts": dict(sorted(labels.items())),
        "labels_other_than_A_E": nontrivial,
        "negative_result_only_A_E": only_ae,
        "negative_result_note": (
            "RESULTADO NEGATIVO (R1, M§20): solo aparecen las fases triviales A/E; la funcional S0 "
            "con mu=0 no produce estructura intermedia en esta malla. Registrado sin reajustes (M§44)."
            if only_ae
            else f"Aparecen fases distintas de A/E: {nontrivial}; ver mapa por punto."
        ),
        "transition_alpha_hat_by_gamma": _transition(probs, scan),
        "points": [
            {
                "alpha_hat": k[0],
                "gamma_hat": k[1],
                "n": v["n"],
                "counts": {lab: c for lab, c in v["counts"].items() if c > 0},
                "ci": {lab: v["ci"][lab] for lab, c in v["counts"].items() if c > 0},
            }
            for k, v in probs.items()
        ],
        "runs": [run_summary(r) for r in results],
    }


def _transition(probs: dict[tuple[float, float], dict[str, Any]], scan: ScanConfig) -> dict[str, Any]:
    """Para cada gamma_hat: menor alpha_hat (reducido) donde P(E)>=0.5 (None si no hay)."""
    out: dict[str, Any] = {}
    for g in scan.gamma_hat:
        first = None
        for a in sorted(scan.alpha_hat):
            v = probs.get((round(a, 9), round(g, 9)))
            if v is not None and v["probabilities"]["E"] >= 0.5:
                first = a
                break
        out[f"{g:g}"] = first
    return out


def run(
    cfg: OmegaConfig,
    out_dir: Path,
    git_commit: str = "unknown",
    *,
    reduced_scan: ScanConfig | None = None,
    raw_scan: ScanConfig | None = None,
) -> Path:
    """Ejecuta el Experimento 3 (malla reducida y cruda); devuelve out_dir/summary.json."""
    cfg = with_experiment(cfg, EXPERIMENT_ID)
    out_dir = Path(out_dir)
    if reduced_scan is None:
        reduced_scan = cfg.scan if cfg.scan is not None and cfg.scan.scaling == "reduced" else default_reduced_scan()
    if raw_scan is None:
        raw_scan = cfg.scan if cfg.scan is not None and cfg.scan.scaling == "raw" else default_raw_scan()
    if reduced_scan.scaling != "reduced" or raw_scan.scaling != "raw":
        raise ValueError("reduced_scan debe ser 'reduced' y raw_scan 'raw'")
    grids = {
        "reduced": _scan_grid(cfg, reduced_scan, out_dir, "reduced", git_commit),
        "raw": _scan_grid(cfg, raw_scan, out_dir, "raw", git_commit),
    }
    raw_hats = sorted(
        {
            round(raw_to_reduced(point_params(a, g, cfg.init.n, raw_scan), cfg.init.n)[0], 6)
            for a in raw_scan.alpha_hat
            for g in raw_scan.gamma_hat
        }
    )
    summary = {
        "experiment": EXPERIMENT_ID,
        "description": "mapa de fases (alpha_hat, gamma_hat); mallas reducida y cruda (M§25)",
        "n_nodes": cfg.init.n,
        "replicates": cfg.seeds.replicates,
        "raw_grid_equivalent_alpha_hat": raw_hats,
        "negative_result_only_A_E": all(g["negative_result_only_A_E"] for g in grids.values()),
        "grids": grids,
    }
    path = out_dir / "summary.json"
    path.write_text(json.dumps(to_jsonable(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path


def test_default_grids() -> None:
    r = default_reduced_scan()
    assert r.scaling == "reduced" and r.alpha_hat[0] == 0.0 and r.alpha_hat[-1] == 4.0
    assert len(r.alpha_hat) == 65
    assert 2.0 in r.alpha_hat and 1.98 in r.alpha_hat and 2.02 in r.alpha_hat
    assert len(set(r.alpha_hat)) == len(r.alpha_hat) and list(r.alpha_hat) == sorted(r.alpha_hat)
    assert r.gamma_hat == (0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
    w = default_raw_scan()
    assert w.scaling == "raw" and len(w.alpha_hat) == 11 and w.alpha_hat[1] == 0.1


def test_smoke(tmp_path: Path) -> None:
    cfg = default_config(16, master_entropy=103, replicates=2)
    cfg = dataclasses.replace(cfg, dynamics=DynamicsConfig(max_steps=20_000))
    red = ScanConfig(alpha_hat=(0.0, 1.0, 3.0), gamma_hat=(0.0, 1.0))
    raw = ScanConfig(alpha_hat=(0.0, 0.5), gamma_hat=(0.0, 0.5), scaling="raw")
    path = run(cfg, tmp_path, "test", reduced_scan=red, raw_scan=raw)
    data = json.loads(path.read_text(encoding="utf-8"))
    g = data["grids"]
    assert g["reduced"]["n_runs"] == 12 and g["raw"]["n_runs"] == 8
    assert g["reduced"]["label_counts"] == {"A": 8, "E": 4}
    assert g["raw"]["label_counts"] == {"A": 4, "E": 4}  # alpha/beta=0.5 -> alpha_hat=3.5: E
    assert data["negative_result_only_A_E"] is True
    assert "RESULTADO NEGATIVO" in g["reduced"]["negative_result_note"]
    assert g["reduced"]["transition_alpha_hat_by_gamma"] == {"0": 3.0, "1": 3.0}
    assert len(list((tmp_path / "passports_reduced").glob("*.json"))) == 6
    assert len(list((tmp_path / "passports_raw").glob("*.json"))) == 4
    with pytest.raises(ValueError):
        run(cfg, tmp_path, reduced_scan=raw, raw_scan=raw)


@pytest.mark.slow
def test_full_run(tmp_path: Path) -> None:
    cfg = default_config(200, master_entropy=20240901, replicates=10)
    data = json.loads(run(cfg, tmp_path).read_text(encoding="utf-8"))
    assert data["grids"]["reduced"]["n_points"] == 65 * 8
