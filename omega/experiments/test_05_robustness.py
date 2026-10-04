"""Experimento 5 (M§26, M§29-M§33, ANALYSIS §4.1): robustez del punto base cfg.functional.

Secciones: semillas, tamanos {100,200,300}, +-10% en (alpha_hat, gamma_hat), ablacion
(-aT / -aT+bS_dens / -aT+gS_deg / completo), invariancia por permutacion (<= 1e-9), nulos
(pesos barajados y recableado con grados preservados), sensibilidad a w_min, histeresis por
continuacion (ascendente y descendente) y confirmacion de F (confirm_geometric). Los umbrales son
los preregistrados; ningun ajuste tras ver resultados (M§44).

Bloques de semilla (point_index de seed_key): 0 base; 1.. tamanos; 100.. vecinos; 200.. ablacion;
300 permutacion; 301 nulos; 400.. histeresis.
"""

from __future__ import annotations

import dataclasses
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from omega.config.convert import raw_to_reduced
from omega.config.seeds import make_rng, seed_key
from omega.config.settings import DynamicsConfig, FunctionalParams, OmegaConfig
from omega.io.passport import save_run
from omega.network.initialization import random_uniform_weights
from omega.network.topology import adjacency
from omega.phases.classification import classify_run
from omega.phases.scan import (
    default_config,
    observe,
    phase_probabilities,
    point_params,
    run_summary,
    simulate,
    to_jsonable,
    with_experiment,
    with_params,
)
from omega.phases.stability import (
    confirm_geometric,
    degree_preserving_null,
    hysteresis_gap,
    hysteresis_sweep,
    permutation_invariance,
    random_permutation,
    seed_summary,
    shuffled_weight_null,
    wmin_candidate_flags,
)
from omega.types import FloatArray, RunResult, RunStatus

EXPERIMENT_ID = 5
BETA_OFF = 1e-9  # beta debe ser > 0: "sin S_dens" se aproxima con beta ~ 0 (el clip acota el problema)
DEFAULT_SIZES = (100, 200, 300)
DEFAULT_ALPHA_PATH = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)


def _runs(cfg: OmegaConfig, p: FunctionalParams, point: int, reps: int) -> list[RunResult]:
    pcfg = with_params(cfg, p)
    return [simulate(pcfg, p, seed_key(cfg.seeds, point, rep)) for rep in range(reps)]


def _deff(r: RunResult) -> float:
    return float(r.observables.geometry.d_eff.value)


def _group_summary(rs: Sequence[RunResult]) -> dict[str, Any]:
    m, s, k = seed_summary(rs, _deff)
    probs = phase_probabilities(rs)
    return {
        "d_eff_mean": m,
        "d_eff_std": s,
        "d_eff_n": k,
        "mean_weight": seed_summary(rs, lambda r: r.observables.topology.mean_weight)[0],
        "phase_probabilities": {f"{a:g},{g:g}": v for (a, g), v in probs.items()},
        "labels": [r.assessment.label.value for r in rs],
        "f_candidate_fraction": sum(1 for r in rs if r.assessment.flags["F_candidate"]) / len(rs),
    }


def _loop(
    w_first: FloatArray, path: Sequence[float], g_hat: float, cfg: OmegaConfig
) -> dict[str, Any]:
    """Ida por `path` desde w_first y vuelta por continuacion; brecha de <W> y de G."""
    fwd = hysteresis_sweep(w_first, path, g_hat, cfg)
    bwd = hysteresis_sweep(fwd[-1].trajectory.w_final, list(path)[::-1], g_hat, cfg)
    gap_w = hysteresis_gap(fwd, bwd, lambda r: r.observables.topology.mean_weight)
    gap_g = hysteresis_gap(fwd, bwd, lambda r: r.observables.topology.giant_fraction)
    return {
        "path": list(path),
        "forward_labels": [r.assessment.label.value for r in fwd],
        "backward_labels_traversal_order": [r.assessment.label.value for r in bwd],
        "gap_mean_weight": gap_w,
        "gap_giant_fraction": gap_g,
        "max_abs_gap_mean_weight": float(np.max(np.abs(gap_w))),
        "hysteresis_detected": bool(np.max(np.abs(gap_w)) > 1e-6),
        "forward_converged": all(r.trajectory.status is RunStatus.CONVERGED for r in fwd),
        "backward_converged": all(r.trajectory.status is RunStatus.CONVERGED for r in bwd),
    }


def run(
    cfg: OmegaConfig,
    out_dir: Path,
    git_commit: str = "unknown",
    *,
    sizes: Sequence[int] = DEFAULT_SIZES,
    alpha_path: Sequence[float] = DEFAULT_ALPHA_PATH,
) -> Path:
    """Ejecuta el Experimento 5 sobre el punto cfg.functional; devuelve out_dir/summary.json."""
    cfg = with_experiment(cfg, EXPERIMENT_ID)
    out_dir = Path(out_dir)
    n, reps = cfg.init.n, cfg.seeds.replicates
    p0 = cfg.functional
    a_hat, g_hat = raw_to_reduced(p0, n)
    delta = cfg.phases.param_perturbation

    # 1. semillas (punto base)
    base = _runs(cfg, p0, 0, reps)
    pcfg0 = with_params(cfg, p0)
    for rep, r in enumerate(base):
        save_run(r, pcfg0, out_dir / "passports", f"exp05_base_rep{rep:03d}", git_commit)

    # 2. tamanos
    by_size: dict[int, list[RunResult]] = {}
    for i, ns in enumerate(sizes):
        c = dataclasses.replace(cfg, init=dataclasses.replace(cfg.init, n=ns))
        by_size[ns] = _runs(c, point_params(a_hat, g_hat, ns), 1 + i, reps)

    # 3. +-10 %
    nb_points = [
        ("alpha_up", a_hat * (1 + delta), g_hat), ("alpha_down", a_hat * (1 - delta), g_hat),
        ("gamma_up", a_hat, g_hat * (1 + delta)), ("gamma_down", a_hat, g_hat * (1 - delta)),
    ]
    neighbors = {
        name: _runs(cfg, point_params(a, g, n), 100 + j, reps) for j, (name, a, g) in enumerate(nb_points)
    }

    # 4. ablacion: mismo alpha crudo; beta -> BETA_OFF y/o gamma -> 0
    variants: dict[str, FunctionalParams] = {
        "-aT": FunctionalParams(alpha=p0.alpha, beta=BETA_OFF, gamma=0.0),
        "-aT+bS_dens": FunctionalParams(alpha=p0.alpha, beta=p0.beta, gamma=0.0),
        "-aT+gS_deg": FunctionalParams(alpha=p0.alpha, beta=BETA_OFF, gamma=p0.gamma),
        "full": p0,
    }
    ablation = {
        name: _group_summary(_runs(cfg, p, 200 + j, reps)) for j, (name, p) in enumerate(variants.items())
    }

    # 5-7. permutacion, nulos y w_min sobre el estado final de la replica 0 (y el de cada replica para nulos)
    w_final = base[0].trajectory.w_final
    perm = random_permutation(n, make_rng(seed_key(cfg.seeds, 300, 0)))
    perm_out = permutation_invariance(w_final, perm, cfg)
    perm_base_w0 = permutation_invariance(base[0].w0, perm, cfg)

    shuffled_f: list[bool] = []
    rewired_f: list[bool] = []
    for rep, r in enumerate(base):
        w = r.trajectory.w_final
        rng = make_rng(seed_key(cfg.seeds, 301, rep))
        sh = shuffled_weight_null(w, cfg.graph.w_min, rng)
        n_edges = int(adjacency(w, cfg.graph.w_min).sum() // 2)
        rw = degree_preserving_null(w, cfg.graph.w_min, 10 * n_edges, rng)
        shuffled_f.append(bool(classify_run(observe(sh, cfg), RunStatus.CONVERGED, cfg.phases).flags["F_candidate"]))
        rewired_f.append(bool(classify_run(observe(rw, cfg), RunStatus.CONVERGED, cfg.phases).flags["F_candidate"]))
    wmin_flags = wmin_candidate_flags(w_final, cfg)

    # 8. histeresis (ascendente desde W0 aleatoria en alpha_min; descendente desde alpha_max)
    path = tuple(float(a) for a in alpha_path)
    w0_rand = random_uniform_weights(n, make_rng(seed_key(cfg.seeds, 400, 0)))
    hyst = {
        "ascending_first": _loop(w0_rand, path, g_hat, cfg),
        "descending_first": _loop(w0_rand, path[::-1], g_hat, cfg),
    }

    # 9. F confirmada
    assessment = confirm_geometric(
        base, by_size, list(neighbors.values()), cfg.phases,
        permutation_max_diff=perm_out["max_abs_diff"], wmin_flags=wmin_flags, null_candidates=shuffled_f,
    )

    summary = {
        "experiment": EXPERIMENT_ID,
        "description": "robustez: semillas, tamanos, +-10%, ablacion, permutacion, nulos, histeresis (M§29-M§33)",
        "n_nodes": n,
        "replicates": reps,
        "point": {"alpha_hat": a_hat, "gamma_hat": g_hat, "params_raw": dataclasses.asdict(p0)},
        "seeds": {**_group_summary(base), "runs": [run_summary(r) for r in base]},
        "sizes": {str(ns): _group_summary(rs) for ns, rs in by_size.items()},
        "perturbations": {name: _group_summary(rs) for name, rs in neighbors.items()},
        "ablation": ablation,
        "ablation_beta_off": BETA_OFF,
        "permutation": {
            "final_state": perm_out,
            "initial_state": perm_base_w0,
            "exact_within_1e-9": bool(perm_out["max_abs_diff"] <= 1e-9 and perm_base_w0["max_abs_diff"] <= 1e-9),
        },
        "nulls": {
            "shuffled_weights_f_candidate": shuffled_f,
            "degree_preserving_f_candidate": rewired_f,
        },
        "w_min_f_candidate_flags": list(wmin_flags),
        "hysteresis": hyst,
        "f_confirmation": {"label": assessment.label, "flags": dict(assessment.flags)},
    }
    path_out = out_dir / "summary.json"
    path_out.write_text(json.dumps(to_jsonable(summary), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return path_out


def test_smoke(tmp_path: Path) -> None:
    n = 16
    cfg = default_config(n, master_entropy=105, replicates=2)
    cfg = dataclasses.replace(
        cfg, dynamics=DynamicsConfig(max_steps=20_000), functional=point_params(3.0, 1.0, n)
    )
    data = json.loads(
        run(cfg, tmp_path, "test", sizes=(12, 16), alpha_path=(0.0, 1.5, 3.0)).read_text(encoding="utf-8")
    )
    assert data["seeds"]["labels"] == ["E", "E"]
    assert set(data["sizes"]) == {"12", "16"}
    assert set(data["perturbations"]) == {"alpha_up", "alpha_down", "gamma_up", "gamma_down"}
    assert set(data["ablation"]) == {"-aT", "-aT+bS_dens", "-aT+gS_deg", "full"}
    assert data["ablation"]["-aT"]["labels"] == ["E", "E"]
    assert data["permutation"]["exact_within_1e-9"] is True
    assert len(data["nulls"]["shuffled_weights_f_candidate"]) == 2
    assert len(data["w_min_f_candidate_flags"]) == 5
    h = data["hysteresis"]
    assert h["ascending_first"]["hysteresis_detected"] is False
    assert h["descending_first"]["hysteresis_detected"] is True  # biestabilidad: E en [1,3] frente a A
    assert data["f_confirmation"]["flags"]["F_confirmed"] is False
    assert data["f_confirmation"]["label"] != "F"
    assert len(list((tmp_path / "passports").glob("*.json"))) == 2


@pytest.mark.slow
def test_full_run(tmp_path: Path) -> None:
    cfg = default_config(200, master_entropy=20240901, replicates=30)
    cfg = dataclasses.replace(cfg, functional=point_params(1.0, 1.0, 200))
    data = json.loads(run(cfg, tmp_path).read_text(encoding="utf-8"))
    assert set(data["sizes"]) == {"100", "200", "300"}


