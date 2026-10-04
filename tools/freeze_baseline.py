"""Congela el baseline dorado de Omega-1.0 (OMEGA_1_1_DESIGN §2 'Golden baseline').

Se ejecuta UNA vez sobre el arbol intacto de 9dfbea7 y escribe `tests/expected_baseline/*.json`.
`tests/test_baseline_golden.py` reutiliza `compute_fixtures` y `compute_manifest` para recalcular.
Uso: python tools/freeze_baseline.py
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from omega.config.convert import reduced_to_raw  # noqa: E402
from omega.config.seeds import SeedKey, make_rng, seed_key  # noqa: E402
from omega.config.settings import DynamicsConfig, OmegaConfig  # noqa: E402
from omega.dynamics.evolution import evolve  # noqa: E402
from omega.experiments.reference_graphs import periodic_lattice  # noqa: E402
from omega.io.passport import array_digest, build_passport  # noqa: E402
from omega.network.initialization import random_uniform_weights  # noqa: E402
from omega.network.weights import upper_triangle  # noqa: E402
from omega.phases.scan import default_config, observe, simulate  # noqa: E402
from omega.types import DimensionEstimate, RunResult  # noqa: E402

BASELINE_COMMIT = "9dfbea7d6be28a8f17fbf336579c23628ac5669d"
OUT_DIR = ROOT / "tests" / "expected_baseline"
FIXTURE_NAMES = (
    "empty_phase",
    "complete_phase",
    "uniform_barrier",
    "geometry_reference",
    "passport_keys",
)
MANIFEST_NAME = "omega10_manifest"


def plain(obj: Any) -> Any:
    """Tipos JSON puros; floats no finitos -> cadenas 'nan'/'inf'/'-inf' (comparables exactamente)."""
    if isinstance(obj, dict):
        return {str(k): plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return plain(obj.tolist())
    if isinstance(obj, np.generic):
        return plain(obj.item())
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else ("nan" if math.isnan(obj) else ("inf" if obj > 0 else "-inf"))
    return obj


def _run_record(r: RunResult) -> dict[str, Any]:
    t, top = r.trajectory, r.observables.topology
    act = t.scalars["action"]
    return {
        "status": t.status.value,
        "steps": t.steps,
        "dt": t.dt,
        "tau": t.tau,
        "label": r.assessment.label.value,
        "flags": {k: bool(v) for k, v in r.assessment.flags.items()},
        "frac_at_zero": top.frac_at_zero,
        "frac_at_one": top.frac_at_one,
        "mean_strength": top.mean_strength,
        "action": [float(act[i]) for i in (0, 1, 2, -1)],
        "max_dw_last": float(t.scalars["max_dw"][-1]),
        "w0_digest": array_digest(upper_triangle(r.w0)),
        "w_final_upper": [float(x) for x in upper_triangle(t.w_final)],
    }


def _phase_fixture(alpha_hat: float, gamma_hat: float) -> dict[str, Any]:
    n = 40
    cfg = default_config(n, master_entropy=20261004, experiment_id=900, replicates=3)
    p = reduced_to_raw(alpha_hat, gamma_hat, n)
    runs = [_run_record(simulate(cfg, p, seed_key(cfg.seeds, 0, rep))) for rep in range(3)]
    return {"n": n, "alpha_hat": alpha_hat, "gamma_hat": gamma_hat, "runs": runs}


def _uniform_fixture() -> dict[str, Any]:
    n = 40
    cfg = default_config(n, master_entropy=20261004, experiment_id=950, replicates=3)
    dyn = DynamicsConfig()
    homogeneous: list[dict[str, Any]] = []
    ones = np.ones((n, n)) - np.eye(n)
    for a_hat in (1.0, 1.5, 2.0):
        p = reduced_to_raw(a_hat, 0.0, n)
        w_star = 1.0 / a_hat
        for tag, w0v in (("below", w_star - 0.01), ("above", min(w_star + 0.01, 1.0))):
            traj = evolve(w0v * ones, p, dyn)
            up = upper_triangle(traj.w_final)
            act = traj.scalars["action"]
            homogeneous.append(
                {
                    "alpha_hat": a_hat,
                    "start": tag,
                    "w0": w0v,
                    "steps": traj.steps,
                    "status": traj.status.value,
                    "w_final_mean": float(up.mean()),
                    "w_final_min": float(up.min()),
                    "w_final_max": float(up.max()),
                    "action": [float(act[0]), float(act[1]), float(act[-1])],
                }
            )
    random_start = []
    for a_hat in (1.9, 2.1):
        r = simulate(cfg, reduced_to_raw(a_hat, 0.0, n), seed_key(cfg.seeds, 0, 0))
        random_start.append({"alpha_hat": a_hat, "label": r.assessment.label.value})
    return {"n": n, "gamma_hat": 0.0, "homogeneous": homogeneous, "random_start": random_start}


def _dim(e: DimensionEstimate) -> dict[str, Any]:
    return {
        "value": e.value,
        "stderr": e.stderr,
        "status": e.status,
        "window": None if e.window is None else list(e.window),
        "plateau": e.plateau,
        "method": e.method,
    }


def _observe_record(w: Any, cfg: OmegaConfig) -> dict[str, Any]:
    obs = observe(w, cfg)
    geo = obs.geometry
    return {
        "topology": dataclasses.asdict(obs.topology),
        "path_length": geo.path_length,
        "path_length_hops": geo.path_length_hops,
        "diameter": geo.diameter,
        "d_eff": _dim(geo.d_eff),
        "d_eff_ball": _dim(geo.d_eff_ball),
        "d_eff_hops": _dim(geo.d_eff_hops),
        "d_s": _dim(geo.d_s),
    }


def _geometry_fixture() -> dict[str, Any]:
    lattice = periodic_lattice((7, 7, 7))
    key = SeedKey(20261004, (901, 0, 0))
    w_rand = random_uniform_weights(60, make_rng(key))
    return {
        "lattice_7x7x7": _observe_record(lattice, default_config(343)),
        "uniform_n60": _observe_record(w_rand, default_config(60)),
    }


def _passport_fixture() -> dict[str, Any]:
    n = 40
    cfg = default_config(n, master_entropy=20261004, experiment_id=900, replicates=3)
    r = simulate(cfg, reduced_to_raw(1.0, 1.0, n), seed_key(cfg.seeds, 0, 0))
    pp = build_passport(cfg, r, {"python": "x"}, "unknown")
    return {"top_level": sorted(pp), "observables": sorted(pp["observables"])}


def compute_fixtures() -> dict[str, dict[str, Any]]:
    """Recalcula todos los fixtures (sin el manifiesto) como JSON puro."""
    raw = {
        "empty_phase": _phase_fixture(1.0, 1.0),
        "complete_phase": _phase_fixture(3.0, 1.0),
        "uniform_barrier": _uniform_fixture(),
        "geometry_reference": _geometry_fixture(),
        "passport_keys": _passport_fixture(),
    }
    out: dict[str, dict[str, Any]] = {k: plain(v) for k, v in raw.items()}
    return out


def manifest_paths(root: Path = ROOT) -> list[str]:
    """Archivos de Omega-1.0 cubiertos: omega/**/*.py, tests/*.py y pyproject.toml."""
    files = {p for p in (root / "omega").rglob("*.py") if "__pycache__" not in p.parts}
    files |= set((root / "tests").glob("*.py"))
    files.add(root / "pyproject.toml")
    return sorted(p.relative_to(root).as_posix() for p in files)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compute_manifest(root: Path = ROOT) -> dict[str, str]:
    return {rel: sha256_file(root / rel) for rel in manifest_paths(root)}


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fixtures = compute_fixtures()
    payloads: dict[str, dict[str, Any]] = dict(fixtures)
    payloads[MANIFEST_NAME] = {"sha256": compute_manifest()}
    for name, body in payloads.items():
        doc = {"generated_from_commit": BASELINE_COMMIT, **body}
        text = json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=False, allow_nan=False)
        (OUT_DIR / f"{name}.json").write_text(text + "\n", encoding="utf-8")
        print(f"escrito {name}.json")


if __name__ == "__main__":
    main()
