"""C0-4b (Rev. 3 §9): discriminador por perfil de abundancias de intervalos.

Constantes congeladas en docs/OMEGA_2_PROPUESTA_CAUSAL.md §9. No ajustar.
Modos: smoke (N pequeno, sin veredicto) y full (N = 1024, arbol limpio, veredicto calculado).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import interval_abundances, mm_dimension  # noqa: E402
from omega.causality.orders import kr_order, thick_layer_order, transitive_percolation  # noqa: E402
from omega.causality.profile import K_MAX, d_ref, profile_distance, reference_abundances  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 301
ENT_REF, ENT_CAL, ENT_OFF = 20261015, 20261016, 20261017
N_REF = N_CAL = 40
SEEDS = 20
CAL_FACTOR, CAL_PCT = 1.5, 95.0
PASS_POS, MAX_CTRL = 0.90, 0.10

Gen = Callable[[np.random.Generator], npt.NDArray[np.bool_]]


def rng_for(entropy: int, cell: int, seed: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, cell, seed))))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def references(n: int, reps: int) -> dict[int, npt.NDArray[np.float64]]:
    refs = {d: reference_abundances(d, n, K_MAX, [rng_for(ENT_REF, d, s) for s in range(reps)]) for d in (3, 4)}
    refs[1] = reference_abundances(1, n, K_MAX)
    refs[2] = reference_abundances(2, n, K_MAX)
    return refs


def classify(c: npt.NDArray[np.bool_], refs: dict[int, npt.NDArray[np.float64]], delta: dict[int, float]) -> dict[str, Any]:
    d_hat = mm_dimension(c)
    dr = d_ref(d_hat)
    if dr is None:
        return {"d_hat": d_hat, "d_ref": None, "D": None, "flag": False}
    dist = profile_distance(interval_abundances(c, K_MAX), refs[dr])
    return {"d_hat": d_hat, "d_ref": dr, "D": dist, "flag": bool(dist <= delta[dr])}


def calibrate(n: int, reps: int, refs: dict[int, npt.NDArray[np.float64]]) -> tuple[dict[int, float], dict[int, list[float]]]:
    raw: dict[int, list[float]] = {}
    delta: dict[int, float] = {}
    for d in (2, 3, 4):
        raw[d] = [
            profile_distance(interval_abundances(causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_CAL, d, s))), K_MAX), refs[d])
            for s in range(reps)
        ]
        delta[d] = CAL_FACTOR * float(np.percentile(raw[d], CAL_PCT))
    delta[1] = delta[2]
    return delta, raw


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], required=True)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    head, dirty = git_state()
    if args.mode == "full" and dirty:
        sys.exit("modo full exige arbol git limpio")
    n, reps, seeds = (1024, N_REF, SEEDS) if args.mode == "full" else (128, 4, 3)
    out = args.out or Path("results") / ("c0_4b_validation" if args.mode == "full" else "c0_4b_smoke")
    out.mkdir(parents=True, exist_ok=True)
    summary: dict[str, Any] = {"phase": "preregistered", "mode": args.mode, "commit": head, "dirty": dirty, "N": n}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    t0 = time.time()
    refs = references(n, reps)
    delta, cal_raw = calibrate(n, reps, refs)
    cells: dict[str, Gen] = {
        "sprinkle_d2": lambda r: causal_matrix(sprinkle_diamond(n, 2, r)),
        "sprinkle_d3": lambda r: causal_matrix(sprinkle_diamond(n, 3, r)),
        "sprinkle_d4": lambda r: causal_matrix(sprinkle_diamond(n, 4, r)),
        "kr": lambda r: kr_order(n, r),
        "percolation_p05": lambda r: transitive_percolation(n, 0.05, r),
        "percolation_p01": lambda r: transitive_percolation(n, 0.01, r),
        "thick_layers_8": lambda r: thick_layer_order(n, 8, 0.5, r),
    }
    rows: dict[str, list[dict[str, Any]]] = {}
    for i, (name, gen) in enumerate(cells.items()):
        rows[name] = [classify(gen(rng_for(ENT_OFF, 100 + i, s)), refs, delta) for s in range(seeds)]
    frac = {k: float(np.mean([r["flag"] for r in v])) for k, v in rows.items()}
    pos = [k for k in rows if k.startswith("sprinkle")]
    ctrl = [k for k in rows if not k.startswith("sprinkle")]
    with (out / "runs.jsonl").open("w") as fh:
        fh.write(json.dumps({"delta": delta, "calibration_D": cal_raw, "references": {k: v.tolist() for k, v in refs.items()}, "rows": rows}) + "\n")
    summary.update({"phase": "final", "seconds": round(time.time() - t0, 1), "delta": delta, "flag_fraction": frac})
    if args.mode == "full":
        summary["verdict_C0-4b"] = all(frac[k] >= PASS_POS for k in pos) and all(frac[k] <= MAX_CTRL for k in ctrl)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
