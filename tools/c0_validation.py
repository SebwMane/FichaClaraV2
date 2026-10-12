"""C0 (Rev. 3 §4): validacion de instrumentos causales. No simula dinamica ni dice nada sobre Omega.

Modos:
  smoke : mallas minimas, sin decision (comprueba que todo corre).
  pilot : entropia 20261012, declarada; solo informa d_MM - d_mid para fijar tau de C0-4.
  full  : entropia 20261011; criterios congelados en docs/OMEGA_2_PROPUESTA_CAUSAL.md; exige arbol limpio.
Semillas: PCG64(SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, celda, semilla))).
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import (  # noqa: E402
    between_counts,
    bd_action_2d,
    midpoint_dimension,
    mm_dimension,
)
from omega.causality.orders import kr_order, transitive_percolation  # noqa: E402
from omega.causality.sprinkling import (  # noqa: E402
    causal_matrix,
    diamond_volume,
    interval_volume,
    sprinkle_diamond,
)

EXPERIMENT_ID = 300
ENTROPY = {"full": 20261011, "pilot": 20261012, "smoke": 20261013}
PERC_P = 0.05


def rng_for(entropy: int, cell: int, seed: int) -> np.random.Generator:
    ss = np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, cell, seed))
    return np.random.Generator(np.random.PCG64(ss))


def git_state() -> tuple[str, bool]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True).stdout.strip())
    return head, dirty


def c0_1(entropy: int, ns: list[int], seeds: int) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for d in (2, 3, 4):
        for n in ns:
            vals = [mm_dimension(causal_matrix(sprinkle_diamond(n, d, rng_for(entropy, 10 * d + 1, s)))) for s in range(seeds)]
            out[f"d{d}_N{n}"] = {"mean": float(np.mean(vals)), "sd": float(np.std(vals, ddof=1)), "values": vals}
    return out


def c0_2(entropy: int, n: int, seeds: int, pairs: int) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for d in (2, 3, 4):
        obs_tot = exp_tot = 0.0
        for s in range(seeds):
            rng = rng_for(entropy, 20 * d + 2, s)
            x = sprinkle_diamond(n, d, rng)
            c = causal_matrix(x)
            b = between_counts(c)
            ii, jj = np.nonzero(c)
            pick = rng.choice(len(ii), size=min(pairs, len(ii)), replace=False)
            for k in pick:
                i, j = int(ii[k]), int(jj[k])
                dt = x[j, 0] - x[i, 0]
                tau = math.sqrt(max(dt * dt - float(np.sum((x[j, 1:] - x[i, 1:]) ** 2)), 0.0))
                obs_tot += float(b[i, j])
                exp_tot += (n - 2) * interval_volume(tau, d) / diamond_volume(d)
        out[f"d{d}"] = {"ratio": obs_tot / exp_tot, "observed": obs_tot, "expected": exp_tot}
    return out


def c0_3(entropy: int, ns: list[int], seeds: int) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for n in ns:
        vals = [bd_action_2d(causal_matrix(sprinkle_diamond(n, 2, rng_for(entropy, 3000 + n, s)))) for s in range(seeds)]
        out[f"N{n}"] = {"mean": float(np.mean(vals)), "se": float(np.std(vals, ddof=1) / math.sqrt(len(vals))), "values": vals}
    return out


def c0_4(entropy: int, n: int, seeds: int) -> dict[str, Any]:
    out: dict[str, Any] = {}
    gens = {
        "sprinkle_d2": lambda r: causal_matrix(sprinkle_diamond(n, 2, r)),
        "sprinkle_d3": lambda r: causal_matrix(sprinkle_diamond(n, 3, r)),
        "sprinkle_d4": lambda r: causal_matrix(sprinkle_diamond(n, 4, r)),
        "kr": lambda r: kr_order(n, r),
        "percolation": lambda r: transitive_percolation(n, PERC_P, r),
    }
    for cell, (name, gen) in enumerate(gens.items()):
        rows = []
        for s in range(seeds):
            c = gen(rng_for(entropy, 4000 + cell, s))
            dmm, dmid = mm_dimension(c), midpoint_dimension(c)
            rows.append({"d_mm": dmm, "d_mid": dmid, "gap": abs(dmm - dmid), "S_bd": bd_action_2d(c)})
        out[name] = rows
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "pilot", "full"], required=True)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    mode = args.mode
    head, dirty = git_state()
    if mode == "full" and dirty:
        sys.exit("modo full exige arbol git limpio")
    out = args.out or Path("results") / ("c0_validation" if mode == "full" else f"c0_{mode}")
    out.mkdir(parents=True, exist_ok=True)
    ent = ENTROPY[mode]
    summary: dict[str, Any] = {"phase": "preregistered", "mode": mode, "commit": head, "dirty": dirty, "entropy": ent}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    t0 = time.time()
    if mode == "smoke":
        res = {"c0_1": c0_1(ent, [128], 3), "c0_2": c0_2(ent, 128, 2, 50), "c0_3": c0_3(ent, [64], 10), "c0_4": c0_4(ent, 128, 2)}
    elif mode == "pilot":
        res = {"c0_4": c0_4(ent, 1024, 10)}
    else:
        res = {
            "c0_1": c0_1(ent, [256, 1024], 20),
            "c0_2": c0_2(ent, 1024, 20, 200),
            "c0_3": c0_3(ent, [64, 256, 1024], 200),
            "c0_4": c0_4(ent, 1024, 20),
        }
    with (out / "runs.jsonl").open("w") as fh:
        fh.write(json.dumps(res) + "\n")
    summary.update({"phase": "final", "seconds": round(time.time() - t0, 1)})
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
