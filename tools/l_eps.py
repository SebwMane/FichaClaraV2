"""L-epsilon (Rev. 3 §18): dispersion de S_eps/eps^2 en 2-ordenes con eps = K/N. Constantes congeladas."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import between_counts, smearing_weight  # noqa: E402
from omega.causality.orders import random_2order  # noqa: E402

EXPERIMENT_ID = 307
ENT = 20261037
KS = (4, 8)
NS = (64, 128, 256, 512, 1024)
M = 200
D_MAX = 3 / math.pi
OUT = Path("results") / "l_eps"


def main() -> None:
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    d: dict[int, dict[int, float]] = {k: {} for k in KS}
    for i, n in enumerate(NS):
        vals: dict[int, list[float]] = {k: [] for k in KS}
        for s in range(M):
            rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENT, spawn_key=(EXPERIMENT_ID, i, s))))
            b = between_counts(random_2order(n, rng))
            v = b[b >= 0].astype(np.float64)
            for k in KS:
                eps = k / n
                s_eps = eps * (n - 2 * eps * float(np.sum(smearing_weight(v, eps))))
                vals[k].append(s_eps / eps**2)
        for k in KS:
            d[k][n] = float(np.std(vals[k], ddof=1))
    fits: dict[str, Any] = {}
    for k in KS:
        x = np.log(NS)
        y = np.log([d[k][n] for n in NS])
        a = np.vstack([x, np.ones_like(x)]).T
        coef = np.linalg.lstsq(a, y, rcond=None)[0]
        r = y - a @ coef
        se = math.sqrt(float(r @ r) / (len(x) - 2) * float(np.linalg.inv(a.T @ a)[0, 0]))
        fits[str(k)] = {"eta": float(coef[0]), "se": se, "D": {str(n): d[k][n] for n in NS}}
    compat = [k for k in KS if d[k][1024] <= D_MAX and fits[str(k)]["eta"] + 2 * fits[str(k)]["se"] <= 0]
    incompat = all(fits[str(k)]["eta"] - 2 * fits[str(k)]["se"] > 0 for k in KS)
    verdict = "COMPATIBLE" if compat else ("INCOMPATIBLE (estructural)" if incompat else "INDETERMINADO")
    summary = {"commit": head, "seconds": round(time.time() - t0, 1), "fits": fits, "L-eps": verdict}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
