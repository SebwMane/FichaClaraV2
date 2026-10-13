"""Consejo Rev. 7: H1 (propagador), H2 (tramo), H3 (borde) y D-3D. Partes: h1_<N> | h3 | d3 | final."""

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
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.estimators import between_counts  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402
from omega.causality.torus import slab_causal_matrix, sprinkle_slab  # noqa: E402

EXPERIMENT_ID = 315
ENT_H1, ENT_H3, ENT_D3 = 20261060, 20261061, 20261062
OUT = Path("results") / "omega_rev7"
Bool = npt.NDArray[np.bool_]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def links(c: Bool) -> Bool:
    out: Bool = c & (between_counts(c) == 0)
    return out


def beta_fit(x: Bool, n: int, k_ref: int) -> float:
    xf = x.astype(np.float64)
    mu = np.linalg.eigvalsh(0.5j * (xf - xf.T))
    pos = np.sort(mu[mu > 1e-12])[::-1]
    if pos.size <= k_ref:
        return float("nan")
    ranks = np.arange(1, pos.size + 1, dtype=np.float64)
    b, _ = np.polyfit(np.log(ranks[1:k_ref]), np.log(pos[1:k_ref]), 1)
    return float(b)


def tramos(n: int) -> dict[str, int]:
    return {"K_8": max(6, math.ceil(math.sqrt(n) / 8)), "K_4": max(6, math.ceil(math.sqrt(n) / 4)), "K_2": max(6, math.ceil(math.sqrt(n) / 2))}


def measure(c: Bool, n: int) -> dict[str, float]:
    row: dict[str, float] = {}
    lk = links(c)
    for name, mat in (("C", c), ("L", lk)):
        for tk, k in tramos(n).items():
            row[f"{name}_{tk}"] = beta_fit(mat, n, k)
    return row


def load(part: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = json.loads((OUT / f"{part}.json").read_text())["rows"]
    return rows


def med(rows: list[dict[str, Any]], key: str) -> float:
    return float(np.nanmedian([r[key] for r in rows]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    if args.part.startswith("h1_"):
        n = int(args.part[3:])
        for s in range(5):
            rows.append({"N": n, "seed": s, **measure(causal_matrix(sprinkle_diamond(n, 4, rng_for(ENT_H1, n, s))), n)})
    elif args.part == "h3":
        for s in range(5):
            rows.append({"N": 2400, "seed": s, **measure(slab_causal_matrix(sprinkle_slab(2400, rng_for(ENT_H3, s))), 2400)})
    elif args.part == "d3":
        for n in (800, 1600, 2400):
            for s in range(5):
                rows.append({"N": n, "seed": s, **measure(causal_matrix(sprinkle_diamond(n, 3, rng_for(ENT_D3, n, s))), n)})
    elif args.part == "final":
        h1 = [r for n in (800, 1600, 2400, 3200) for r in load(f"h1_{n}")]
        tab = {str(n): {k: med([r for r in h1 if r["N"] == n], k) for k in h1[0] if k not in ("N", "seed")} for n in (800, 1600, 2400, 3200)}
        bL, bC = tab["3200"]["L_K_4"], tab["3200"]["C_K_4"]
        h1_ok = abs(bL + 0.5) <= 0.10 and abs(bL + 0.5) < abs(bC + 0.5)
        spread = {f"{n}_{m}": max(tab[n][f"{m}_K_8"], tab[n][f"{m}_K_4"], tab[n][f"{m}_K_2"]) - min(tab[n][f"{m}_K_8"], tab[n][f"{m}_K_4"], tab[n][f"{m}_K_2"]) for n in tab for m in ("C", "L")}
        h2_ok = any(v > 0.10 for v in spread.values())
        h3 = load("h3")
        slab = {m: med(h3, f"{m}_K_4") for m in ("C", "L")}
        diam = {m: tab["2400"][f"{m}_K_4"] for m in ("C", "L")}
        h3_ok = any(abs(slab[m] - diam[m]) > 0.10 for m in ("C", "L"))
        d3 = load("d3")
        d3tab = {str(n): {m: med([r for r in d3 if r["N"] == n], f"{m}_K_4") for m in ("C", "L")} for n in (800, 1600, 2400)}
        out = {"commit": head, "H1_apoyada": h1_ok, "H2_apoyada": h2_ok, "H3_apoyada": h3_ok, "tabla_4D": tab,
               "dispersion_por_tramo": spread, "bloque_vs_diamante_2400": {"bloque": slab, "diamante": diam}, "D3": d3tab}
        (OUT / "summary.json").write_text(json.dumps(out, indent=2))
        print(json.dumps({k: out[k] for k in ("H1_apoyada", "H2_apoyada", "H3_apoyada")}))
        return
    else:
        sys.exit("parte desconocida")
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "rows": rows}))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
