"""F3 y F4 (docs/OMEGA_F3F4_PREREGISTRO.md). Partes: f3cal | f3_<estructura> | f4_<d> | final."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from f1f2 import random_dorder  # noqa: E402

from omega.causality.entropy import knee_from_spectrum  # noqa: E402
from omega.causality.estimators import between_counts  # noqa: E402
from omega.causality.orders import kr_order, thick_layer_order, transitive_percolation  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402
from omega.causality.torus import slab_causal_matrix, sprinkle_slab  # noqa: E402

EXPERIMENT_ID = 316
ENT_CAL, ENT_TEST, ENT_F4 = 20261063, 20261064, 20261065
OUT = Path("results") / "omega_f3f4"
Bool = npt.NDArray[np.bool_]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def links(c: Bool) -> Bool:
    out: Bool = c & (between_counts(c) == 0)
    return out


def spectrum(x: Bool) -> npt.NDArray[np.float64]:
    xf = x.astype(np.float64)
    mu: npt.NDArray[np.float64] = np.asarray(np.linalg.eigvalsh(0.5j * (xf - xf.T)), dtype=np.float64)
    return mu


def sig(c: Bool) -> dict[str, Any]:
    n = c.shape[0]
    kl, bl = knee_from_spectrum(spectrum(links(c)), n)
    return {"rodilla_L": kl, "beta_L": bl}


def de_sitter_4d(n: int, r: np.random.Generator) -> Bool:
    pts: list[npt.NDArray[np.float64]] = []
    total = 0
    while total < n:
        cand = sprinkle_diamond(4 * n, 4, r)
        eta = cand[:, 0] - 1.0
        keep = cand[r.random(cand.shape[0]) < (0.5 / np.abs(eta)) ** 4]
        pts.append(keep)
        total += keep.shape[0]
    allp = np.concatenate(pts)[:n]
    return causal_matrix(allp)


F3_STRUCTS: dict[str, tuple[str, Callable[[int, np.random.Generator], Bool]]] = {
    "diamante_4d": ("positiva", lambda n, r: causal_matrix(sprinkle_diamond(n, 4, r))),
    "bloque_4d": ("positiva", lambda n, r: slab_causal_matrix(sprinkle_slab(n, r))),
    "de_sitter_4d": ("positiva", de_sitter_4d),
    "kr": ("negativa", lambda n, r: kr_order(n, r)),
    "orden_d4": ("negativa", random_dorder(4)),
    "percolacion_005": ("negativa", lambda n, r: transitive_percolation(n, 0.05, r)),
    "capas_gruesas": ("negativa", lambda n, r: thick_layer_order(n, 8, 0.5, r)),
}


def windows(cal: list[dict[str, Any]]) -> dict[str, tuple[float, float]]:
    b = np.array([r["beta_L"] for r in cal])
    k = np.array([r["rodilla_L"] for r in cal if r["rodilla_L"] is not None])
    mb = float(np.median(b))
    mad = float(np.median(np.abs(b - mb)))
    hw = max(0.05, 6 * mad)
    mk = float(np.median(k)) if k.size else float("nan")
    return {"beta": (mb - hw, mb + hw), "rodilla": (0.5 * mk, 2 * mk)}


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
    if args.part == "f3cal":
        for n in (1600, 2400):
            for s in range(5):
                rows.append({"N": n, "seed": s, **sig(causal_matrix(sprinkle_diamond(n, 4, rng_for(ENT_CAL, n, s))))})
    elif args.part.startswith("f3_"):
        name = args.part[3:]
        role, gen = F3_STRUCTS[name]
        k = list(F3_STRUCTS).index(name)
        for n in (1600, 2400):
            for s in range(5):
                rows.append({"estructura": name, "papel": role, "N": n, "seed": s, **sig(gen(n, rng_for(ENT_TEST, k, n, s)))})
    elif args.part.startswith("f4_"):
        d = int(args.part[3:])
        for n in (400, 800, 1600, 2400):
            for s in range(5):
                c = causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_F4, d, n, s)))
                _, bc = knee_from_spectrum(spectrum(c), n)
                rows.append({"d": d, "N": n, "seed": s, **sig(c), "beta_C": bc})
    elif args.part == "final":
        cal = json.loads((OUT / "f3cal.json").read_text())["rows"]
        win = windows(cal)
        res3: dict[str, Any] = {}
        for name, (role, _) in F3_STRUCTS.items():
            rr = json.loads((OUT / f"f3_{name}.json").read_text())["rows"]
            b = float(np.nanmedian([r["beta_L"] for r in rr]))
            ks = [r["rodilla_L"] for r in rr if r["rodilla_L"] is not None]
            km = float(np.median(ks)) if ks else None
            comp = bool(win["beta"][0] <= b <= win["beta"][1] and km is not None and win["rodilla"][0] <= km <= win["rodilla"][1])
            res3[name] = {"papel": role, "beta_L": b, "rodilla_L": km, "con_rodilla": len(ks), "compatible": comp,
                          "correcto": comp if role == "positiva" else not comp}
        f3 = "PASA" if all(v["correcto"] for v in res3.values()) else "FALLA"
        res4: dict[str, Any] = {}
        for d in (2, 3, 4, 5):
            rr = json.loads((OUT / f"f4_{d}.json").read_text())["rows"]
            mL = {n: float(np.nanmedian([r["beta_L"] for r in rr if r["N"] == n])) for n in (400, 800, 1600, 2400)}
            mC = {n: float(np.nanmedian([r["beta_C"] for r in rr if r["N"] == n])) for n in (400, 800, 1600, 2400)}
            res4[str(d)] = {"beta_L": mL, "beta_C": mC, "Q1_estable": bool(abs(mL[2400] - mL[800]) <= 0.05)}
        q2 = bool(res4["3"]["beta_L"][2400] > res4["4"]["beta_L"][2400] > res4["5"]["beta_L"][2400])
        out = {"commit": head, "F3": f3, "ventanas": win, "F3_detalle": res3, "F4": res4, "Q2_monotonia": q2}
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print(json.dumps({"F3": f3, "Q1": {d: v["Q1_estable"] for d, v in res4.items()}, "Q2": q2}))
        return
    else:
        sys.exit("parte desconocida")
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "rows": rows}, default=str))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
