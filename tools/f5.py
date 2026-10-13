"""F5 (docs/OMEGA_CONSEJO_REV8_F5.md §3 y §6). Partes: ref_<d> | cal_<d> | test_<estructura> | final."""

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
from f3f4 import de_sitter_4d, links, spectrum, windows  # noqa: E402

from omega.causality.entropy import knee_from_spectrum  # noqa: E402
from omega.causality.estimators import interval_abundances, midpoint_dimension, mm_dimension  # noqa: E402
from omega.causality.orders import kr_order, thick_layer_order, transitive_percolation  # noqa: E402
from omega.causality.profile import K_MAX, profile_distance  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402
from omega.causality.torus import slab_causal_matrix, sprinkle_slab  # noqa: E402

EXPERIMENT_ID = 317
ENT_REF, ENT_CAL, ENT_TEST = 20261068, 20261066, 20261067
NS = (1600, 2400)
DIMS = (2, 3, 4, 5)
F3_WIN = {"beta": (-0.363, -0.263), "rodilla": (1.511, 6.046)}
OUT = Path("results") / "omega_f5"
Bool = npt.NDArray[np.bool_]

TEST: dict[str, tuple[str, int, Callable[[int, np.random.Generator], Bool]]] = {
    "diamante_4d": ("positiva", 4, lambda n, r: causal_matrix(sprinkle_diamond(n, 4, r))),
    "bloque_4d": ("positiva", 4, lambda n, r: slab_causal_matrix(sprinkle_slab(n, r))),
    "de_sitter_4d": ("positiva", 4, de_sitter_4d),
    "diamante_3d": ("positiva", 3, lambda n, r: causal_matrix(sprinkle_diamond(n, 3, r))),
    "producto_3": ("negativa", 0, random_dorder(3)),
    "producto_4": ("negativa", 0, random_dorder(4)),
    "producto_5": ("negativa", 0, random_dorder(5)),
    "kr": ("negativa", 0, lambda n, r: kr_order(n, r)),
    "percolacion_005": ("negativa", 0, lambda n, r: transitive_percolation(n, 0.05, r)),
    "capas_gruesas": ("negativa", 0, lambda n, r: thick_layer_order(n, 8, 0.5, r)),
}


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def measure(c: Bool) -> dict[str, Any]:
    n = c.shape[0]
    kl, bl = knee_from_spectrum(spectrum(links(c)), n)
    return {"N": n, "beta_L": bl, "rodilla_L": kl, "d_mm": mm_dimension(c), "d_mid": midpoint_dimension(c),
            "abund": [int(v) for v in interval_abundances(c, K_MAX)]}


def load(part: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = json.loads((OUT / f"{part}.json").read_text())["rows"]
    return rows


def d_spec(beta: float, curve: dict[int, float]) -> float:
    ds = list(DIMS)
    m = [curve[d] for d in ds]
    for i in range(len(ds) - 1):
        lo, hi = m[i], m[i + 1]
        if (beta <= lo and beta >= hi) or i == len(ds) - 2 or (i == 0 and beta > lo):
            return ds[i] + (beta - lo) / (hi - lo) * (ds[i + 1] - ds[i])
    return float("nan")


def mad(x: list[float]) -> float:
    a = np.array(x, dtype=float)
    return float(np.median(np.abs(a - np.median(a))))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    if "R8-1 | **Firmada" not in (Path("docs") / "OMEGA_CONSEJO_REV8_F5.md").read_text():
        sys.exit("F5 requiere la firma R8-1")
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    if args.part.startswith("ref_"):
        d = int(args.part[4:])
        for n in NS:
            ab = [interval_abundances(causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_REF, d, n, s))), K_MAX) for s in range(20)]
            rows.append({"d": d, "N": n, "ref": np.mean(np.array(ab, dtype=float), axis=0).tolist()})
    elif args.part.startswith("cal_"):
        d = int(args.part[4:])
        for n in NS:
            for s in range(10):
                rows.append({"d": d, "seed": s, **measure(causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_CAL, d, n, s))))})
    elif args.part.startswith("test_"):
        name = args.part[5:]
        role, dtrue, gen = TEST[name]
        k = list(TEST).index(name)
        for n in NS:
            for s in range(5):
                rows.append({"estructura": name, "papel": role, "d_verdadera": dtrue, "seed": s, **measure(gen(n, rng_for(ENT_TEST, k, n, s)))})
    elif args.part == "final":
        refs = {(r["d"], r["N"]): np.array(r["ref"]) for d in DIMS for r in load(f"ref_{d}")}
        cal = {d: load(f"cal_{d}") for d in DIMS}
        curve = {n: {d: float(np.nanmedian([r["beta_L"] for r in cal[d] if r["N"] == n])) for d in DIMS} for n in NS}

        def dist(r: dict[str, Any], d: int) -> float:
            return profile_distance(np.array(r["abund"], dtype=float), refs[(d, r["N"])])

        delta = {(d, n): 1.5 * float(np.percentile([dist(r, d) for r in cal[d] if r["N"] == n], 95)) for d in DIMS for n in NS}
        win = {d: (F3_WIN if d == 4 else {k: tuple(v) for k, v in windows(cal[d]).items()}) for d in DIMS}
        c4 = cal[4]
        tau_a = max(0.15, 3 * 1.4826 * mad([r["d_mm"] - d_spec(r["beta_L"], curve[r["N"]]) for r in c4]))
        tau_b = max(0.15, 3 * 1.4826 * mad([r["d_mm"] - r["d_mid"] for r in c4 if math.isfinite(r["d_mid"])]))
        res: dict[str, Any] = {}
        for name, (role, dtrue, _) in TEST.items():
            rr = load(f"test_{name}")
            dstar = int(min(5, max(2, round(float(np.nanmedian([r["d_mm"] for r in rr]))))))
            b = float(np.nanmedian([r["beta_L"] for r in rr]))
            ks = [r["rodilla_L"] for r in rr if r["rodilla_L"] is not None]
            km = float(np.median(ks)) if ks else None
            w = win[dstar]
            firma = bool(w["beta"][0] <= b <= w["beta"][1] and km is not None and w["rodilla"][0] <= km <= w["rodilla"][1])
            gap_a = float(np.nanmedian([abs(r["d_mm"] - d_spec(r["beta_L"], curve[r["N"]])) for r in rr]))
            ratio_c = float(np.median([dist(r, dstar) / delta[(dstar, r["N"])] for r in rr]))
            gap_b = float(np.nanmedian([abs(r["d_mm"] - r["d_mid"]) for r in rr]))
            cls = {"A": firma and gap_a <= tau_a, "C": firma and ratio_c <= 1.0, "B_info": firma and gap_b <= tau_b}
            ok = {k: (v and dstar == dtrue) if role == "positiva" else (not v) for k, v in cls.items()}
            res[name] = {"papel": role, "d_estrella": dstar, "beta_L": b, "rodilla_L": km, "firma": firma,
                         "brecha_A": gap_a, "razon_C": ratio_c, "brecha_B": gap_b, "clasificado": cls, "correcto": ok}
        verd = {arm: ("PASA" if all(v["correcto"][arm] for v in res.values()) else "FALLA") for arm in ("A", "C", "B_info")}
        out = {"commit": head, "tau_A": tau_a, "tau_B": tau_b, "curva_beta_L": curve, "veredictos": verd, "detalle": res}
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print(json.dumps({"veredictos": verd, "tau_A": round(tau_a, 3), "tau_B": round(tau_b, 3)}))
        return
    else:
        sys.exit("parte desconocida")
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "rows": rows}, default=str))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
