"""F1 y F2 (docs/OMEGA_F1F2_PREREGISTRO.md). Constantes congeladas. Uso: --part f1 | f2_<d> | final"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import numpy.typing as npt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.entropy import knee_from_spectrum, pauli_jordan  # noqa: E402
from omega.causality.orders import layered_2order, transitive_percolation  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 314
ENT_F1, ENT_F2 = 20261058, 20261059
BETA_WIN, KNEE_WIN = (-1.104, -0.904), (0.227, 0.906)
OUT = Path("results") / "omega_f1f2"
Bool = npt.NDArray[np.bool_]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def de_sitter_2d(n: int, r: np.random.Generator) -> Bool:
    pts: list[tuple[float, float]] = []
    while len(pts) < n:
        eta = r.uniform(-1.5, -0.5)
        x = r.uniform(-0.5, 0.5)
        if abs(x) <= 0.5 - abs(eta + 1.0) and r.random() < (1.0 / eta**2) / 4.0:
            pts.append((eta, x))
    return causal_matrix(np.array(pts))


def cylinder_2d(n: int, r: np.random.Generator) -> Bool:
    t = r.uniform(0, 0.4, n)
    x = r.uniform(0, 1, n)
    dt = t[None, :] - t[:, None]
    dx = (x[None, :] - x[:, None] + 0.5) % 1.0 - 0.5
    rel: Bool = (dt > 0) & (dt >= np.abs(dx))
    return rel


def random_dorder(d: int) -> Callable[[int, np.random.Generator], Bool]:
    def gen(n: int, r: np.random.Generator) -> Bool:
        u = r.random((n, d))
        rel: Bool = np.all(u[:, None, :] < u[None, :, :], axis=2)
        return rel
    return gen


F1_STRUCTS: dict[str, tuple[str, Callable[[int, np.random.Generator], Bool]]] = {
    "de_sitter_2d": ("positiva", de_sitter_2d),
    "cilindro_2d": ("positiva", cylinder_2d),
    "capas_200": ("negativa", lambda n, r: layered_2order(n, 200, r)),
    "capas_400": ("negativa", lambda n, r: layered_2order(n, 400, r)),
    "orden_d3": ("negativa", random_dorder(3)),
    "orden_d4": ("negativa", random_dorder(4)),
    "percolacion_002": ("negativa", lambda n, r: transitive_percolation(n, 0.02, r)),
    "percolacion_010": ("negativa", lambda n, r: transitive_percolation(n, 0.10, r)),
}


def signature(c: Bool) -> tuple[float | None, float]:
    mu = np.asarray(np.linalg.eigvalsh(pauli_jordan(c)), dtype=np.float64)
    return knee_from_spectrum(mu, c.shape[0])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.part == "f1":
        rows = []
        for k, (name, (role, gen)) in enumerate(F1_STRUCTS.items()):
            for n in (800, 1600):
                for s in range(5):
                    kn, beta = signature(gen(n, rng_for(ENT_F1, k, n, s)))
                    rows.append({"estructura": name, "papel": role, "N": n, "seed": s, "rodilla": kn, "beta": beta})
        (OUT / "f1.json").write_text(json.dumps({"commit": head, "rows": rows}, default=str))
    elif args.part.startswith("f2_"):
        d = int(args.part[3:])
        rows = []
        for n in (400, 800, 1600, 2400):
            for s in range(5):
                kn, beta = signature(causal_matrix(sprinkle_diamond(n, d, rng_for(ENT_F2, d, n, s))))
                rows.append({"d": d, "N": n, "seed": s, "rodilla": kn, "beta": beta})
        (OUT / f"f2_{d}.json").write_text(json.dumps({"commit": head, "rows": rows}, default=str))
    elif args.part == "final":
        f1 = json.loads((OUT / "f1.json").read_text())["rows"]
        res1: dict[str, Any] = {}
        for name, (role, _) in F1_STRUCTS.items():
            rr = [r for r in f1 if r["estructura"] == name]
            b = float(np.nanmedian([r["beta"] for r in rr]))
            ks = [r["rodilla"] for r in rr if r["rodilla"] is not None]
            kmed = float(np.median(ks)) if ks else None
            comp = bool(BETA_WIN[0] <= b <= BETA_WIN[1] and kmed is not None and KNEE_WIN[0] <= kmed <= KNEE_WIN[1])
            res1[name] = {"papel": role, "beta_mediana": b, "rodilla_mediana": kmed, "con_rodilla": len(ks), "compatible": comp,
                          "correcto": comp if role == "positiva" else not comp}
        f1v = "PASA" if all(v["correcto"] for v in res1.values()) else "FALLA"
        res2: dict[str, Any] = {}
        for d in (2, 3, 4):
            rows2 = json.loads((OUT / f"f2_{d}.json").read_text())["rows"]
            med = {n: float(np.nanmedian([r["beta"] for r in rows2 if r["N"] == n])) for n in (400, 800, 1600, 2400)}
            tgt = -2 / d
            ok = abs(med[2400] - tgt) < abs(med[400] - tgt) and abs(med[2400] - tgt) <= 0.10
            res2[str(d)] = {"objetivo": tgt, "beta_por_N": med, "pasa": bool(ok),
                            "rodilla_2400": float(np.median([r["rodilla"] for r in rows2 if r["N"] == 2400 and r["rodilla"] is not None] or [float("nan")]))}
        f2v = "PASA" if all(v["pasa"] for v in res2.values()) else "FALLA"
        out = {"commit": head, "F1": f1v, "F1_detalle": res1, "F2": f2v, "F2_detalle": res2}
        (OUT / "summary.json").write_text(json.dumps(out, indent=2, default=str))
        print(json.dumps({"F1": f1v, "F2": f2v}, indent=2))
        return
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
