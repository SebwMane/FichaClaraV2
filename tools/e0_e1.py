"""Omega-E (docs/OMEGA_E_PROPUESTA.md): E0 reproduce Sorkin-Yazdi; E1 aplica la ley de area a familias no geometricas.

Partes: e0 | e1_<familia> | final. Constantes congeladas.
"""

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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from omega.causality.entropy import interval_closest, pauli_jordan, sj_both, ssee  # noqa: E402
from omega.causality.orders import kr_order, random_bipartite, thick_layer_order, transitive_closure, transitive_percolation  # noqa: E402
from omega.causality.sprinkling import causal_matrix, sprinkle_diamond  # noqa: E402

EXPERIMENT_ID = 311
ENT_E0, ENT_E1 = 20261047, 20261048
NS_E0_ORIG = (200, 400, 800, 1200, 1600)
NS = (800, 1200, 1600, 2000, 2400)  # E-A1: malla de E0b' y de E1
ENT_E0B = 20261049
SEEDS = 5
A_E0A, A_TRUNC = (0.25, 0.40), (0.23, 0.43)
MIN_REGION = 20
OUT = Path("results") / "omega_e"
Bool = npt.NDArray[np.bool_]


def rng_for(entropy: int, *key: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy, spawn_key=(EXPERIMENT_ID, *key))))


def sandwich(n: int, rng: np.random.Generator) -> Bool:
    a = (n - 2) // 2
    layer = np.concatenate([np.zeros(a, int), np.ones(2, int), np.full(n - a - 2, 2)])
    adj = (layer[:, None] + 1 == layer[None, :]) & (rng.random((n, n)) < 0.5)
    perm = rng.permutation(n)
    return transitive_closure(adj[np.ix_(perm, perm)])


FAMILIES: dict[str, Callable[[int, np.random.Generator], Bool]] = {
    "geometrico": lambda n, r: causal_matrix(sprinkle_diamond(n, 2, r)),
    "kr": lambda n, r: kr_order(n, r),
    "percolacion": lambda n, r: transitive_percolation(n, 0.05, r),
    "capas_gruesas": lambda n, r: thick_layer_order(n, 8, 0.5, r),
    "bicapa": lambda n, r: random_bipartite(n, 0.5, r),
    "sandwich": sandwich,
}


def entropies(c: Bool, region: npt.NDArray[np.int64]) -> tuple[float, float]:
    n = c.shape[0]
    idelta = pauli_jordan(c)
    w_full, idt, w_t = sj_both(idelta, math.sqrt(n) / (4 * math.pi))
    return ssee(idelta, w_full, region, truncate=False), ssee(idt, w_t, region, truncate=True)


def fits(n: list[float], s: list[float]) -> dict[str, float]:
    x = np.array(n, dtype=float)
    y = np.array(s, dtype=float)

    def r2(xx: np.ndarray[Any, Any]) -> tuple[float, float, float]:
        a = np.vstack([xx, np.ones_like(xx)]).T
        coef = np.linalg.lstsq(a, y, rcond=None)[0]
        res = y - a @ coef
        tot = float(np.sum((y - y.mean()) ** 2)) or 1e-300
        return float(coef[0]), float(coef[1]), 1 - float(res @ res) / tot

    al, bl, r2l = r2(x)
    ag, bg, r2g = r2(np.log(np.sqrt(x) / (4 * math.pi)))
    return {"a_lin": al, "b_lin": bl, "r2_lin": r2l, "a_log": ag, "b_log": bg, "r2_log": r2g}


def part_e0(ns: tuple[int, ...] = NS_E0_ORIG, entropy: int = ENT_E0) -> dict[str, Any]:
    rows = []
    for i, n in enumerate(ns):
        for s in range(SEEDS):
            x = sprinkle_diamond(n, 2, rng_for(entropy, i, s))
            u, v = x[:, 0] + x[:, 1], x[:, 0] - x[:, 1]
            region = np.flatnonzero((np.abs(u) <= 0.25) & (np.abs(v) <= 0.25)).astype(np.int64)
            s_un, s_tr = entropies(causal_matrix(x), region)
            rows.append({"N": n, "n": int(region.size), "S_untr": s_un, "S_tr": s_tr})
    return {"rows": rows}


def part_e1(fam: str, ns: tuple[int, ...]) -> dict[str, Any]:
    gen = FAMILIES[fam]
    k = list(FAMILIES).index(fam)
    rows = []
    for n in ns:
        i = NS.index(n)
        for s in range(SEEDS):
            c = gen(n, rng_for(ENT_E1, k, i, s))
            region = interval_closest(c, n // 4)
            s_un, s_tr = entropies(c, region)
            rows.append({"N": n, "n": int(region.size), "S_untr": s_un, "S_tr": s_tr})
    return {"rows": rows}


def area_flag(rows: list[dict[str, Any]]) -> tuple[bool, dict[str, float]]:
    f = fits([r["n"] for r in rows], [r["S_tr"] for r in rows])
    return bool(f["r2_log"] > f["r2_lin"] and A_TRUNC[0] <= f["a_log"] <= A_TRUNC[1]), f


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True)
    args = ap.parse_args()
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--", ".", ":!results"], capture_output=True, text=True).stdout.strip():
        sys.exit("exige arbol git limpio")
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if args.part == "final":
        e0 = json.loads((OUT / "e0.json").read_text())["data"]["rows"]
        f_un = fits([r["n"] for r in e0], [r["S_untr"] for r in e0])
        e0a = bool(A_E0A[0] <= f_un["a_lin"] <= A_E0A[1] and f_un["r2_lin"] > f_un["r2_log"])
        e0b_orig, _ = area_flag(e0)
        e0b_rows = [r for p in ("e0b_a", "e0b_b") for r in json.loads((OUT / f"{p}.json").read_text())["data"]["rows"]]
        e0b, f_tr = area_flag(e0b_rows)  # E-A1a
        summary: dict[str, Any] = {"commit": head, "E0a": e0a, "E0b_original_fallido": e0b_orig, "E0b": e0b, "fit_untr": f_un, "fit_tr": f_tr}
        if e0a and e0b:
            fam_res = {}
            for fam in FAMILIES:
                rows = [r for p in ("a", "b", "c") for r in json.loads((OUT / f"e1_{fam}_{p}.json").read_text())["data"]["rows"]]
                big = [r["n"] for r in rows if r["N"] == NS[-1]]
                if max(big) < MIN_REGION:
                    fam_res[fam] = {"ley_area": False, "nota": "no aplicable: sin regiones extensas", "n_max_1600": max(big)}
                else:
                    flag, f = area_flag(rows)
                    fam_res[fam] = {"ley_area": flag, "fit": f, "n_medio_1600": float(np.mean(big))}
            summary["E1_familias"] = fam_res
            if not fam_res["geometrico"]["ley_area"]:
                summary["E1"] = "INDETERMINADO"
            elif any(v["ley_area"] for k2, v in fam_res.items() if k2 != "geometrico"):
                summary["E1"] = "NO DISCRIMINA"
            else:
                summary["E1"] = "DISCRIMINA"
        else:
            summary["E1"] = "NO EJECUTADO (E0 fallo)"
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(json.dumps(summary, indent=2))
        return
    if args.part == "e0":
        data = part_e0()
    elif args.part == "e0b_a":
        data = part_e0(NS[:3], ENT_E0B)
    elif args.part == "e0b_b":
        data = part_e0(NS[3:], ENT_E0B)
    elif args.part.startswith("e1_") and args.part.endswith("_a"):
        data = part_e1(args.part.removeprefix("e1_").removesuffix("_a"), NS[:3])
    elif args.part.startswith("e1_") and args.part.endswith("_b"):
        data = part_e1(args.part.removeprefix("e1_").removesuffix("_b"), NS[3:4])
    elif args.part.startswith("e1_") and args.part.endswith("_c"):
        data = part_e1(args.part.removeprefix("e1_").removesuffix("_c"), NS[4:])
    else:
        sys.exit("parte desconocida")
    (OUT / f"{args.part}.json").write_text(json.dumps({"commit": head, "seconds": round(time.time() - t0, 1), "data": data}))
    print(args.part, "ok", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    main()
