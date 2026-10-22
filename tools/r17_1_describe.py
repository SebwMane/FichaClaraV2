"""R17-1, descriptivo (no oficial, sin evidencia; N <= 64): S_SJ^max frente a S_SJ de E1' (docs/OMEGA_CONSEJO_REV17_ES0.md §5.1).

Familias por N en {32, 64}: (a) 2-ordenes al azar (sprinkling 2D), 50; (b) trampa de ES-0: ordenes de altura 2
(`random_bipartite`, p = 1/2), 50; (c) KR de 3 capas (`kr_order`, p = 1/2), 50, y la cadena total (determinista, 1).
Semillas: SeedSequence(20261101, spawn_key=(326, familia, N, replica)); familia 0 = 2D, 1 = altura 2, 2 = KR.
Salida: JSON en el directorio pasado como argumento (por defecto el scratchpad de la sesion) y una tabla breve.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

from omega.causality.mcmc_sj import s_sj
from omega.causality.orders import kr_order, random_2order, random_bipartite
from omega.causality.sj_max import s_sj_max_info

ENTROPY = 20261101
EXPERIMENT_ID = 326
NS = (32, 64)
REPS = 50
FAMILIES = {"2D": 0, "altura2": 1, "KR3": 2, "cadena": 3}
DEFAULT_OUT = Path("/tmp/claude-0/-home-user-FichaClaraV2/90c7523c-b9d2-5eee-afba-260b3848fae0/scratchpad/r17_1")


def rng_for(fam: int, n: int, rep: int) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence(ENTROPY, spawn_key=(EXPERIMENT_ID, fam, n, rep))))


def make(fam: str, n: int, rng: np.random.Generator) -> npt.NDArray[np.bool_]:
    if fam == "2D":
        return random_2order(n, rng)
    if fam == "altura2":
        return random_bipartite(n, 0.5, rng)
    if fam == "KR3":
        return kr_order(n, rng, 0.5)
    return np.triu(np.ones((n, n), dtype=bool), k=1)


def summary(vals: list[float]) -> dict[str, float]:
    a = np.asarray(vals, dtype=np.float64)
    if a.size == 0:
        return {}
    return {"mean": float(a.mean()), "sd": float(a.std(ddof=1)) if a.size > 1 else 0.0, "min": float(np.min(a)),
            "median": float(np.median(a)), "max": float(np.max(a))}


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    records: list[dict[str, Any]] = []
    for fam, fid in FAMILIES.items():
        for n in NS:
            reps = 1 if fam == "cadena" else REPS
            for rep in range(reps):
                c = make(fam, n, rng_for(fid, n, rep))
                info = s_sj_max_info(c)
                s1, nr1 = s_sj(c)
                records.append({"family": fam, "N": n, "rep": rep, "spawn_key": [EXPERIMENT_ID, fid, n, rep],
                                "S_E1p": s1, "nR_E1p": nr1, **{f"max_{k}": v for k, v in info.items()}})
    groups: dict[str, Any] = {}
    rows = []
    for fam in FAMILIES:
        for n in NS:
            rs = [r for r in records if r["family"] == fam and r["N"] == n]
            both = [r for r in rs if r["nR_E1p"] > 0]
            corr = float("nan")
            if len(both) > 2:
                x = np.array([r["S_E1p"] for r in both])
                y = np.array([r["max_S"] for r in both])
                if x.std() > 0 and y.std() > 0:
                    corr = float(np.corrcoef(x, y)[0, 1])
            g = {
                "count": len(rs),
                "E1p_region_frac": len(both) / len(rs),
                "S_E1p": summary([r["S_E1p"] for r in both]),
                "nR_E1p": summary([r["nR_E1p"] for r in both]),
                "S_max": summary([r["max_S"] for r in rs]),
                "nR_max": summary([r["max_nR"] for r in rs]),
                "modes_max": summary([r["max_modes"] for r in rs]),
                "height_R": summary([r["max_height"] for r in rs]),
                "rel_frac_R": summary([r["max_rel_frac_R"] for r in rs]),
                "min1_frac": sum(r["max_min1"] for r in rs) / len(rs),
                "zero_mode_frac": sum(r["max_modes"] == 0 for r in rs) / len(rs),
                "S_max_zero_frac": sum(abs(r["max_S"]) < 1e-12 for r in rs) / len(rs),
                "pearson_S_E1p_vs_S_max": corr,
            }
            groups[f"{fam}_N{n}"] = g
            rows.append((fam, n, g))
    meta = {"official": False, "spec": "docs/OMEGA_CONSEJO_REV17_ES0.md §5.1, §6", "entropy": ENTROPY,
            "spawn_key": "(326, familia, N, replica)", "families": FAMILIES, "reps": REPS,
            "seconds": time.perf_counter() - t0}
    (out / "r17_1_describe.json").write_text(json.dumps({"meta": meta, "groups": groups, "records": records}, indent=1))
    hdr = (f"{'familia':8} {'N':>3} {'n':>3} {'E1p%':>5} {'nR_E1p':>7} {'S_E1p':>14} {'nR_max':>13} {'S_max':>14} "
           f"{'modos':>6} {'alt':>5} {'relR':>5} {'min1%':>5} {'r':>6}")
    print(hdr)
    for fam, n, g in rows:
        def ms(d: dict[str, float]) -> str:
            return f"{d['mean']:.3f}±{d['sd']:.3f}" if d else "-"
        nr = g["nR_max"]
        nr_s = f"{nr['median']:.0f} [{nr['min']:.0f},{nr['max']:.0f}]"
        nre = f"{g['nR_E1p']['median']:.0f}" if g["nR_E1p"] else "-"
        print(f"{fam:8} {n:>3} {g['count']:>3} {100 * g['E1p_region_frac']:>5.0f} {nre:>7} {ms(g['S_E1p']):>14} "
              f"{nr_s:>13} {ms(g['S_max']):>14} {g['modes_max']['median']:>6.0f} {g['height_R']['median']:>5.0f} "
              f"{g['rel_frac_R']['mean']:>5.2f} {100 * g['min1_frac']:>5.0f} {g['pearson_S_E1p_vs_S_max']:>6.2f}")
    print(f"JSON: {out / 'r17_1_describe.json'}  ({meta['seconds']:.1f} s)")


if __name__ == "__main__":
    main()
