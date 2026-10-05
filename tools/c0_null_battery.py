"""Omega-C0, control D-2 (enmienda C0-A7): D_L y D_s de los finales candidatos frente a su recableado con grados fijos."""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
from collections import defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from c0_battery import size_exponent  # noqa: E402
from c0_dynamics import _atomic_write, read_rows  # noqa: E402

from omega.c0 import references as R  # noqa: E402
from omega.c0.locality import degree_preserving_rewire, giant_component, mean_hop_connected, strong_support  # noqa: E402
from omega.config.settings import SpectralConfig  # noqa: E402
from omega.experiments.v11.gate import head_commit  # noqa: E402
from omega.geometry.spectral import spectral_dimension  # noqa: E402

OUT = ROOT / "results" / "c0_battery"


def _one(args: tuple[dict[str, Any], str]) -> dict[str, Any]:
    r, path = args
    z = np.load(path)
    w = np.asarray(z[z.files[0]], dtype=np.float64)
    a = strong_support(w)
    g = giant_component(a)
    ag = np.asarray(a[np.ix_(g, g)], dtype=np.bool_)
    rng = R.rng_from_key((20261005, 90, r["n"], r["cell_idx"], "UER".index(r["init"]), r["seed"]))
    m = int(ag.sum() // 2)
    an = degree_preserving_rewire(ag, 10 * m, rng)
    gn = giant_component(an)
    agn = np.asarray(an[np.ix_(gn, gn)], dtype=np.bool_)
    ds = spectral_dimension(agn.astype(np.float64), SpectralConfig())
    return {"cell_idx": r["cell_idx"], "init": r["init"], "n": r["n"], "seed": r["seed"],
            "mean_hop": r["mean_hop"], "D_s": r["D_s"], "null_mean_hop": mean_hop_connected(agn),
            "null_D_s": ds.value, "null_giant_frac": gn.size / g.size}


def main() -> int:
    finals = read_rows(OUT / "finals.jsonl")
    scal = json.loads((OUT / "summary.json").read_text())["size_scaling"]
    keep = {(s["cell_idx"], s["init"]) for s in scal if s["D_L"] is not None}
    tasks = []
    for r in finals:
        if (r["cell_idx"], r["init"]) not in keep:
            continue
        if r["n"] == 216:
            path = ROOT / "runs" / "c0" / "out" / f"{r['id']}.npz"
        else:
            path = ROOT / "runs" / "c0_redteam" / "out" / f"{r['id']}.npz"
        tasks.append((r, str(path)))
    with mp.get_context("fork").Pool(4) as pool:
        rows = pool.map(_one, tasks)
    by: dict[tuple[int, str], dict[str, dict[int, list[float]]]] = defaultdict(lambda: {"f": defaultdict(list), "n": defaultdict(list), "ds": defaultdict(list), "nds": defaultdict(list)})
    for x in rows:
        k = (x["cell_idx"], x["init"])
        by[k]["f"][x["n"]].append(x["mean_hop"])
        by[k]["n"][x["n"]].append(x["null_mean_hop"])
        if x["n"] == 343:
            by[k]["ds"][343].append(x["D_s"])
            by[k]["nds"][343].append(x["null_D_s"])
    table = []
    for (c, i), d in sorted(by.items()):
        dl, dln = size_exponent(d["f"]), size_exponent(d["n"])
        def med(v: list[float]) -> float | None:
            u = [t for t in v if t is not None and np.isfinite(t)]
            return float(np.median(u)) if u else None
        sep = dl is not None and dln is not None and np.isfinite(dln) and dl < dln - 0.5
        table.append({"cell_idx": c, "init": i, "D_L": dl, "D_L_null": dln, "D_s_343": med(d["ds"][343]),
                      "D_s_343_null": med(d["nds"][343]), "separated": bool(sep)})
    out = {"step": "c0_null_battery", "code_commit": head_commit(), "table": table,
           "n_separated": sum(t["separated"] for t in table), "n_pairs": len(table)}
    _atomic_write(OUT / "null_control.json", json.dumps(out, indent=2, default=str))
    for t in table:
        print(t["cell_idx"], t["init"], {k: (round(v, 2) if isinstance(v, float) else v) for k, v in t.items() if k not in ("cell_idx", "init")})
    print("separados:", out["n_separated"], "/", out["n_pairs"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
